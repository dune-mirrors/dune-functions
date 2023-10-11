// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LINEARTRANSFORMEDLOCALFINITEELEMENT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LINEARTRANSFORMEDLOCALFINITEELEMENT_HH

#include <array>
#include <dune/common/std/type_traits.hh>
#include <dune/functions/functionspacebases/globalvaluedlocalfiniteelement.hh>
#include <dune/functions/functionspacebases/transformatorconcepts.hh>
#include <dune/istl/scaledidmatrix.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <type_traits>
#include <vector>

namespace Dune::Functions
{
/**
 * \brief This set of classes models a global valued Finite Element similar to
 * the classes in
 * dune/functions/functionspacebases/globalvaluedlocalfiniteelement.hh
 *
 * Focus is on linear transformations of the basisfunctions, as is the case for
 * various Elements like Hermite, Morley and Argyris. The Transformation is
 * objectified to implemenent caching of the transformation. This leads to some
 * changes in the interface in contrast to range-space transformations like
 * Piola-Transformations.
 * The main contrast is, that the transformations here take the concrete
 * structure of finite elements in account and are thus specific to each finite
 * element, unlike the Piola- Transformations, which can be applied to a set of
 * admissible finite elements. Generally, we cannot use the reference dofs for
 * interpolation, so two different approaches are offered.
 *
 * In particular, the interfaces differ from the "GlobalValued"-classes in the
 * following points: LinearTransformator interface:
 *      - The LinearTransformator classes implement a linear transformation from
 * one set of basisfunction onto another, that means, the transformation is
 * invariant to differentiation, and the need for methods like "applyJacobian"
 * vanishes.
 *      - The LinearTransformator classes are not static.
 *      - The LinearTransformator object are included in binding routine.
 *      - Interpolation: Instead of the inner class LocalValuedFunction, two
 * options are possible:
 *        1. affine (Piola) interpolation equivalent FEs:
 *          The LinearTransformator offers a method applyInverse, which is then
 * wrapped to a LocalValuedFunction.
 *        2. non intepolation equivalent FEs:
 *          We cannot use the reference DOFs to interpolate into the physical FE
 * space. Hence the LinearTransformators need to export an inner class
 *          GlobalValuedInterpolation, which fulfils the LocalInterpolation
 * interface and fill the coefficient vector for the transformed FE.
 *      - Exports type ElementInformation, that stores elementspecific data,
 * like orientation or direction of vertices. Should implement the methods
 * isDirichlet(LocalKey)/ isClamped(LocalKey) which return a boolean, which true
 * iff the localKey is intended to be used to interpolate Dirichlet/Clamped BC.
 *    LinearTransformedLocalBasis:
 *      - holds a pointer to the reference basis
 *      - higher derivatives are implemented, possible throwing errors if
 * reference basis does not implement at least the partial method of
 * corresponding degree
 *
 *    LinearTransformedLocalInterpolation:
 *      - provided by Transformator class as exported type GlobalInterpolation
 * or via CoefficientTransformingGlobalValuedInterpolation
 *
 *    LinearTransformedLocalCoefficients:
 *      - wraps the LocalCoefficients of the localvalued finite element
 *      - can be bound to ElementInformation
 *      - provides an additional method isDirichlet(size_t i)/isClamped(size_t
 * i) which return boolean indicating whether the ith local DOF is to be used in
 * Dirichlet Interpolation
 *
 *    LinearTransformedLocalFiniteElement:
 *      - holds an instance of the Transformator class
 *      - calls transformator.bind(...) whenever bound to an element
 *      - similar to GlobalValuedLocalFiniteElement with adapted typedefs and
 * additional binding operations
 *
 */

namespace Impl
{
// forward declaration
template<class Transformator, class LocalValuedLFE, class Element>
class TransformedLocalFiniteElement;

/** \brief Implementation of a dune-localfunctions LocalBasis that applies a
 * linear transformation
 *
 * \tparam Transformator The transformation that is to be applied
 * \tparam LocalValuedLocalBasis The local-valued LocalBasis that is getting
 * transformed \tparam Element The element that the global-valued FE lives on
 */
template<class Transformator, class LocalValuedLFE, class Element>
class TransformedLocalBasis
{
    friend class TransformedLocalFiniteElement<Transformator, LocalValuedLFE, Element>;
    using LocalValuedLocalBasis = typename LocalValuedLFE::Traits::LocalBasisType;

  public:
    using Traits = typename LocalValuedLocalBasis::Traits;
    using LocalValuedRangeType = typename Traits::RangeType;
    using LocalValuedJacobianType = typename Traits::JacobianType;
    using LocalValuedHessianType = typename Impl::HessianType<LocalValuedLocalBasis>::type;

    TransformedLocalBasis(Transformator const &transformator) : transformator_(&transformator) {}

  private:
    /** Bind the Basis to an element
     */
    void bind(Element const &, LocalValuedLocalBasis const &localValuedLocalBasis)
    {
      localValuedLocalBasis_ = &localValuedLocalBasis;
    }

  public:
    /** \brief Number of shape functions
        This does generally not equal the size of the underlying local finite
       element.
     */
    auto size() const
    {
      if constexpr (models<Concept::SizeProvidingTransformator<Element, LocalValuedLFE>,
                           Transformator>())
        return transformator_->size();
      else
        return localValuedLocalBasis_->size();
    }

    //! \brief Evaluate all shape functions
    void evaluateFunction(const typename Traits::DomainType &x,
                          std::vector<typename Traits::RangeType> &out) const
    {

      thread_local std::vector<LocalValuedRangeType> rangeBuffer_;
      rangeBuffer_.resize(localValuedLocalBasis_->size());
      localValuedLocalBasis_->evaluateFunction(x, rangeBuffer_);
      out.resize(size());
      transformator_->transform(rangeBuffer_, out);
    }

    /** \brief Evaluate Jacobian of all shape functions
     *
     * \param x Point in the reference element where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
    void evaluateJacobian(const typename Traits::DomainType &x,
                          std::vector<typename Traits::JacobianType> &out) const
    {

      thread_local std::vector<LocalValuedJacobianType> jacobianBuffer_;
      jacobianBuffer_.resize(localValuedLocalBasis_->size());
      localValuedLocalBasis_->evaluateJacobian(x, jacobianBuffer_);
      out.resize(size());
      transformator_->transform(jacobianBuffer_, out);
    }

    /** \brief Evaluate Hessian of all shape functions
     *   \note Sfinae protected: this method is only available, if the wrapped
     * LocalBasis 1. exports a <HessianType> and 2. provides a evaluateHessian
     * method with a corresponding signature \param x Point in the reference
     * element where to evaluation the Hessians \param[out] out The Hessians of
     * all shape functions at the point x
     */
    template<class TT = LocalValuedLocalBasis,
             std::enable_if_t<models<Concept::H2Basis, TT>(), int> = 0>
    void evaluateHessian(const typename Traits::DomainType &x,
                         std::vector<typename TT::HessianType> &out) const
    {

      thread_local std::vector<typename TT::HessianType> hessianBuffer_;
      static_assert(std::is_same_v<TT, LocalValuedLocalBasis>);
      hessianBuffer_.resize(localValuedLocalBasis_->size());
      localValuedLocalBasis_->evaluateHessian(x, hessianBuffer_);
      out.resize(size());
      transformator_->transform(hessianBuffer_, out);
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     *
     * \param order Order of the partial derivatives, in the classic multi-index
     * notation \param in Position where to evaluate the derivatives \param[out]
     * out The desired partial derivatives
     */
    void partial(std::array<unsigned int, Traits::dimDomain> const &order,
                 const typename Traits::DomainType &x,
                 std::vector<typename Traits::RangeType> &out) const
    {

      thread_local std::vector<LocalValuedRangeType> rangeBuffer_;
      rangeBuffer_.resize(localValuedLocalBasis_->size());
      localValuedLocalBasis_->partial(order, x, rangeBuffer_);
      out.resize(size());
      transformator_->transform(rangeBuffer_, out);
    }

    //! \brief Polynomial order of the shape functions
    auto order() const { return localValuedLocalBasis_->order(); }

    Transformator const &transformator() const { return *transformator_; }

    LocalValuedLocalBasis const &cacheable() const { return *localValuedLocalBasis_; }

  private:
    Transformator const *transformator_;
    LocalValuedLocalBasis const *localValuedLocalBasis_;
};

/** \brief Implementation of a dune-localfunctions LocalInterpolation
 *    that accepts global-valued functions
 *
 * \tparam Transformator The transformation (e.g., Piola) that transforms from local to global
 * values \tparam LocalValuedLocalInterpolation The local-valued LocalInterpolation that is used for
 * the actual interpolation \tparam Element The element that the global-valued FE lives on
 */
template<class Transformator, class LocalValuedLocalInterpolation, class Element>
class RangeSpaceTransformingGlobalValuedLocalInterpolation
{
  public:
    /** \brief Bind the local interpolation object to a particular grid element
     */
    void bind(Element const &element,
              LocalValuedLocalInterpolation const &localValuedLocalInterpolation)
    {
      localValuedLocalInterpolation_ = &localValuedLocalInterpolation;
      element_ = &element;
    }

    template<typename F, typename C>
    void interpolate(const F &f, std::vector<C> &out) const
    {
      using LocalCoordinate = typename Element::Geometry::LocalCoordinate;
      typename Transformator::template LocalValuedFunction<F, LocalCoordinate, Element>
          localValuedFunction(f, *element_);
      localValuedLocalInterpolation_->interpolate(localValuedFunction, out);
    }

  private:
    LocalValuedLocalInterpolation const *localValuedLocalInterpolation_;
    Element const *element_;
};

/**
 * @brief This class is a generic implementation of the GlobalValuedInterpolation
 * interface. After binding, it uses an invertible Basis Transformator and
 * the reference LocalInterpolation and to construct a GlobalValuedLocalInterpolation.
 *
 * @tparam Transformator
 * @tparam LocalInterpolation
 * @tparam Element
 */
template<class Transformator, class LocalFE, class Element>
class CoefficientTransformingGlobalValuedLocalInterpolation
{
    friend class TransformedLocalFiniteElement<Transformator, LocalFE, Element>;
    using LocalInterpolation = typename LocalFE::Traits::LocalInterpolationType;

  public:
    CoefficientTransformingGlobalValuedLocalInterpolation(Transformator const &t)
        : transformator_(&t)
    {
    }

  private:
    void bind(Element const &e, LocalInterpolation const &lI)
    {
      localInterpolation_ = &lI;
      element_ = &e;
    }

  public:
    template<class F, class C>
    void interpolate(const F &f, std::vector<C> &out) const
    {
      localInterpolation_->interpolate(f, out);
      transformator_->applyInverse(out);
    }

  private:
    Transformator const *transformator_;
    LocalInterpolation const *localInterpolation_;
    Element const *element_;
};

// The transformator class either exports a LocalValuedLocalFunction, a
// GlobalValuedInterpolation, or implements an inverse transformation Primary
// Template
template<class Transformator, class LocalFE, class Element, class Enabled = bool>
struct GlobalValuedInterpolationType;

// Spezialization for Transformator class that implements a "applyInverse"
// method. This should only be used for Finite elements which are interpolation
// affine equivalent.
template<class Transformator, class LocalFE, class Element>
struct GlobalValuedInterpolationType<
    Transformator, LocalFE, Element,
    std::enable_if_t<
        models<Concept::InterpolationEquivalentTransformator<Element, LocalFE>, Transformator>(),
        bool>> {
    using type =
        CoefficientTransformingGlobalValuedLocalInterpolation<Transformator, LocalFE, Element>;
};

// Specialization for Transformator classes which implement a
// GlobalValuedInterpolation class. This approach is suitable for non (!)
// interpolation affine equivalent FEs.
template<class Transformator, class LocalFE, class Element>
struct GlobalValuedInterpolationType<
    Transformator, LocalFE, Element,
    std::enable_if_t<
        models<Concept::InterpolationProvidingTransformator<Element, LocalFE>, Transformator>(),
        bool>> {
    using type = typename Transformator::GlobalValuedInterpolation;
};

// Specialization for Transformator classes which implement a
// LocalValuedFunction class. This approach is suitable for FEs with a
// Rangespace transforming pullback like Piola transformations.
template<class Transformator, class LocalFE, class Element>
struct GlobalValuedInterpolationType<
    Transformator, LocalFE, Element,
    std::enable_if_t<models<Concept::RangeSpaceTransformator<Element, LocalFE>, Transformator>(),
                     bool>> {
    using type = RangeSpaceTransformingGlobalValuedLocalInterpolation<
        Transformator, typename LocalFE::Traits::LocalInterpolation, Element>;
};

/** \brief Associations of the transformed degrees of freedom to subentities of
 * the reference simplex. \note We assume here, that the transformation does not
 * change the LocalCoefficients, i.e. it does not change to which subentity a
 * DoF is associated with. However, it may reduce the size (e.g. for the reduced
 * Hermite Element / DKT).
 *
 * \tparam dim Dimension of the reference simplex
 */
// TODO discuss whether one can include the transformedindexbasis feature in
// here as well
template<class Transformator, class LocalValuedLFE, class Element>
class TransformedLocalCoefficients
{
    friend class TransformedLocalFiniteElement<Transformator, LocalValuedLFE, Element>;
    using LocalValuedLocalCoefficients = typename LocalValuedLFE::Traits::LocalCoefficientsType;
    using LocalValuedLocalBasis = typename LocalValuedLFE::Traits::LocalBasisType;

  public:
    using size_type = typename LocalValuedLocalCoefficients::size_type;

  public:
    TransformedLocalCoefficients(Transformator const &transformator)
        : transformator_(&transformator)
    {
    }

  private:
    /** Bind the Coefficients to an element
     */
    void bind(Element const &, LocalValuedLocalCoefficients const &localCoefficients)
    {
      localValuedLocalCoefficients_ = &localCoefficients;
    }

  public:
    //! number of coefficients
    size_type size() const
    {
      if constexpr (models<Concept::SizeProvidingTransformator<Element, LocalValuedLFE>,
                           Transformator>())
        return transformator_->size();
      else
        return localValuedLocalCoefficients_->size();
    }

    //! get i'th index
    LocalKey const &localKey(size_type i) const
    {
      return localValuedLocalCoefficients_->localKey(i);
    }

  private:
    Transformator const *transformator_;
    LocalValuedLocalCoefficients const *localValuedLocalCoefficients_;
};

/** \brief LocalFiniteElement implementation that uses values defined wrt
 * particular grid elements
 *
 * \tparam Transformator Class implementing linear transformations
 * \tparam LocalValuedLFE LocalFiniteElement implementation whose values are to
 * be transformed
 * \tparam Element Element onto which the FE is transformed
 */
template<class Transformator, class LocalValuedLFE, class Element>
class TransformedLocalFiniteElement
{

    using LocalBasis = TransformedLocalBasis<Transformator, LocalValuedLFE, Element>;
    using LocalInterpolation =
        typename GlobalValuedInterpolationType<Transformator, LocalValuedLFE, Element>::type;
    using LocalCoefficients = TransformedLocalCoefficients<Transformator, LocalValuedLFE, Element>;

  public:
    /** \brief Export number types, dimensions, etc.
     */
    using Traits = LocalFiniteElementTraits<LocalBasis, LocalCoefficients, LocalInterpolation>;

    // Variadic Args, which are forwarded to the constructur of the Generic
    // Reference Element.
    template<class GlobalState, class... Args>
    TransformedLocalFiniteElement(GlobalState &&globalState, Args &&...args)
        : transformator_(std::forward<GlobalState>(globalState)),
          localValuedLFE_(std::forward<Args>(args)...), globalValuedLocalBasis_(transformator_),
          globalValuedLocalInterpolation_(transformator_),
          globalValuedLocalCoefficients_(transformator_)

    {
    }
    /**
     * \brief Binding routine. This binds all other objects of the
     * LocalFiniteElement interface, as well as the Transformator.
     *
     * \param localValuedLFE
     * \param element
     * \param elementInfo
     */

    void bind(Element const &element)
    {
      element_ = &element;
      // bind the transformator, i.e. also the global State
      transformator_.bind(*element_);
      // To allow chaining, we bind the localValued fe as well
      if constexpr (models<Concept::BindableTo<Element>, LocalValuedLFE>())
        localValuedLFE_.bind(*element_);
      globalValuedLocalBasis_.bind(*element_, localValuedLFE_.localBasis());
      globalValuedLocalCoefficients_.bind(*element_, localValuedLFE_.localCoefficients());
      globalValuedLocalInterpolation_.bind(*element_, localValuedLFE_.localInterpolation());
    }

    /** \brief Returns the local basis, i.e., the set of shape functions
     */
    const typename Traits::LocalBasisType &localBasis() const { return globalValuedLocalBasis_; }

    /** \brief Returns the assignment of the degrees of freedom to the element
     * subentities
     */
    const typename Traits::LocalCoefficientsType &localCoefficients() const
    {
      return globalValuedLocalCoefficients_;
    }

    /** \brief Returns object that evaluates degrees of freedom
     */
    const typename Traits::LocalInterpolationType &localInterpolation() const
    {
      return globalValuedLocalInterpolation_;
    }

    /** \brief The number of shape functions */
    std::size_t size() const
    {
      if constexpr (models<Concept::SizeProvidingTransformator<Element, LocalValuedLFE>,
                           Transformator>())
        return transformator_.size();
      else
        return localValuedLFE_.size();
    }

    /** \brief The reference element that the local finite element is defined on
     */
    GeometryType type() const { return localValuedLFE_.type(); }

  private:
    Transformator transformator_;
    LocalValuedLFE localValuedLFE_;
    typename Traits::LocalBasisType globalValuedLocalBasis_;
    typename Traits::LocalInterpolationType globalValuedLocalInterpolation_;
    typename Traits::LocalCoefficientsType globalValuedLocalCoefficients_;
    Element const *element_;
};
} // namespace Impl

template<typename GV, class Transformator, typename LocalValuedLocalFiniteElement>
class TransformedNode : public LeafBasisNode
{

    static constexpr unsigned int dim = GV::dimension;

  public:
    using size_type = std::size_t;
    using Element = typename GV::template Codim<0>::Entity;

    using FiniteElement =
        Impl::TransformedLocalFiniteElement<Transformator, LocalValuedLocalFiniteElement, Element>;

    template<class GlobalState>
    TransformedNode(GlobalState const &globalState) : finiteElement_(globalState)
    {
      this->setSize(finiteElement_.size());
    }

    ~TransformedNode() {}

    //! Return current element, throw if unbound
    Element const &element() const { return element_; }

    /** \brief Return the LocalFiniteElement for the element we are bound to
     *
     * The LocalFiniteElement implements the corresponding interfaces of the
     * dune-localfunctions module
     */
    FiniteElement const &finiteElement() const { return finiteElement_; }

    //! Bind to element.
    void bind(Element const &e)
    {
      element_ = e;
      finiteElement_.bind(element_);
      this->setSize(finiteElement_.size());
    }

    unsigned int order() const { return finiteElement_.localBasis().order(); }

  protected:
    FiniteElement finiteElement_;
    Element element_;
};

} // namespace Dune::Functions
#endif
