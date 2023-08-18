// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LINEARTRANSFORMEDLOCALFINITEELEMENT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LINEARTRANSFORMEDLOCALFINITEELEMENT_HH
#include "dune/localfunctions/common/localkey.hh"
#include <dune/common/std/type_traits.hh>
#include <dune/functions/functionspacebases/globalvaluedlocalfiniteelement.hh>
#include <dune/istl/scaledidmatrix.hh>
#include <type_traits>
#include <dune/functions/functionspacebases/transformatorconcepts.hh>

namespace Dune::Functions {
namespace Impl {
/**
 * \brief This set of classes models a global valued Finite Element similar to the classes in
 * dune/functions/functionspacebases/globalvaluedlocalfiniteelement.hh
 *
 * Focus is on linear transformations of the basisfunctions, as is the case for various Elements
 * like Hermite, Morley and Argyris.
 * The Transformation is objectified to implemenent caching of the transformation.
 * This leads to some changes in the interface in contrast to range-space transformations like
 * Piola-Transformations.
 * The main contrast is, that the transformations here take the concrete structure of finite
 * elements in account and are thus specific to each finite element, unlike the Piola-
 * Transformations, which can be applied to a set of admissible finite elements.
 * Generally, we cannot use the reference dofs for interpolation, so two different approaches are
 * offered.
 *
 * In particular, the interfaces differ from the "GlobalValued"-classes in the following points:
 *    LinearTransformator interface:
 *      - The LinearTransformator classes implement a linear transformation from one set of
 *        basisfunction onto another, that means, the transformation is invariant to
 *        differentiation, and the need for methods like "applyJacobian" vanishes.
 *      - The LinearTransformator classes are not static.
 *      - The LinearTransformator object are included in binding routine.
 *      - Interpolation: Instead of the inner class LocalValuedFunction, two options are possible:
 *        1. affine (Piola) interpolation equivalent FEs:
 *          The LinearTransformator offers a method applyInverse, which is then wrapped to a
 *          LocalValuedFunction.
 *        2. non intepolation equivalent FEs:
 *          We cannot use the reference DOFs to interpolate into the physical FE space.
 *          Hence the LinearTransformators need to export an inner class
 *          GlobalValuedInterpolation, which fulfils the LocalInterpolation interface and fill
 *          the coefficient vector for the transformed FE.
 *      - Exports type ElementInformation, that stores elementspecific data, like orientation
 *        or direction of vertices. Should implement the methods isDirichlet(LocalKey)/
 *        isClamped(LocalKey) which return a boolean, which true iff the localKey is intended
 *        to be used to interpolate Dirichlet/Clamped BC.
 *    LinearTransformedLocalBasis:
 *      - holds a pointer to the reference basis
 *      - higher derivatives are implemented, possible throwing errors if reference basis does
 *        not implement at least the partial method of corresponding degree
 *
 *    LinearTransformedLocalInterpolation:
 *      - provided by Transformator class as exported type GlobalInterpolation or via
 *        DefaultGlobalValuedInterpolation
 *
 *    LinearTransformedLocalCoefficients:
 *      - wraps the LocalCoefficients of the localvalued finite element
 *      - can be bound to ElementInformation
 *      - provides an additional method isDirichlet(size_t i)/isClamped(size_t i) which return
 *        boolean indicating whether the ith local DOF is to be used in Dirichlet Interpolation
 *
 *    LinearTransformedLocalFiniteElement:
 *      - holds an instance of the Transformator class
 *      - calls transformator.bind(...) whenever bound to an element
 *      - similar to GlobalValuedLocalFiniteElement with adapted typedefs and additional
 *        binding operations
 *
 */
/** \brief Example of LinearTransformation Interface
 *  \tparam F Field type
 *  \tparam size number of basisfunctions
 */

template <class F, unsigned int size>
struct IdentityTransformation {

    template <class Element>
    class ElementInformation {};

    IdentityTransformation() : m_(1) {}

    template <class LocalBasis, class Element>
    void bind(LocalBasis const &lB, Element e, ElementInformation<Element> eInfo) {
      // Adapt values
    }

    template <typename Values>
    void apply(Values &values) {
      // Values tmp = values;
      // m_.mv(tmp, values);
    }

    template <typename Values>
    void applyInverse(Values &values) {
      // Values tmp = values;
      // m_.mv(tmp, values);
    }

    /**Alternative to applyInverse(values)
    This approach can be used for finite elements, for which the linear transformed FE is not
    interpolation equivalent to the generic FE*/
    // template <class LocalBasis, class Element>
    // class GlobalInterpolation
    // {
    //   using size_type = std::size_t;
    //   using LocalCoordinate = typename LocalBasis::Traits::DomainType;

    // public:
    //   GlobalInterpolation() {}

    //   /**
    //    * \brief Binding routine. Collects tangentials and normals
    //    * \param element
    //    * \param elementInfo
    //    */
    //   void bind(const Element &element, const ElementInformation<Element> &elementInfo)
    //   {
    //     // set up interpolation
    //   }

    //   /** \brief Evaluate a given function at the Lagrange nodes
    //    *
    //    * \tparam F Type of function to evaluate
    //    * \tparam C Type used for the values of the function
    //    * \param[in] ff Function to evaluate
    //    * \param[out] out Array of function values
    //    */
    //   template <typename Func, typename C>
    //   void interpolate(const Func &ff, std::vector<C> &out) const
    //   {
    //     // This is specific to the finite element to be wrapped
    //   }
    // };

    ScaledIdentityMatrix<F, size> m_;
};
} // namespace Impl

namespace Impl {

/** \brief Implementation of a dune-localfunctions LocalBasis that applies a linear
 * transformation
 *
 * \tparam Transformator The transformation that is to be applied
 * \tparam LocalValuedLocalBasis The local-valued LocalBasis that is getting transformed
 * \tparam Element The element that the global-valued FE lives on
 */
template <class LinearTransformator, class LocalValuedLocalBasis, class Element>
class LinearTransformedLocalBasis {
    using ElementInformation = typename LinearTransformator::template ElementInformation<Element>;

  public:
    using Traits = typename LocalValuedLocalBasis::Traits;

    /** \brief Bind the local basis to a particular grid element
     */
    void bind(const LocalValuedLocalBasis &localValuedLocalBasis, LinearTransformator const &t,
              const Element &element, const ElementInformation &elementInfo) {
      transformator_ = &t;
      localValuedLocalBasis_ = &localValuedLocalBasis;
      element_ = element;
      elementInformation_ = &elementInfo;
    }

    /** \brief Number of shape functions
     */
    auto size() const { return localValuedLocalBasis_->size(); }

    //! \brief Evaluate all shape functions
    void evaluateFunction(const typename Traits::DomainType &x,
                          std::vector<typename Traits::RangeType> &out) const {
      localValuedLocalBasis_->evaluateFunction(x, out);

      transformator_->apply(out);
    }

    /** \brief Evaluate Jacobian of all shape functions
     *
     * \param x Point in the reference element where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
    void evaluateJacobian(const typename Traits::DomainType &x,
                          std::vector<typename Traits::JacobianType> &out) const {
      localValuedLocalBasis_->evaluateJacobian(x, out);

      transformator_->apply(out);
    }

    /** \brief Evaluate Hessian of all shape functions
     *   \note Sfinae protected: this method is only available, if the wrapped LocalBasis 1.
     * exports a <HessianType> and 2. provides a evaluateHessian method with a corresponding
     * signature
     * \param x Point in the reference element where to evaluation the Hessians
     * \param[out] out The Hessians of all shape functions at the point x
     */
    template <class TT = LocalValuedLocalBasis,
              std::enable_if_t<models<Concept::H2Basis, TT>(), int> = 0>
    void evaluateHessian(const typename Traits::DomainType &x,
                         std::vector<typename TT::HessianType> &out) const {

      static_assert(std::is_same_v<TT, LocalValuedLocalBasis>);
      localValuedLocalBasis_->evaluateHessian(x, out);
      transformator_->apply(out);
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     *
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out The desired partial derivatives
     */
    void partial(const std::array<unsigned int, Traits::dimDomain> &order,
                 const typename Traits::DomainType &x,
                 std::vector<typename Traits::RangeType> &out) const {
      localValuedLocalBasis_->partial(order, x, out);
      transformator_->apply(out);
    }

    //! \brief Polynomial order of the shape functions
    auto order() const { return localValuedLocalBasis_->order(); }

    LinearTransformator const &transformator() { return *transformator_; }

  private:
    const LocalValuedLocalBasis *localValuedLocalBasis_;
    Element element_;
    const ElementInformation *elementInformation_;
    const LinearTransformator *transformator_;
};

/**
 * \brief Class representing Linear transformed coefficients. This class mainly forwards the
 * LocalCoefficient class of the reference finite element.
 * We assume here that the transformation does not change the location of a DOF.
 * This class extends the LocalCoefficient interface by providing a binding methods,
 * such that it can access the ElementInformation object, and provides methods to determine
 * whether the ith degree of freedom is to be used when strongly incorporation Dirichlet or
 * Clamped boundary conditions. This information is obtained from the ElementInformation object.
 * \tparam LocalCoefficients
 * \tparam ElementInformation
 */
template <class LocalCoefficients, class ElementInformation>
class LinearTransformedLocalCoefficients : public LocalCoefficients {

  public:
    // TODO maybe variadic forward constructor
    LinearTransformedLocalCoefficients() : LocalCoefficients(), elementInfo_(nullptr) {}
    void bind(ElementInformation const &elementInfo) { elementInfo_ = &elementInfo; }

    /**
     * \brief Whether or not to use the ith dof when strongly incorporating Dirichlet Conditions.
     * If the ElementInformation class does not implement a corresponding methods, this method
     * will always return true. Also note that this methods does not check whether this is actually a boundary dof.
     *
     * \param i
     * \return true
     * \return false
     */
    bool isDirichlet(std::size_t i) const {
      return isDirichletImpl(*elementInfo_, i, PriorityTag<42>{});
    }

    /**
     * \brief Whether or not to use the ith dof when strongly incorporating Clamped Conditions. If
     * the ElementInformation class does not implement a corresponding methods, this method will
     * always return true. Also note that this methods does not check whether this is actually a
     * boundary dof.
     *
     * \param i
     * \return true
     * \return false
     */
    bool isClamped(std::size_t i) const {
      return isClampedImpl(*elementInfo_, i, PriorityTag<42>{});
    }

    ElementInformation const &getElementInformation() const { return *elementInfo_; }

  private:
    bool isDirichletImpl(ElementInformation const &elInfo, std::size_t i,
                         PriorityTag<0> tag) const {
      return true;
    }

    // template ElInfo to enable sfinae
    template <class ElInfo>
    bool isDirichletImpl(
        ElInfo const &elInfo, std::size_t i, PriorityTag<1> tag,
        decltype((std::declval<ElInfo>().isDirichlet(std::declval<LocalKey>()),
                  true)) = true) const {
      return elInfo.isDirichlet(LocalCoefficients::localKey(i));
    }

    bool isClampedImpl(ElementInformation const &e, std::size_t i, PriorityTag<0> tag) const {
      return true;
    }
    template <class ElInfo, decltype((std::declval<ElInfo>().isClamped(
                                          std::declval<LocalKey>()),
                                      true)) = true>
    bool isClampedImpl(ElInfo const &elInfo, std::size_t i, PriorityTag<1> tag) const {
      return elInfo.isClamped(LocalCoefficients::localKey(i));
    }

    ElementInformation const *elementInfo_;
};

/**
 * @brief This class is a generic implementation of the LocalValuedFunction interface. Upon binding, it wraps a InterpolationEquivalentTransformator and uses the reference LocalInterpolation and the transformators applyinverse method to construct a GlobalValuedLocalInterpolation.
 *
 * @tparam Transformator
 * @tparam LocalInterpolation
 * @tparam Element
 */
template <class Transformator, class LocalInterpolation, class Element>
class DefaultGlobalValuedInterpolation {
    Transformator const *transformator_;
    LocalInterpolation const *localInterpolation;
    Element const *element_;

  public:
    template <class ElementInfo>
    void bind(Element const &e, ElementInfo const &eInfo, Transformator const &t,
              LocalInterpolation const &lI) {
      transformator_ = &t;
      localInterpolation = &lI;
      element_ = &e;
    }

    template <class F, class C>
    void interpolate(const F &f, std::vector<C> &out) const {
      localInterpolation->interpolate(f, out);
      transformator_->applyInverse(out);
    }
};

// The transformator class either exports a GlobalValuedInterpolation or implements an inverse
// transformation
template <class Transformator, class LocalInterpolation, class LocalBasis, class Element,
          class Enabled = bool>
struct GlobalValuedInterpolationType;

// Spezialization for Transformator class that implements a "applyInverse" method. This should only
// be used for Finite elements which are interpolation affine equivalent.
template <class Transformator, class LocalInterpolation, class LocalBasis, class Element>
struct GlobalValuedInterpolationType<
    Transformator, LocalInterpolation, LocalBasis, Element,
    std::enable_if_t<
        models<Concept::InterpolationEquivalentTransformator<Element, LocalBasis>, Transformator>(),
        bool>> {
    using type = typename Dune::Functions::Impl::DefaultGlobalValuedInterpolation<
        Transformator, LocalInterpolation, Element>;
};

// Specialization for Transformator classes which implement a GlobalValuedInterpolation class. This
// approach is suitable for non interpolation affine equivalent FEs.
template <class Transformator, class LocalInterpolation, class LocalBasis, class Element>
struct GlobalValuedInterpolationType<
    Transformator, LocalInterpolation, LocalBasis, Element,
    std::enable_if_t<
        models<Concept::InterpolationProvidingTransformator<Element, LocalBasis>,
               Transformator>(),
        bool>> {
    using type = typename Transformator::template GlobalValuedInterpolation<LocalBasis, Element>;
};

/** \brief LocalFiniteElement implementation that uses values defined wrt particular grid
 * elements
 *
 * \tparam Transformator Class implementing linear transformations
 *  \tparam LocalValuedLFE LocalFiniteElement implementation whose values are to be
 * transformed
 * \tparam Element Element onto which the FE is transformed
 */
template <class LinearTransformator, class LocalValuedLFE, class Element>
class LinearTransformedLocalFiniteElement {
    using LocalBasis =
        LinearTransformedLocalBasis<LinearTransformator,
                                    typename LocalValuedLFE::Traits::LocalBasisType, Element>;
    using LocalInterpolation = typename GlobalValuedInterpolationType<
        LinearTransformator, typename LocalValuedLFE::Traits::LocalInterpolationType,
        typename LocalValuedLFE::Traits::LocalBasisType, Element>::type;

    using ElementInformation = typename LinearTransformator::template ElementInformation<Element>;
    using LocalCoefficients =
        LinearTransformedLocalCoefficients<typename LocalValuedLFE::Traits::LocalCoefficientsType,
                                           ElementInformation>;

  public:
    /** \brief Export number types, dimensions, etc.
     */
    using Traits = LocalFiniteElementTraits<LocalBasis, LocalCoefficients, LocalInterpolation>;

    /**
     * \brief Binding routine. This binds all other objects of the LocalFiniteElement interface,
     * as well as the Transformator.
     *
     * \param localValuedLFE
     * \param element
     * \param elementInfo
     */
    void bind(const LocalValuedLFE &localValuedLFE, const Element &element,
              const ElementInformation &elementInfo) {
      element_ = element;
      transformator_.bind(element_, elementInfo);
      globalValuedLocalBasis_.bind(localValuedLFE.localBasis(), transformator_, element_,
                                   elementInfo);
      bindGlobalValuedInterpolation(localValuedLFE, element, elementInfo);
      globalValuedLocalCoefficients_.bind(elementInfo);

      localValuedLFE_ = &localValuedLFE;
    }

    void bindGlobalValuedInterpolation(const LocalValuedLFE &localValuedLFE, const Element &element,
                                       const ElementInformation &elementInfo) {
      if constexpr (models<Concept::InterpolationEquivalentTransformator<Element, LocalBasis>,
                           LinearTransformator>())
        globalValuedLocalInterpolation_.bind(element, elementInfo, transformator_,
                                             localValuedLFE.localInterpolation());
      else
        globalValuedLocalInterpolation_.bind(element_, elementInfo);
    }

    /** \brief Returns the local basis, i.e., the set of shape functions
     */
    const typename Traits::LocalBasisType &localBasis() const { return globalValuedLocalBasis_; }

    /** \brief Returns the assignment of the degrees of freedom to the element subentities
     */
    const typename Traits::LocalCoefficientsType &localCoefficients() const {
      return globalValuedLocalCoefficients_;
    }

    /** \brief Returns object that evaluates degrees of freedom
     */
    const typename Traits::LocalInterpolationType &localInterpolation() const {
      return globalValuedLocalInterpolation_;
    }

    /** \brief The number of shape functions */
    std::size_t size() const { return localValuedLFE_->size(); }

    /** \brief The reference element that the local finite element is defined on
     */
    GeometryType type() const { return localValuedLFE_->type(); }

  private:
    LinearTransformator transformator_;
    typename Traits::LocalBasisType globalValuedLocalBasis_;
    typename Traits::LocalInterpolationType globalValuedLocalInterpolation_;
    typename Traits::LocalCoefficientsType globalValuedLocalCoefficients_;
    const LocalValuedLFE *localValuedLFE_;
    Element element_;
};
} // namespace Impl
} // namespace Dune::Functions
#endif
