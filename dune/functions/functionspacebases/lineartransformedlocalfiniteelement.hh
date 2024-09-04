// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LINEARTRANSFORMEDLOCALFINITEELEMENT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LINEARTRANSFORMEDLOCALFINITEELEMENT_HH

#include <array>
#include <dune/common/std/type_traits.hh>
#include <dune/functions/functionspacebases/globalvaluedlocalfiniteelement.hh>
#include <dune/functions/functionspacebases/transformatorconcepts.hh>
#include <dune/functions/functionspacebases/nodes.hh>
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
 * element, unlike the Piola-transformations, which can be applied to a set of
 * admissible finite elements.
 *
 * Some details:
 * A transformation can reduce the size of a finite element (like in the case of the Bell element)
 * In such case, the transformator has to implement a size() method, which will be
 * used accordingly by all objects. Currently, this implies that the size() methods of the
 * TransformedLocal... classes are runtime members.
 * Generally, we cannot use the reference dofs for interpolation, so three different
 * approaches are offered:
 * - Wrapping the function with the inverse of a rangespace transformation and
 *    using the reference dofs (e.g. for Piola)
 * - Using the reference dofs and the inverse of the basis transformation (For
 *    interpolation affine equivalent FEs)
 * - Using a independent global interpolation provided by the Transformator.
 *
 * TransformedLocalBasis:
 *  - Maintains the level of derivatives implemented by the reference basis.
 *  - Offers access the the wrapped basis and the transformator object.
 *
 * TransformedLocalCoefficients:
 *  - wraps the LocalCoefficients of the localvalued finite element
 *  - can be bound to ElementInformation
 *
 * TransformedLocalFiniteElement:
 *  - owns transformator, generic reference FE, and transformed LocalBasis/Interpolation/Coefficients
 *  - Binds all its bindable subobjects whenever bound to an element
 * TransformedNode:
 *  - Generic Implementation.
 *
 * Limitations:
 * You cannot chain TransformedFiniteElements, but you can chain the transformator
 * Currently the concept check only allows for Rangespace transformations that maintain the type
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
 * transformed
 * \tparam Element The element that the global-valued FE lives on
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
      if constexpr (models<Concept::Impl::SizeProvidingTransformator<Element, LocalValuedLFE>,
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
      transformator_->transform(rangeBuffer_, out, x);
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
      transformator_->transform(jacobianBuffer_, out, x);
    }

    /** \brief Evaluate Hessian of all shape functions
     *   \note Sfinae protected: this method is only available, if the wrapped
     * LocalBasis 1. exports a <HessianType> and 2. provides a evaluateHessian
     * method with a corresponding signature \param x Point in the reference
     * element where to evaluation the Hessians \param[out] out The Hessians of
     * all shape functions at the point x
     */
    template<class TT = LocalValuedLocalBasis,
             std::enable_if_t<models<Concept::Impl::H2Basis, TT>(), int> = 0>
    void evaluateHessian(const typename Traits::DomainType &x,
                         std::vector<typename TT::HessianType> &out) const
    {

      thread_local std::vector<typename TT::HessianType> hessianBuffer_;
      static_assert(std::is_same_v<TT, LocalValuedLocalBasis>);
      hessianBuffer_.resize(localValuedLocalBasis_->size());
      localValuedLocalBasis_->evaluateHessian(x, hessianBuffer_);
      out.resize(size());
      transformator_->transform(hessianBuffer_, out, x);
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     *
     * \param order Order of the partial derivatives, in the classic multi-index
     * notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out The desired partial derivatives
     */
    void partial(std::array<unsigned int, Traits::dimDomain> const &order,
                 const typename Traits::DomainType &x,
                 std::vector<typename Traits::RangeType> &out) const
    {

      thread_local std::vector<LocalValuedRangeType> rangeBuffer_;
      rangeBuffer_.resize(localValuedLocalBasis_->size());
      localValuedLocalBasis_->partial(order, x, rangeBuffer_);
      out.resize(size());
      transformator_->transform(rangeBuffer_, out, x);
    }

    //! \brief Polynomial order of the shape functions
    auto order() const { return localValuedLocalBasis_->order(); }

    Transformator const &transformator() const { return *transformator_; }

    LocalValuedLocalBasis const &cacheable() const { return *localValuedLocalBasis_; }

  private:
    Transformator const *transformator_;
    LocalValuedLocalBasis const *localValuedLocalBasis_;
};

/** \brief Associations of the transformed degrees of freedom to subentities of
 * the reference simplex. \note We assume here, that the transformation does not
 * change the LocalCoefficients, i.e. it does not change to which subentity a
 * DoF is associated with. However, it may reduce the size (e.g. for the reduced
 * Hermite Element / DKT).
 *
 * \tparam Transformator. The Transformation to apply to the finite element.
 * \tparam LocalValuedLFE. The local finite element to be transformed
 * \tparam Element. The Gridelement type
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
      if constexpr (models<Concept::Impl::SizeProvidingTransformator<Element, LocalValuedLFE>,
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
    static_assert(not Impl::IsTransformedLocalFiniteElement<LocalValuedLFE>::value);

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
          localValuedLFE_(std::forward<Args>(args)...),
          globalValuedLocalBasis_(transformator_),
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
      // In particular, this allows to wrap a FEVariant & FEMap combiniation
      if constexpr (models<Concept::Impl::BindableTo<Element>, LocalValuedLFE>())
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
      if constexpr (models<Concept::Impl::SizeProvidingTransformator<Element, LocalValuedLFE>,
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
    Element const *element_ = nullptr;
};

template<typename GV, class Transformator, typename LocalValuedLocalFiniteElement>
class TransformedNode : public LeafBasisNode
{
    static_assert(not IsTransformedLocalFiniteElement<LocalValuedLocalFiniteElement>::value);
    static constexpr unsigned int dim = GV::dimension;

  public:
    using size_type = std::size_t;
    using Element = typename GV::template Codim<0>::Entity;

    using FiniteElement =
        Impl::TransformedLocalFiniteElement<Transformator, LocalValuedLocalFiniteElement, Element>;

    template<class... GlobalState>
    TransformedNode(GlobalState&&... globalState)
    : finiteElement_(std::forward<GlobalState>(globalState)...)
    {
      // finiteElement_ is not bound yet, i.e. it might not have a size
    }

    //TODO Copy Constr, Move, assignment, etc

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

    //! The order of the local basis.
    unsigned int order() const { return finiteElement_.localBasis().order(); }

  protected:
    FiniteElement finiteElement_;
    Element element_;
};
} // namespace Impl

} // namespace Dune::Functions
#endif
