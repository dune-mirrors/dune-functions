#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMATORCONCEPTS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMATORCONCEPTS_HH

#include <dune/common/concept.hh>
#include <dune/functions/analyticfunctions/polynomial.hh>
#include <dune/functions/common/functionconcepts.hh>

#include <type_traits>
#include <vector>
/** \brief Several Concepts used throughout lineartransformedlocalfiniteelement.hh */
namespace Dune::Functions {
namespace Concept::Impl {


/** \brief A basis is an H2Basis, if it allows evaluation of the Hessians and export a corresponding Type via its Traits*/
struct H2Basis {
  template <class LB>
  auto require(LB &&lb) -> decltype(
    lb.evaluateHessian(std::declval<typename LB::Traits::DomainType>(),
      std::declval<std::vector<typename LB::Traits::Hessiantype>>())
  );
};
} // namespace Concept

// Structure to deduce the type of Hessians. Defaults to void if not exported by
// LocalBasis
namespace Impl {
template <class LVLB, class Enabled = bool>
struct HessianType {
  using type = void;
};

template <class LVLB>
struct HessianType<LVLB,
                   std::enable_if_t<models<Concept::Impl::H2Basis, LVLB>(), bool>> {
  using type = typename LVLB::Traits::HessianType;
};

  // forward declaration
template<class Transformator, class LocalValuedLFE, class Element>
class TransformedLocalFiniteElement;

template<class FE>
struct IsTransformedLocalFiniteElement
{
  static constexpr bool value = false;
};

template<class ... Args>
struct IsTransformedLocalFiniteElement<TransformedLocalFiniteElement<Args...>>
{
  static constexpr bool value = true;
};

} // namespace Impl

// Transformator concepts
namespace Concept::Impl {

/** \brief An Object is binable to an Element if it offers and bind(Element) method. */
template <class Element>
struct BindableTo {
  template <class Object>
  auto require(Object &&object)
      -> decltype(
        object.bind(std::declval<Element>())
  );
};

/** \brief A class fullfils the TransformatorBase concept for an Element and a LFE, if it is
    - Bindable to Element
    - offers transform(...) methods for the LFEs RangeType and JacobianType*/
template <class Element, class LocalFiniteElement>
struct TransformatorBase : Refines<BindableTo<Element>>{
  using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;
  using RangeType = typename LocalBasis::Traits::RangeType;
  using JacobianType = typename LocalBasis::Traits::JacobianType;
  using LocalCoordinate = typename Element::Geometry::LocalCoordinate;
  template <class T>
  auto require(T &&t)
      -> decltype(t.transform(std::declval<std::vector<RangeType> const &>(),
                              std::declval<std::vector<RangeType> &>(), std::declval<LocalCoordinate const&>()),
                  t.transform(std::declval<std::vector<JacobianType> const &>(),
                      std::declval<std::vector<JacobianType> &>(), std::declval<LocalCoordinate const&>() )
  );
};

/** A Transformator is a SizeProvidingTransformator if it offers a size() method.*/
template <class Element, class LocalFiniteElement>
struct SizeProvidingTransformator
    : Refines<TransformatorBase<Element, LocalFiniteElement>> {
  template <class T> auto require(T &&t) -> decltype(t.size());
};

/** A Transformator is a RangeSpaceTransformator if it exports a template class `LocalValuedFunction<F,LC>`.
    This class wraps a given function of type `F` and applies the inverse of the Transformation during interpolation.
    This approach is suitable for Piola-transformed finite elements.*/
template <class Element, class LocalFiniteElement>
struct RangeSpaceTransformator
    : Dune::Concept::Refines<TransformatorBase<Element, LocalFiniteElement>> {

  using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;
  using RangeType =
      typename LocalBasis::Traits::RangeType; // Range type of the underlying fe

  using LocalCoordinate = typename Element::Geometry::LocalCoordinate;
  using F = DifferentiableFunction<RangeType(LocalCoordinate)>;
  template <class T>
  auto require(T &&t)
      -> decltype(std::declval<typename T::template LocalValuedFunction<
                      F, LocalCoordinate>>()
                  //     ,
                  // // Here we assume that the Transformation does not change the
                  // // Type. This might not be true, e.g. on surfaces the Piola transform does just that.
                  // requireConcept<Function<RangeType(LocalCoordinate)>>(
                  //     std::declval<typename T::template LocalValuedFunction<
                  //         F, LocalCoordinate>>())
                          );
};

/** A Transformator is a InterpolationEquivalentTransformator if it offers an applyInverse method.
    Mathematically this corresponds to Finite Elements, where the global Dofs and the push-forwards of local Dofs
    span the same space in the infinite dimensional dual.
    In this case, we can use the inverse transformation to create a global interpolation.
    This approach is suitable for finite elements where the Dofs are only function evaluations or full gradients, but no single partial derivative a some point.
    */
template <class Element, class LocalFiniteElement>
struct InterpolationEquivalentTransformator
    : Dune::Concept::Refines<TransformatorBase<Element, LocalFiniteElement>> {
  template <class T>
  auto require(T &&t)-> decltype(
    t.applyInverse(std::declval<std::vector<double> &>(),
                  std::declval<std::vector<double> &>())
  );
};

/** A Transformator is an InterpolationProvidingTransformator if it exports a class `GlobalValuedInterpolation`, which is bindable to Element.
    This approach is suitable for all finite elements in principle.*/

template <class Element, class LocalFiniteElement>
struct InterpolationProvidingTransformator
    : Dune::Concept::Refines<TransformatorBase<Element, LocalFiniteElement>> {
  using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;
  using LocalInterpolation = typename LocalFiniteElement::Traits::LocalInterpolationType;

  using LocalCoordinate = typename Element::Geometry::LocalCoordinate;
  using RangeType = typename LocalBasis::Traits::RangeType;
  using F = DifferentiableFunction<RangeType(LocalCoordinate)>;

  template <class T>
  auto require(T &&t)-> decltype(
    requireType<typename T::GlobalValuedInterpolation>(),
    std::declval<typename T::GlobalValuedInterpolation>().bind(
          std::declval<Element>(), std::declval<LocalInterpolation>())
    // I would like to test this, but it somehow fails
    // , std::declval<typename T::template GlobalValuedInterpolation<LocalBasis, Element>>().interpolate(
    //       std::declval<F>(), std::declval<std::vector<double>&>())
  );
};

} // namespace Concept::Impl

namespace Impl{

/** \brief Implementation of a dune-localfunctions LocalInterpolation
 *    that accepts global-valued functions
 * \tparam Transformator The transformation (e.g., Piola) that transforms from local to global
 * values
 * \tparam LocalValuedLocalInterpolation The local-valued LocalInterpolation that is used for
 * the actual interpolation
 * \tparam Element The element that the global-valued FE lives on
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
 * \brief This class is a generic implementation of the GlobalValuedInterpolation
 * interface. After binding, it uses an invertible Basis Transformator and
 * the reference LocalInterpolation and to construct a GlobalValuedLocalInterpolation.
 *
 * \tparam Transformator
 * \tparam LocalInterpolation
 * \tparam Element
 */
template<class Transformator, class LocalFE, class Element>
class InterpolationEquivalentGlobalValuedLocalInterpolation
{
    friend class TransformedLocalFiniteElement<Transformator, LocalFE, Element>;
    using LocalInterpolation = typename LocalFE::Traits::LocalInterpolationType;

  public:
    InterpolationEquivalentGlobalValuedLocalInterpolation(Transformator const &t)
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

/**  Structure to deduce Traits::LocalInterpolationType of a transformed finite element
 * The transformator class either exports a LocalValuedLocalFunction, a
 * GlobalValuedInterpolation, or implements an inverse transformation. This class selects the corresponding type.
 * Primary Template - not implemented
 */
template<class Transformator, class LocalFE, class Element, class Enabled = bool>
struct GlobalValuedInterpolationType;

// Spezialization for Transformator class that implements a "applyInverse"
// method. This should only be used for Finite elements which are interpolation
// affine equivalent.
template<class Transformator, class LocalFE, class Element>
struct GlobalValuedInterpolationType<
    Transformator, LocalFE, Element,
    std::enable_if_t<
        models<Concept::Impl::InterpolationEquivalentTransformator<Element, LocalFE>, Transformator>(),
        bool>> {
    using type =
        InterpolationEquivalentGlobalValuedLocalInterpolation<Transformator, LocalFE, Element>;
};

// Specialization for Transformator classes which implement a
// GlobalValuedInterpolation class. This approach is suitable for non (!)
// interpolation affine equivalent FEs.
template<class Transformator, class LocalFE, class Element>
struct GlobalValuedInterpolationType<
    Transformator, LocalFE, Element,
    std::enable_if_t<
        models<Concept::Impl::InterpolationProvidingTransformator<Element, LocalFE>, Transformator>(),
        bool>> {
    using type = typename Transformator::GlobalValuedInterpolation;
};

// Specialization for Transformator classes which implement a
// LocalValuedFunction class. This approach is suitable for FEs with a
// Rangespace transforming pullback like Piola transformations.
template<class Transformator, class LocalFE, class Element>
struct GlobalValuedInterpolationType<
    Transformator, LocalFE, Element,
    std::enable_if_t<models<Concept::Impl::RangeSpaceTransformator<Element, LocalFE>, Transformator>(),
                     bool>> {
    using type = RangeSpaceTransformingGlobalValuedLocalInterpolation<
        Transformator, typename LocalFE::Traits::LocalInterpolationType, Element>;
};
}// namespace Impl
} // namespace Dune::Functions
#endif
