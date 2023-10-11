#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMATORCONCEPTS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMATORCONCEPTS_HH

#include <dune/common/concept.hh>
#include <dune/functions/analyticfunctions/polynomial.hh>
#include <dune/functions/common/functionconcepts.hh>

#include <type_traits>
#include <vector>

namespace Dune::Functions {
namespace Concept {

template <class Element> struct BindableTo {
  template <class Object>
  auto require(Object &&object)
      -> decltype(
        object.bind(std::declval<Element>())
  );
};

struct H2Basis {
  template <class LB>
  auto require(LB &&lb) -> decltype(
    lb.evaluateHessian(std::declval<typename LB::Traits::DomainType>(),
      std::declval<std::vector<typename LB::Traits::Hessiantype>>())
  );
};
} // namespace Concept

// Structure to decude the type of Hessians. Defaults to void if not exported by
// LocalBasis
namespace Impl {
template <class LVLB, class Enabled = bool> struct HessianType {
  using type = void;
};

template <class LVLB>
struct HessianType<LVLB,
                   std::enable_if_t<models<Concept::H2Basis, LVLB>(), bool>> {
  using type = typename LVLB::Traits::HessianType;
};
} // namespace Impl

// Transformator concepts
namespace Concept {

template <class Element, class LocalFiniteElement>
struct TransformatorBase {
  using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;
  using JacobianType = typename LocalBasis::Traits::JacobianType;
  template <class T>
  auto require(T &&t)
      -> decltype(t.bind(std::declval<Element>()),
                  t.transform(std::declval<std::vector<double> const &>(),
                              std::declval<std::vector<double> &>()),
                  t.transform(std::declval<std::vector<JacobianType> const &>(),
                      std::declval<std::vector<JacobianType> &>()),
                  t.makeGlobalValuedInterpolation()
  );
};

template <class Element, class LocalFiniteElement>
struct SizeProvidingTransformator
    : Refines<TransformatorBase<Element, LocalFiniteElement>> {
  template <class T> auto require(T &&t) -> decltype(t.size());
};


// The following concepts are not enforced
template <class Element, class LocalFiniteElement>
struct InterpolationEquivalentTransformator
    : Dune::Concept::Refines<TransformatorBase<Element, LocalFiniteElement>> {
  template <class T>
  auto require(T &&t)-> decltype(
    t.applyInverse(std::declval<std::vector<double> &>(),
                  std::declval<std::vector<double> &>())
  );
};

template <class Element, class LocalFiniteElement>
struct InterpolationProvidingTransformator
    : Dune::Concept::Refines<TransformatorBase<Element, LocalFiniteElement>> {
  using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;
  using LocalInterpolation = typename LocalFiniteElement::Traits::LocalInterpolationType;

  using LocalCoordinate = typename Element::Geometry::LocalCoordinate;
  using RangeType = typename LocalBasis::Traits::RangeType;

  template <class T>
  auto require(T &&t)-> decltype(
    requireType<typename T::GlobalValuedInterpolation>(),
    std::declval<typename T::GlobalValuedInterpolation>().bind(
          std::declval<Element>(), std::declval<LocalInterpolation>())
                  /*,
                  std::declval<typename T::template
                  GlobalValuedInterpolation<LocalBasis, Element>>().interpolate(
                    std::declval<DifferentiableFunction<RangeType(LocalCoordinate)>>(),
                       std::declval<std::vector<double>&>())*/
  );
};

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
      -> decltype(std::declval<typename T::template LocalFunction<
                      F, LocalCoordinate, Element>>(),
                  // Here we assume that the Transformation does not change the
                  // Type. This might not be true, e.g. on surfaces
                  requireConcept<Function<RangeType(LocalCoordinate)>>(
                      std::declval<typename T::template LocalFunction<
                          F, LocalCoordinate, Element>>()));
};
} // namespace Concept
} // namespace Dune::Functions
#endif
