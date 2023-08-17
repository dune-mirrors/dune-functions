#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMATORCONCEPTS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMATORCONCEPTS_HH

#include "dune/functions/common/functionconcepts.hh"
#include <dune/common/concept.hh>
#include <dune/functions/analyticfunctions/polynomial.hh>
#include <type_traits>
#include <vector>

namespace Dune::Functions::Concept {
struct H2Basis {
    template <class LB>
    auto require(LB &&lb) -> decltype(lb(std::declval<typename LB::Traits::DomainType>(),
                                         std::declval<std::vector<typename LB::Hessiantype>>()));
};

template <class Element, class LocalBasis>
struct LinearTransformatorBase {
    template <class T>
    auto require(T &&t)
        -> decltype(
          t.bind(std::declval<Element>(), std::declval<typename T::template ElementInformation<Element>>()),
          t.apply(std::declval<std::vector<double> &>()),
          t.apply(std::declval<std::vector<typename LocalBasis::Traits::JacobianType>&>())
          );
};

template <class Element, class LocalBasis>
struct InterpolationEquivalentTransformator
    : Dune::Concept::Refines<LinearTransformatorBase<Element, LocalBasis>> {
    template <class T>
    auto require(T &&t) -> decltype(
      t.applyInverse(std::declval<std::vector<double>&>())); // TODO what is the correct type here?
};

template <class Element, class LocalBasis>
struct InterpolationProvidingTransformator
    : Dune::Concept::Refines<LinearTransformatorBase<Element, LocalBasis>> {
    template <class T>
    auto require(T &&t) -> decltype(
          std::declval<typename T::template GlobalValuedInterpolation<LocalBasis, Element>>(),
          std::declval<typename T::template GlobalValuedInterpolation<LocalBasis, Element>>().bind(
            std::declval<Element>(), std::declval<typename T::template ElementInformation<Element>>()),
          std::declval<typename T::template GlobalValuedInterpolation<LocalBasis, Element>>().interpolate(std::declval<DifferentiableFunction<typename LocalBasis::Traits::RangeType(typename Element::Geometry::LocalCoordinate)>>(),
               std::declval<std::vector<double>&>())
          );
};
} // namespace Dune::Functions::Concept

#endif