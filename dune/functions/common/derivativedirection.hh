// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DERIVATIVEDIRECTION_HH
#define DUNE_FUNCTIONS_COMMON_DERIVATIVEDIRECTION_HH

#include <dune/functions/common/type_traits.hh>
#include "concept.hh"

namespace Dune {
namespace Functions {

template<int D>
struct DerivativeDirection :
    public std::integral_constant<int, D>
{};

  namespace derivativeDirection
  {
    extern DerivativeDirection<1> _d1;
    extern DerivativeDirection<2> _d2;
    extern DerivativeDirection<3> _d3;
    extern DerivativeDirection<4> _d4;
    extern DerivativeDirection<5> _d5;
    extern DerivativeDirection<6> _d6;
    extern DerivativeDirection<7> _d7;
    extern DerivativeDirection<8> _d8;
    extern DerivativeDirection<9> _d9;
    extern DerivativeDirection<10> _d10;
    extern DerivativeDirection<10> _dN;
  }

  /**
   * A concept describing types that have a derivative(f,dir) method found by ADL
   */
  template<int D>
  struct HasFreeDerivative
  {
    template<class F>
    auto require(F&& f, DerivativeDirection<D>&& d) -> decltype(
      derivative(f,d)
      );
  };

  /**
   * A concept describing types that have a derivative(f) method found by ADL
   */
  struct HasFreeDerivativeWithoutDirection
  {
    template<class F>
    auto require(F&& f) -> decltype(
      derivative(f)
      );
  };

  #warning add SFINAE to remove this from the lookup table iff Function takes only one parameter
  template<typename Function,
           typename std::enable_if<
             Concept::models< HasFreeDerivativeWithoutDirection, Function>()
             and
             not Concept::models< HasFreeDerivative<1>, Function>(), int>::type = 0>
  auto derivative(const Function & f, DerivativeDirection<1>)
    -> decltype(derivative(f))
  {
    return derivative(f);
  }

  // #warning add SFINAE to remove this from the lookup table iff Function takes multiple parameter
  // template<typename Function>
  // auto derivative(const Function & f)
  //   -> decltype(derivative(f, derivativeDirection::_d1))
  // {
  //   return derivative(f, derivativeDirection::_d1);
  // }

}
}

#endif // DUNE_FUNCTIONS_COMMON_DERIVATIVEDIRECTION_HH
