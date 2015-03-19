// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DERIVATIVEDIRECTION_HH
#define DUNE_FUNCTIONS_COMMON_DERIVATIVEDIRECTION_HH

#include <dune/functions/common/type_traits.hh>

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

#ifdef DOXYGEN
  template<int D>
  DerivativeType derivative(const Function & f, DerivativeDirection<D>);
#endif

}
}

#endif // DUNE_FUNCTIONS_COMMON_DERIVATIVEDIRECTION_HH
