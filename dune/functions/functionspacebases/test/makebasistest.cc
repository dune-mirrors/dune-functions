// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <array>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/backends/istlvectorbackend.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;

  // Generate grid for testing
  const int dim = 2;
  using Grid = Dune::YaspGrid<dim>;
  Dune::FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{10, 10}};
  Grid grid(l,elements);

  using namespace Dune::Functions::BasisFactory;

  const int N = 10;
  const int M = 10;

  auto gridView = grid.leafGridView();
  auto basis = makeBasis(gridView,
      power<N>(
        power<M>(
          composite(
            lagrange<3>(),
            lagrange<1>(),
            flatLexicographic()),
          flatLexicographic()),
        blockedInterleaved())
      );

  test.subTest(checkBasis(basis, EnableContinuityCheck()));

  {
    [[maybe_unused]] auto& firstLagrangeFactor = basis.preBasis().subPreBasis().subPreBasis().subPreBasis(Dune::Indices::_0);
    [[maybe_unused]] auto& secondLagrangeFactor = basis.preBasis().subPreBasis().subPreBasis().subPreBasis<1>();
  }

  using Vector = std::vector<Dune::FieldVector<double,N>>;

  Vector x;

  auto f = [](const auto& x){
    std::array<std::array<Dune::FieldVector<double,2>, M>, N> y;
    for(auto& yi : y)
      for(auto& yij : yi)
        yij = 1.0;
    return y;
  };

  Dune::Functions::interpolate(basis, x, f);

  for(const auto& xi : x)
    for(const auto& xij : xi)
      test.require(std::abs(xij - 1.0) < 1e-10)
        << "Coefficient of interpolated 1-function does not match";

  {
    auto gridView = grid.leafGridView();
    auto basis = makeBasis(gridView,
        power<2>(
          composite(
            power<1>(power<1>(lagrange<1>())),
            power<2>(lagrange<1>()),
            power<3>(lagrange<1>())
          ),
          flatInterleaved()
        )
      );
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  {
    auto gridView = grid.leafGridView();
    auto basis = makeBasis(gridView,
        power<2>(
          composite(
            power<1>(power<1>(lagrange<1>())),
            power<2>(lagrange<1>()),
            power<3>(lagrange<1>())
          ),
          flatLexicographic()
        )
      );
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  return test.exit();
}
