// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/morleybasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

using namespace Dune;
using namespace Dune::Functions;

int main(int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test("morley");

  using namespace Dune::Functions::BasisFactory;

  { // 2d
    using Grid = UGGrid<2>;
    // using Grid = YaspGrid<2>;
    auto grid = StructuredGridFactory<Grid>::createSimplexGrid({0., 0.}, {1., 1.}, {{3, 3}});

    auto gridView = grid->leafGridView();
    std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1)
              << " facettes and " << gridView.size(2) << " vertices" << std::endl;
    using GridView = decltype(gridView);
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, morley());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test.subTest(checkBasis(basis, EnableVertexContinuityCheck(),
                              EnableNormalDifferentiabilityAtMidpointsCheck(),
                              CheckLocalFiniteElementFlag()));
    }
  }

  return test.exit();
}
