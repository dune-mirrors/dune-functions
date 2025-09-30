// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>
#include <algorithm>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/functions/common/differentiablefunctionfromcallables.hh>
#include <dune/functions/functionspacebases/arnoldwintherbasis.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/albertagrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/printgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/utility/parmetisgridpartitioner.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

using namespace Dune;
using namespace Dune::Functions;


int main(int argc, char *argv[]) {
  const MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);
  auto mpiSize = mpiHelper.size();
  auto rank = mpiHelper.rank();
  Dune::TestSuite test("arnold-winther");

  using namespace Dune::Functions::BasisFactory;
  using Grid = UGGrid<2>;
  // Test with a single triangle
  {
    auto gridFactory = GridFactory<Grid>();
    gridFactory.insertVertex({0., 0.});
    gridFactory.insertVertex({1., 0.});
    gridFactory.insertVertex({0., 1.});
    // gridFactory.insertVertex({0.2, 1.});

    gridFactory.insertElement(GeometryTypes::simplex(2), {0, 1, 2});

    auto grid = gridFactory.createGrid();
    auto gridView = grid->leafGridView();
    std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1)
              << " facettes and " << gridView.size(2) << " vertices" << std::endl;
    using GridView = decltype(gridView);
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, arnoldWinther());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test.subTest(
          checkBasis(basis, EnableNormalContinuityCheck(), CheckLocalFiniteElementFlag<0>()));
    }
  }
  // Test with parallelogram
  {
    auto gridFactory = GridFactory<Grid>();
    gridFactory.insertVertex({0., 0.});
    gridFactory.insertVertex({1., 0.});
    gridFactory.insertVertex({1.2, 1.});
    gridFactory.insertVertex({0.2, 1.});

    gridFactory.insertElement(GeometryTypes::simplex(2), {0, 1, 2});
    gridFactory.insertElement(GeometryTypes::simplex(2), {0, 2, 3});

    auto grid = gridFactory.createGrid();
    auto gridView = grid->leafGridView();
    std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1)
              << " facettes and " << gridView.size(2) << " vertices" << std::endl;
    using GridView = decltype(gridView);
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, arnoldWinther());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test.subTest(
          checkBasis(basis, EnableNormalContinuityCheck(), CheckLocalFiniteElementFlag<0>()));
    }
  }
  // Test with square
  {
    auto grid = StructuredGridFactory<Grid>::createSimplexGrid({0., 0.}, {1., 1.}, {{2, 2}});
    auto gridView = grid->leafGridView();
    std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1)
              << " facettes and " << gridView.size(2) << " vertices" << std::endl;
    printGrid(*grid, mpiHelper, "testGrid");
    // using GridView = decltype(gridView);
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, arnoldWinther());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test.subTest(
          checkBasis(basis, EnableNormalContinuityCheck(), CheckLocalFiniteElementFlag<0>()));
    }
  }
  return test.exit();
}
