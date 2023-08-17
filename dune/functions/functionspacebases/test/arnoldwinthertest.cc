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

std::vector<unsigned int> &permute(std::vector<unsigned int> &indices) {
  static unsigned int counter = 0;
  for (auto i = 0u; i < counter % indices.size(); ++i)
    std::next_permutation(indices.begin(), indices.end());
  return indices;
}

template <int dim>
auto createPerturbedGrid(unsigned int elementsPerDim, int rank, double relativeShift = 0.01,
                         bool perturbBoundary = false) {
  // std::srand(std::time(nullptr)); // use current time as seed for random generator
  using namespace Dune::Functions::BasisFactory;
  if constexpr (dim == 1) {
    using Grid = OneDGrid;

    auto gridFactory = GridFactory<Grid>();
    if (rank == 0) {
      Dune::FieldVector<double, dim> start = {0.};
      Dune::FieldVector<double, dim> end = {1.};

      gridFactory.insertVertex(start);
      Dune::FieldVector<double, dim> shiftDirection = 1.;
      double shiftFactor = relativeShift / elementsPerDim;

      auto shiftPerColumn = 1. / elementsPerDim * (end - start);
      for (unsigned int i = 1; i < elementsPerDim; ++i, shiftDirection *= -1) {
        auto newVertex = start + i * shiftPerColumn + shiftDirection * shiftFactor;
        gridFactory.insertVertex(newVertex);
        std::vector<unsigned int> indices = {i - 1, i};

        gridFactory.insertElement(GeometryTypes::simplex(dim), permute(indices));
      }
      gridFactory.insertVertex(end);
      gridFactory.insertElement(GeometryTypes::simplex(dim), {elementsPerDim, elementsPerDim - 1});
    }
    return gridFactory.createGrid();
  } else if constexpr (dim == 2) {
    using Grid = UGGrid<dim>;

    auto gridFactory = GridFactory<Grid>();
    if (rank == 0) {
      Dune::FieldVector<double, dim> lowerLeft = {0., 0.};
      Dune::FieldVector<double, dim> lowerRight = {1., 0.};
      Dune::FieldVector<double, dim> upperLeft = {0., 1.};
      Dune::FieldVector<double, dim> upperRight = {1., 1.};
      auto shiftPerRow = 1. / elementsPerDim * (upperLeft - lowerLeft);
      Dune::FieldVector<double, dim> shiftDirection = {1., 1.};
      for (unsigned int row = 0; row < elementsPerDim + 1; ++row) {
        double shiftFactor = relativeShift / elementsPerDim;
        double rotation = 10.; // so long as it is not related to pi

        Dune::FieldMatrix<double, dim, dim> rotationMatrix = {
            {std::cos(rotation), -std::sin(rotation)}, {std::sin(rotation), std::cos(rotation)}};
        Dune::FieldVector<double, dim> tmp;

        rotationMatrix.mv(shiftDirection, tmp);
        shiftDirection = tmp;
        Dune::FieldVector<double, dim> start = lowerLeft + row * shiftPerRow;
        Dune::FieldVector<double, dim> end = lowerRight + row * shiftPerRow;
        auto shiftPerColumn = 1. / elementsPerDim * (end - start);
        if (perturbBoundary)
          gridFactory.insertVertex(start + shiftDirection * shiftFactor);
        else
          gridFactory.insertVertex(start);

        for (std::size_t i = 1; i < elementsPerDim;
             ++i, rotationMatrix.mv(shiftDirection, tmp), shiftDirection = tmp) {
          auto newVertex = start + i * shiftPerColumn;
          if (perturbBoundary || (row > 0 && row < elementsPerDim))
            newVertex += shiftDirection * shiftFactor;
          gridFactory.insertVertex(newVertex);
        }
        if (perturbBoundary)
          end += shiftDirection * shiftFactor;
        gridFactory.insertVertex(end);
      }
      // insert elements
      // This for loop iterates over element rows
      for (unsigned int row = 0; row < elementsPerDim; ++row) {

        for (unsigned int col = 0; col < elementsPerDim; ++col) {
          // lower triangle of rect
          std::vector<unsigned int> indices = {row * (elementsPerDim + 1) + col,
                                               (row + 1) * (elementsPerDim + 1) + col,
                                               (row + 1) * (elementsPerDim + 1) + (col + 1)};
          gridFactory.insertElement(GeometryTypes::simplex(dim), permute(indices));
          indices = {row * (elementsPerDim + 1) + col, (row) * (elementsPerDim + 1) + (col + 1),
                     (row + 1) * (elementsPerDim + 1) + (col + 1)};
          gridFactory.insertElement(GeometryTypes::simplex(dim), permute(indices));
        }
      }
    }
    return gridFactory.createGrid();
  }
  // else if constexpr (dim == 3)
  // {
  // }
  DUNE_THROW(Dune::NotImplemented, "perturbed grid for desired dimension " + std::to_string(dim));
}

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
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, arnoldWinther());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test.subTest(
          checkBasis(basis, EnableNormalContinuityCheck(), CheckLocalFiniteElementFlag<0>()));
    }
  }
  // Test with perturbed grid
  {
    auto grid = createPerturbedGrid<2>(4, 0, 0.2, true);
    grid->globalRefine(1);
    auto gridView = grid->leafGridView();
    std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1)
              << " facettes and " << gridView.size(2) << " vertices" << std::endl;
    using GridView = decltype(gridView);
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, arnoldWinther());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test.subTest(checkBasis(
          basis, EnableNormalContinuityCheck(),
          CheckLocalFiniteElementFlag<0>())); // TODO once GVLFE is fixed to treat partials of
                                              // tensors correctly, the template param here should
                                              // be set to <1>
    }
  }

  return test.exit();
}
