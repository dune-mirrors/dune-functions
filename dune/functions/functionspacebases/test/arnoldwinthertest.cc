// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>
#include <algorithm>

#include <dune/common/test/testsuite.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/rangegenerators.hh>
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
#include <stdexcept>
#include <string>

using namespace Dune;
using namespace Dune::Functions;

template<class Basis, class Interpolation>
Dune::TestSuite testDeltaProperty(Basis const& basis, Interpolation const& interpolation){
  Dune::TestSuite test("Local Test");

  using Traits = typename Basis::Traits;
  double eps = 1e-12;
  for (auto i : Dune::range(basis.size()))
  {
    auto f  =[i, &b = basis](auto const& x){
      std::vector<typename Traits::RangeType> values(b.size());
      b.evaluateFunction(x, values);
      return values[i];
    };
    std::vector<double> coeffs(basis.size());
    interpolation.interpolate(f, coeffs);

    for (auto j : Dune::range(basis.size())){
      test.check(std::abs(int(i==j) - coeffs[j])< eps,"Delta check for functional " + std::to_string(j)+ " on shapefunction " + std::to_string(i) + " returned " + std::to_string(coeffs[j])+ ", but "+std::to_string(int(i==j)) + " was expected");
    }
  }
  return test;
}

int main(int argc, char *argv[]) {
  const MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);
  auto mpiSize = mpiHelper.size();
  auto rank = mpiHelper.rank();
  Dune::TestSuite test("arnold-winther");
  std::cout<<"Testing AW reference finite element"<<std::endl;
  // first test the plain reference basis and interpolation
  using Basis = Dune::Functions::Impl::ArnoldWintherReferenceLocalBasis<double, double>;
  Basis basis;
  Dune::Functions::Impl::ArnoldWintherReferenceLocalInterpolation<double, double> interpolation;
  test.subTest(testDeltaProperty(basis, interpolation));

  std::cout<<"Testing AW finite element on grid with one element"<<std::endl;
  // Second test with transfromed basis and global interpolation
  using namespace Dune::Functions::BasisFactory;
  using Grid = UGGrid<2>;
  // test on a Grid with one triangle
  {
    auto gridFactory = GridFactory<Grid>();
    gridFactory.insertVertex({0., 0.});
    gridFactory.insertVertex({1., 0.});
    gridFactory.insertVertex({0., 1.});

    gridFactory.insertElement(GeometryTypes::simplex(2), {0,2, 1});

    auto grid = gridFactory.createGrid();
    auto gridView = grid->leafGridView();

    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, arnoldWinther());

      test.subTest(
          checkBasis(basis, EnableNormal_VectorContinuityCheck(), CheckLocalFiniteElementFlag<0>()));

      std::cout<<"Edge orientations: \n";
      for (auto && bitset :  basis.preBasis().data_)
          std::cout<<bitset<<std::endl;
    }
  }

  std::cout<<"Testing AW finite element on grid with two elements"<<std::endl;
   // Test with parallelogram
  {
    auto gridFactory = GridFactory<Grid>();
    gridFactory.insertVertex({0., 0.});
    gridFactory.insertVertex({1.1, 0.});
    gridFactory.insertVertex({1.1, 1.});
    gridFactory.insertVertex({0., 1.});

    gridFactory.insertElement(GeometryTypes::simplex(2), {0, 1, 2});
    gridFactory.insertElement(GeometryTypes::simplex(2), {2, 0, 3});

    auto grid = gridFactory.createGrid();

    {
      auto gridView = grid->leafGridView();
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, arnoldWinther());
      test.subTest(checkBasis(basis, EnableNormal_VectorContinuityCheck(), CheckLocalFiniteElementFlag<0>()));
    }

    grid->globalRefine(1);

    {
      auto gridView = grid->leafGridView();
      Dune::printGrid(gridView.grid(), Dune::MPIHelper::instance(), "grid");
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, arnoldWinther());
      test.subTest(checkBasis(basis, EnableNormal_VectorContinuityCheck(), CheckLocalFiniteElementFlag<0>()));
      std::cout<<"Edge orientations: \n";
      for (auto && bitset :  basis.preBasis().data_)
          std::cout<<bitset<<std::endl;

    }

  }
  return test.exit();
}
