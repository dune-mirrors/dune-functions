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

template<class T>
Dune::TestSuite testBasis(T arg, const MPIHelper &mpiHelper){
  Dune::TestSuite test("Global test");

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

    using GridView = decltype(gridView);
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, arnoldWinther());

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

    using GridView = decltype(gridView);
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, arnoldWinther());

      test.subTest(
          checkBasis(basis, EnableNormalContinuityCheck(), CheckLocalFiniteElementFlag<0>()));
    }
  }
  // Test with square
  {
    auto grid = StructuredGridFactory<Grid>::createSimplexGrid({0., 0.}, {1., 1.}, {{2, 2}});
    auto gridView = grid->leafGridView();

    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, arnoldWinther());

      test.subTest(
          checkBasis(basis, EnableNormalContinuityCheck(), CheckLocalFiniteElementFlag<0>()));
    }
  }
  return test;
}

template<class Basis, class Interpolation>
Dune::TestSuite testFiniteElement(Basis const& basis, Interpolation const& interpolation){
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
  test.subTest(testFiniteElement(basis, interpolation));

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

    gridFactory.insertElement(GeometryTypes::simplex(2), {1, 2, 0});

    auto grid = gridFactory.createGrid();
    auto gridView = grid->leafGridView();

    using GridView = decltype(gridView);
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, arnoldWinther());
      auto localView = basis.localView();
      for (auto const& e : elements(gridView))
      {
        localView.bind(e);
        test.subTest(testFiniteElement(localView.tree().finiteElement().localBasis(), localView.tree().finiteElement().localInterpolation()));
      }

      // test.subTest(
      //     checkBasis(basis, EnableNormalContinuityCheck(), CheckLocalFiniteElementFlag<0>()));
      std::cout<<"Edge orientations: \n";
      for (auto && bitset :  basis.preBasis().data_)
          std::cout<<bitset<<std::endl;
      auto& indexSet = gridView.indexSet();
      for (auto&& e : elements(gridView))
      {
        auto ref = referenceElement(e);
        for (auto i : Dune::range(3))
        {
          auto subgeo = e.subEntity<1>(i).geometry();
          std::cout<<"Tangent "<< i<<" : "<<subgeo.corner(1) - subgeo.corner(0)<<std::endl;
          // std::cout<<"dx: "<<subgeo.integrationElement(0)<<" detJ: "<<subgeo.jacobianTransposed(0)*subgeo.jacobian(0)<<std::endl;
          std::cout<<"Corners: "<<indexSet.subIndex(e.subEntity<1>(i),0,2)<<indexSet.subIndex(e.subEntity<1>(i),1,2)<<std::endl;
        }
      }
    }
    // printGrid(*grid, "testGrid");
  }
  return test.exit();
  std::cout<<"Testing AW finite element on grid with two elements"<<std::endl;
   // Test with parallelogram
  {
    auto gridFactory = GridFactory<Grid>();
    gridFactory.insertVertex({0., 0.});
    gridFactory.insertVertex({1., 0.});
    gridFactory.insertVertex({1., 1.});
    gridFactory.insertVertex({0., 1.});

    gridFactory.insertElement(GeometryTypes::simplex(2), {0, 1, 2});
    gridFactory.insertElement(GeometryTypes::simplex(2), {2, 0, 3});

    auto grid = gridFactory.createGrid();
    auto gridView = grid->leafGridView();

    using GridView = decltype(gridView);
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, arnoldWinther());
      auto localView = basis.localView();
      for (auto const& e : elements(gridView))
      {
        localView.bind(e);
        test.subTest(testFiniteElement(localView.tree().finiteElement().localBasis(), localView.tree().finiteElement().localInterpolation()));
      }
      // test.subTest(checkBasis(basis, EnableNormalContinuityCheck(), CheckLocalFiniteElementFlag<0>()));

    int twist = 0;
    auto&& indexSet = gridView.grid().globalIdSet();
    auto edgeIndexToVertices = [](int edge)->std::array<int,2>
    {
      if (edge == 0) return {0,1};
      else if (edge == 1) return {0,2};
      else if (edge == 2) return {1,2};
      else
       throw std::runtime_error("Invalid edge index");
    };

    for (auto&& e : elements(gridView))
      for (auto&& i : intersections(gridView, e))
      {
        auto refElement = referenceElement(e);
        if (i.neighbor())
        {
        // Local vertex indices within the element
        auto localV0 = refElement.subEntity(i.indexInInside(), 1, 0, 2);
        auto localV1 = refElement.subEntity(i.indexInInside(), 1, 1, 2);
        auto globalV0 = indexSet.subId(i.inside(), localV0, 2);
        auto globalV1 = indexSet.subId(i.inside(), localV1, 2);
        bool insideFlip = globalV0 > globalV1;
        localV0 = refElement.subEntity(i.indexInOutside(), 1, 0, 2);
        localV1 = refElement.subEntity(i.indexInOutside(), 1, 1, 2);
        globalV0 = indexSet.subId(i.outside(), localV0, 2);
        globalV1 = indexSet.subId(i.outside(), localV1, 2);
        bool outsideFlip = globalV0 > globalV1;
        // The grid is twisted if inside and outside are flipped differently
        twist += insideFlip != outsideFlip;

        }
      }
    std::cout<<"Grid is twisted "<<twist/2<<" times\n";
    std::cout<<"Edge orientations: \n";
    for (auto && bitset :  basis.preBasis().data_)
      std::cout<<bitset<<std::endl;
    }
  }
  return test.exit();
}
