// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <dune/common/test/testsuite.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/uggrid.hh>

#include <dune/functions/functionspacebases/arnoldwintherbasis.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>

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

template<class Basis>
Dune::TestSuite testEdgeDOFNumbering(Basis const& basis)
{
  Dune::TestSuite test("Arnold-Winther edge DOF numbering");
  auto insideLocalView = basis.localView();
  auto outsideLocalView = basis.localView();
  const auto& indexSet = basis.gridView().indexSet();
  const auto& idSet = basis.gridView().grid().globalIdSet();

  for (const auto& inside : elements(basis.gridView())) {
    insideLocalView.bind(inside);
    for (const auto& intersection : intersections(basis.gridView(), inside)) {
      if (not intersection.neighbor())
        continue;

      const auto outside = intersection.outside();
      if (indexSet.index(inside) >= indexSet.index(outside))
        continue;

      outsideLocalView.bind(outside);
      const auto insideEdge = intersection.indexInInside();
      const auto outsideEdge = intersection.indexInOutside();
      const auto insideOrientations = basis.preBasis().faceOrientations(inside);
      const auto outsideOrientations = basis.preBasis().faceOrientations(outside);
      const bool insideFlipped = insideOrientations.faceOrientationIndex(insideEdge, 1);
      const bool outsideFlipped = outsideOrientations.faceOrientationIndex(outsideEdge, 1);

      auto referenceMomentVertexId = [&](const auto& element, auto edge,
                                         bool flipped, std::size_t moment) {
        const auto& referenceElement = Dune::referenceElement(element);
        // A flipped edge swaps moment positions, but not the two tensor
        // components stored at each position.
        const auto referenceMoment = flipped ? 1 - moment : moment;
        const auto vertex = referenceElement.subEntity(edge, 1, referenceMoment, 2);
        return idSet.subId(element, vertex, 2);
      };

      for (std::size_t moment = 0; moment < 2; ++moment) {
        test.check(referenceMomentVertexId(inside, insideEdge, insideFlipped, moment)
                   == referenceMomentVertexId(outside, outsideEdge, outsideFlipped, moment))
            << "FaceOrientations assigns different vertices to edge moment " << moment
            << " (inside edge " << insideEdge << ", flip " << insideFlipped
            << "; outside edge " << outsideEdge << ", flip " << outsideFlipped << ")";

        for (std::size_t component = 0; component < 2; ++component) {
          const auto edgeDOF = 2 * moment + component;
          const auto insideLocalDOF = 9 + 4 * insideEdge + edgeDOF;
          const auto outsideLocalDOF = 9 + 4 * outsideEdge + edgeDOF;
          test.check(insideLocalView.index(insideLocalDOF)
                     == outsideLocalView.index(outsideLocalDOF))
              << "Edge moment " << moment << ", tensor component " << component
              << " has inconsistent global basis indices (inside edge " << insideEdge
              << ", flip " << insideFlipped << "; outside edge " << outsideEdge
              << ", flip " << outsideFlipped << ")";
        }
      }
    }
  }
  return test;
}

int main(int argc, char *argv[]) {
  Dune::MPIHelper::instance(argc, argv);
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

      test.subTest(testEdgeDOFNumbering(basis));
      test.subTest(
          checkBasis(basis, EnableNormal_VectorContinuityCheck()));

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
      test.subTest(testEdgeDOFNumbering(basis));
      test.subTest(checkBasis(basis, EnableNormal_VectorContinuityCheck()));
    }

    grid->globalRefine(1);

    {
      auto gridView = grid->leafGridView();
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, arnoldWinther());
      test.subTest(testEdgeDOFNumbering(basis));
      test.subTest(checkBasis(basis, EnableNormal_VectorContinuityCheck()));
    }

  }
  return test.exit();
}
