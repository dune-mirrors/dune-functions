// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <dune/common/test/testsuite.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/geometrygrid.hh>
#include <dune/grid/uggrid.hh>

#include <dune/functions/functionspacebases/arnoldwintherbasis.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;
using namespace Dune::Functions;

class AffineSurfaceCoordinates
    : public Dune::AnalyticalCoordFunction<double, 2, 3, AffineSurfaceCoordinates>
{
public:
  void evaluate(const DomainVector& x, RangeVector& y) const
  {
    y = {x[0] + 0.2 * x[1], 0.3 * x[0] + x[1], 0.4 * x[0] - 0.2 * x[1]};
  }
};

struct NonAffineGeometryStub
{
  using ctype = double;
  using LocalCoordinate = FieldVector<double, 2>;
  using GlobalCoordinate = FieldVector<double, 2>;
  static constexpr int mydimension = 2;
  static constexpr int coorddimension = 2;
  bool affine() const
  {
    return false;
  }
  GeometryType type() const
  {
    return GeometryTypes::simplex(2);
  }
  FieldMatrix<double, 2, 2> jacobian(const LocalCoordinate&) const
  {
    return {{1, 0}, {0, 1}};
  }
  double integrationElement(const LocalCoordinate&) const
  {
    return 1;
  }
};

struct NonAffineElementStub
{
  using Geometry = NonAffineGeometryStub;
  Geometry geometry() const
  {
    return {};
  }
};

template <class Basis, class Interpolation>
Dune::TestSuite testDeltaProperty(Basis const& basis, Interpolation const& interpolation)
{
  Dune::TestSuite test("Local Test");

  using Traits = typename Basis::Traits;
  double eps = 1e-12;
  for (auto i : Dune::range(basis.size())) {
    auto f = [i, &b = basis](auto const& x) {
      std::vector<typename Traits::RangeType> values(b.size());
      b.evaluateFunction(x, values);
      return values[i];
    };
    std::vector<double> coeffs(basis.size());
    interpolation.interpolate(f, coeffs);

    for (auto j : Dune::range(basis.size())) {
      test.check(std::abs(int(i == j) - coeffs[j]) < eps,
                 "Delta check for functional " + std::to_string(j) + " on shape function " +
                     std::to_string(i) + " returned " + std::to_string(coeffs[j]) + ", but " +
                     std::to_string(int(i == j)) + " was expected");
    }
  }
  return test;
}

template <class Basis>
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

      auto referenceMomentVertexId = [&](const auto& element, auto edge, bool flipped,
                                         std::size_t moment) {
        const auto& referenceElement = Dune::referenceElement(element);
        // A flipped edge swaps moment positions, but not the two tensor
        // components stored at each position.
        const auto referenceMoment = flipped ? 1 - moment : moment;
        const auto vertex = referenceElement.subEntity(edge, 1, referenceMoment, 2);
        return idSet.subId(element, vertex, 2);
      };

      for (std::size_t moment = 0; moment < 2; ++moment) {
        test.check(referenceMomentVertexId(inside, insideEdge, insideFlipped, moment) ==
                   referenceMomentVertexId(outside, outsideEdge, outsideFlipped, moment))
            << "FaceOrientations assigns different vertices to edge moment " << moment
            << " (inside edge " << insideEdge << ", flip " << insideFlipped << "; outside edge "
            << outsideEdge << ", flip " << outsideFlipped << ")";

        for (std::size_t component = 0; component < 2; ++component) {
          const auto edgeDOF = 2 * moment + component;
          const auto insideLocalDOF = 9 + 4 * insideEdge + edgeDOF;
          const auto outsideLocalDOF = 9 + 4 * outsideEdge + edgeDOF;
          test.check(insideLocalView.index(insideLocalDOF) ==
                     outsideLocalView.index(outsideLocalDOF))
              << "Edge moment " << moment << ", tensor component " << component
              << " has inconsistent global basis indices (inside edge " << insideEdge << ", flip "
              << insideFlipped << "; outside edge " << outsideEdge << ", flip " << outsideFlipped
              << ")";
        }
      }
    }
  }
  return test;
}

template <class Basis>
Dune::TestSuite testTransformedDeltaProperty(Basis const& basis)
{
  Dune::TestSuite test("Transformed Arnold-Winther delta property");
  auto localView = basis.localView();
  for (const auto& element : elements(basis.gridView())) {
    localView.bind(element);
    test.subTest(testDeltaProperty(localView.tree().finiteElement().localBasis(),
                                   localView.tree().finiteElement().localInterpolation()));
  }
  return test;
}

template <class Basis>
Dune::TestSuite testEdgeInterpolationConsistency(Basis const& basis)
{
  Dune::TestSuite test("Arnold-Winther edge interpolation consistency");
  auto insideLocalView = basis.localView();
  auto outsideLocalView = basis.localView();
  const auto& indexSet = basis.gridView().indexSet();

  auto constantTensor = [](const auto&) {
    return FieldMatrix<double, 2, 2>{{1.0, 0.25}, {0.25, 2.0}};
  };

  for (const auto& inside : elements(basis.gridView())) {
    insideLocalView.bind(inside);
    for (const auto& intersection : intersections(basis.gridView(), inside)) {
      if (not intersection.neighbor())
        continue;
      const auto outside = intersection.outside();
      if (indexSet.index(inside) >= indexSet.index(outside))
        continue;

      outsideLocalView.bind(outside);
      std::vector<double> insideCoefficients;
      std::vector<double> outsideCoefficients;
      insideLocalView.tree().finiteElement().localInterpolation().interpolate(constantTensor,
                                                                              insideCoefficients);
      outsideLocalView.tree().finiteElement().localInterpolation().interpolate(constantTensor,
                                                                               outsideCoefficients);

      const auto insideEdge = intersection.indexInInside();
      const auto outsideEdge = intersection.indexInOutside();
      for (std::size_t edgeDOF = 0; edgeDOF < 4; ++edgeDOF) {
        const auto insideDOF = 9 + 4 * insideEdge + edgeDOF;
        const auto outsideDOF = 9 + 4 * outsideEdge + edgeDOF;
        test.check(std::abs(insideCoefficients[insideDOF] - outsideCoefficients[outsideDOF]) <
                   1e-12)
            << "Edge interpolation differs for local edge DOF " << edgeDOF;
      }
    }
  }
  return test;
}

Dune::TestSuite testNonAffineGeometryRejection()
{
  Dune::TestSuite test("Non-affine Arnold-Winther geometry rejection");
  using FiniteElement =
      Dune::Functions::Impl::ArnoldWintherLocalFiniteElement<NonAffineElementStub, double, double>;
  FiniteElement finiteElement;
  bool rejected = false;
  try {
    finiteElement.bind({}, NonAffineElementStub{});
  } catch (const Dune::NotImplemented&) {
    rejected = true;
  }
  test.check(rejected) << "A non-affine geometry must not use the affine divergence formula";
  return test;
}

template <class GridView>
Dune::TestSuite testEmbeddedSurfaceRejection(const GridView& gridView)
{
  Dune::TestSuite test("Embedded Arnold-Winther surface rejection");
  bool rejected = false;
  try {
    [[maybe_unused]] auto basis = makeBasis(gridView, BasisFactory::arnoldWinther());
  } catch (const Dune::NotImplemented&) {
    rejected = true;
  }
  test.check(rejected) << "Arnold-Winther must reject dimension != dimensionworld when the "
                          "basis is constructed";
  return test;
}

int main(int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test("arnold-winther");
  std::cout << "Testing AW reference finite element" << std::endl;
  // first test the plain reference basis and interpolation
  using Basis = Dune::Functions::Impl::ArnoldWintherReferenceLocalBasis<double, double>;
  Basis basis;
  Dune::Functions::Impl::ArnoldWintherReferenceLocalInterpolation<double, double> interpolation;
  test.subTest(testDeltaProperty(basis, interpolation));

  std::cout << "Testing AW finite element on grid with one element" << std::endl;
  // Second test with transformed basis and global interpolation
  using namespace Dune::Functions::BasisFactory;
  using Grid = UGGrid<2>;
  // test on a Grid with one triangle
  {
    auto gridFactory = GridFactory<Grid>();
    gridFactory.insertVertex({0., 0.});
    gridFactory.insertVertex({1., 0.});
    gridFactory.insertVertex({0., 1.});

    gridFactory.insertElement(GeometryTypes::simplex(2), {0, 2, 1});

    auto grid = gridFactory.createGrid();
    auto gridView = grid->leafGridView();

    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, arnoldWinther());

      test.subTest(testEdgeDOFNumbering(basis));
      test.subTest(checkBasis(basis, EnableNormal_VectorContinuityCheck()));
    }
  }

  std::cout << "Testing AW finite element on grid with two elements" << std::endl;
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
      test.subTest(testTransformedDeltaProperty(basis));
      test.subTest(testEdgeInterpolationConsistency(basis));
      test.subTest(testEdgeDOFNumbering(basis));
      test.subTest(checkBasis(basis, EnableNormal_VectorContinuityCheck()));
    }
  }

  std::cout << "Testing rejection of an embedded surface" << std::endl;
  {
    auto gridFactory = GridFactory<Grid>();
    gridFactory.insertVertex({0., 0.});
    gridFactory.insertVertex({1.1, 0.});
    gridFactory.insertVertex({1.1, 1.});
    gridFactory.insertVertex({0., 1.});
    gridFactory.insertElement(GeometryTypes::simplex(2), {0, 1, 2});
    gridFactory.insertElement(GeometryTypes::simplex(2), {2, 0, 3});
    std::shared_ptr<Grid> hostGrid = gridFactory.createGrid();
    auto coordinates = std::make_shared<AffineSurfaceCoordinates>();
    using SurfaceGrid = GeometryGrid<Grid, AffineSurfaceCoordinates>;
    auto surfaceGrid = std::make_shared<SurfaceGrid>(hostGrid, coordinates);
    test.subTest(testEmbeddedSurfaceRejection(surfaceGrid->leafGridView()));
  }

  test.subTest(testNonAffineGeometryRejection());
  return test.exit();
}
