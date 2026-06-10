// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include "basistest.hh"
#include <config.h>

#include <cmath>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/functions/functionspacebases/hellanherrmannjohnsonbasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/functionspacebases/test/testtransformedlocalbasis.hh>


using namespace Dune;

template<class LocalFiniteElement>
TestSuite testInterpolationDuality(LocalFiniteElement const& finiteElement)
{
  TestSuite test("HHJ interpolation duality");

  auto const& basis = finiteElement.physicalBasis();
  for (std::size_t shapeFunction = 0; shapeFunction < basis.size(); ++shapeFunction) {
    auto function = [&](auto const& x) {
      std::vector<typename LocalFiniteElement::PhysicalBasis::Range> values;
      basis.evaluate(Functions::Derivatives::Value{},x,values);
      return values[shapeFunction];
    };

    std::vector<double> coefficients;
    finiteElement.localInterpolation().interpolate(function,coefficients);
    test.check(coefficients.size() == basis.size(), "interpolation result size");
    for (std::size_t i = 0; i < coefficients.size(); ++i)
      test.check(std::abs(coefficients[i] - (i == shapeFunction ? 1.0 : 0.0)) < 1e-8,
        "interpolation duality");
  }

  return test;
}

template<int k, class GridView>
void testHellanHerrmannJohnsonBasis(TestSuite& test, const GridView& gridView)
{
  std::cout<<"  Testing order: "<< k <<std::endl;
  // Check basis created 'manually'
  {
    Functions::HellanHerrmannJohnsonBasis<GridView,k> basis(gridView);
    test.subTest(checkBasis(basis, EnableNormalNormalContinuityCheck()));

    auto localView = basis.localView();
    localView.bind(*gridView.template begin<0>());
    test.subTest(Functions::Test::testTransformedLocalFiniteElement(
      localView.tree().finiteElement(),
      localView.element(),
      Functions::Derivatives::Value{},
      Functions::Derivatives::DivDiv{}));
    test.subTest(testInterpolationDuality(localView.tree().finiteElement()));
  }

  // Check basis created using basis builder mechanism
  // {
  //   using namespace Functions::BasisFactory;
  //   auto basis = makeBasis(gridView, hhj<k>());
  //   test.subTest(checkBasis(basis, EnableNormalNormalContinuityCheck()));
  // }
}

TestSuite test2d()
{
  TestSuite test("HHJ_2d");

  // Test with pure simplex grid
  std::cout<<"Testing Hellan-Herrmann-Johnson basis in 2D with simplex grid\n";
  {
    auto triangleGrid = Dune::StructuredGridFactory<UGGrid<2>>::createSimplexGrid({0.0,0.0},{1.0,1.0},{1u,1u});
    auto triangleGridView = triangleGrid->leafGridView();
    testHellanHerrmannJohnsonBasis<0>(test, triangleGridView);
    testHellanHerrmannJohnsonBasis<1>(test, triangleGridView);
    testHellanHerrmannJohnsonBasis<2>(test, triangleGridView);
  }

  {
    auto gridFactory = GridFactory<UGGrid<2>>();
    gridFactory.insertVertex({0., 0.});
    gridFactory.insertVertex({1., 0.});
    gridFactory.insertVertex({1., 1.});
    gridFactory.insertVertex({0., 1.});

    gridFactory.insertElement(GeometryTypes::simplex(2), {0, 1, 2});
    gridFactory.insertElement(GeometryTypes::simplex(2), {2, 3,0});

    auto grid = gridFactory.createGrid();
    auto gridView = grid->leafGridView();

    testHellanHerrmannJohnsonBasis<0>(test, gridView);
    testHellanHerrmannJohnsonBasis<1>(test, gridView);
    testHellanHerrmannJohnsonBasis<2>(test, gridView);
    testHellanHerrmannJohnsonBasis<3>(test, gridView);
    // testHellanHerrmannJohnsonBasis<4>(test, gridView);
    // testHellanHerrmannJohnsonBasis<5>(test, gridView);
    // testHellanHerrmannJohnsonBasis<6>(test, gridView);

  }

  return test;
}

int main (int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test;

  test.subTest(test2d());
  // test.subTest(test3d());

  return test.exit();
}
