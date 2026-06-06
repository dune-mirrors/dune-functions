// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/functions/functionspacebases/nedelecbasis.hh>
#include <dune/functions/functionspacebases/raviartthomasbasis.hh>
#include <dune/functions/functionspacebases/test/testtransformedlocalbasis.hh>
#include <dune/functions/functionspacebases/transformed/derivative.hh>


using namespace Dune;


template<class Basis, class Derivative>
TestSuite checkBasisFEs(const Basis& basis, Derivative d) {
  TestSuite test("basis finite elements");
  auto localView = basis.localView();
  for (const auto& element : elements(basis.gridView()))
  {
    localView.bind(element);
    test.subTest(Functions::Test::testTransformedLocalFiniteElement(
      localView.tree().finiteElement(),
      localView.element(),
      Functions::Derivatives::Value{},
      d));

    auto const& finiteElement = localView.tree().finiteElement();
    auto const& physicalBasis = finiteElement.physicalBasis();
    using LegacyJacobian = typename std::decay_t<
      decltype(finiteElement.localBasis())>::Traits::JacobianType;
    std::vector<LegacyJacobian> legacyJacobians;
    finiteElement.localBasis().evaluateJacobian(
      Dune::FieldVector<double,2>{0.25,0.25},legacyJacobians);
    test.check(legacyJacobians.size() == finiteElement.size(),
      "legacy evaluateJacobian compatibility");

    for (std::size_t j = 0; j < physicalBasis.size(); ++j) {
      auto shapeFunction = [&](auto const& x) {
        std::vector<typename std::decay_t<decltype(physicalBasis)>::template DerivativeRange<
          Functions::Derivatives::Value>> values;
        physicalBasis.evaluate(Functions::Derivatives::Value{},x,values);
        return values[j];
      };
      std::vector<double> coefficients;
      finiteElement.localInterpolation().interpolate(shapeFunction,coefficients);
      for (std::size_t i = 0; i < coefficients.size(); ++i)
        test.check(std::abs(coefficients[i]-(i == j ? 1.0 : 0.0)) < 1e-10,
          "physical interpolation is dual to the transformed basis");
    }
  }
  return test;
}


int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);
  TestSuite test;

  using namespace Functions::BasisFactory;

  // Check with UGGrid
  {
    // Generate grid for testing
    const int dim = 2;
    using Grid = UGGrid<dim>;

    const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";

    std::string filename = path + "curved2d.msh";
    std::shared_ptr<Grid> grid = GmshReader<Grid>::read(filename);

    auto gridView = grid->leafGridView();

    ///////////////////////////////////////////////////////////////////////
    //  Test GlobalValuedLocalFiniteElement for a H(div)-conforming space
    //  We use the Raviart-Thomas basis.
    ///////////////////////////////////////////////////////////////////////

    test.subTest(checkBasisFEs(makeBasis(gridView, raviartThomas<0>()), Dune::Functions::Derivatives::Divergence{}));
    test.subTest(checkBasisFEs(makeBasis(gridView, raviartThomas<1>()), Dune::Functions::Derivatives::Divergence{}));

    ///////////////////////////////////////////////////////////////////////
    //  Test GlobalValuedLocalFiniteElement for a H(curl)-conforming space
    //  We use the Nedelec basis of the first kind.
    ///////////////////////////////////////////////////////////////////////

    test.subTest(checkBasisFEs(makeBasis(gridView, nedelec<1,1,double>()), Dune::Functions::Derivatives::Curl{}));
  }

  // Check with YaspGrid
  {
    // Generate grid for testing
    auto grid = YaspGrid<2>({1.0, 1.0}, {5,5});
    auto gridView = grid.leafGridView();

    ///////////////////////////////////////////////////////////////////////
    //  Test GlobalValuedLocalFiniteElement for a H(div)-conforming space
    //  We use the Raviart-Thomas basis.
    ///////////////////////////////////////////////////////////////////////

    test.subTest(checkBasisFEs(makeBasis(gridView, raviartThomas<0>()), Dune::Functions::Derivatives::Divergence{}));
    test.subTest(checkBasisFEs(makeBasis(gridView, raviartThomas<1>()), Dune::Functions::Derivatives::Divergence{}));
    test.subTest(checkBasisFEs(makeBasis(gridView, raviartThomas<2>()), Dune::Functions::Derivatives::Divergence{}));

    ///////////////////////////////////////////////////////////////////////
    //  Test GlobalValuedLocalFiniteElement for a H(curl)-conforming space
    //  We use the Nedelec basis of the first kind.
    ///////////////////////////////////////////////////////////////////////

    test.subTest(checkBasisFEs(makeBasis(gridView, nedelec<1,1,double>()), Dune::Functions::Derivatives::Curl{}));
  }

  return test.exit();

} catch (Exception &e)
{
  std::cerr << "Dune reported error: " << e.what() << std::endl;
  return 1;
}
