// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/functions/functionspacebases/hellanhermannjohnsonbasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;

template<int k, class GridView>
void testHellanHermannJohnsonBasis(TestSuite& test, const GridView& gridView)
{
  std::cout<<"  Testing order: "<< k <<std::endl;

  // Check basis created 'manually'
  {
    Functions::HellanHermannJohnsonBasis<GridView,k> basis(gridView);
    test.subTest(checkBasis(basis, EnableNormalNormalContinuityCheck()));
  }

  // Check basis created using basis builder mechanism
  {
    using namespace Functions::BasisFactory;
    auto basis = makeBasis(gridView, hhj<k>());
    test.subTest(checkBasis(basis, EnableNormalNormalContinuityCheck()));
  }
}


int main (int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test;

  // Test with pure simplex grid
  // (Unfortunately there is no grid implementation available that only supports simplices.)
  std::cout<<"Testing Hellan-Hermann-Johnson basis in 2D with simplex grid\n";
  auto triangleGrid = Dune::StructuredGridFactory<UGGrid<2>>::createSimplexGrid({0.0,0.0},{1.0,1.0},{4u,4u});
  auto triangleGridView = triangleGrid->leafGridView();
  testHellanHermannJohnsonBasis<0>(test, triangleGridView);
  testHellanHermannJohnsonBasis<1>(test, triangleGridView);
  testHellanHermannJohnsonBasis<2>(test, triangleGridView);

  return test.exit();
}
