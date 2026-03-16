// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <iterator>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/printgrid.hh>
#include <dune/functions/functionspacebases/hellanhermannjohnsonbasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>


using namespace Dune;

template<int k, class GridView>
void testHellanHermannJohnsonBasis(TestSuite& test, const GridView& gridView)
{
  std::cout<<"  Testing order: "<< k <<std::endl;
  Dune::printGrid(gridView.grid(), Dune::MPIHelper::instance(), "grid");
  // Check basis created 'manually'
  {
    Functions::HellanHermannJohnsonBasis<GridView,k> basis(gridView);
    // test.subTest(checkBasis(basis, EnableNormalNormalContinuityCheck()));

    auto f = [](auto const& x) {
      return Dune::FieldMatrix<double,2,2>({
        {x[0], 2.0},
        {2.0, 2.0*x[1]}
      });
      // return Dune::FieldMatrix<double,2,2>({
      //   {2*x[0]*x[1],      x[0]*x[0] + x[1]},
      //   {x[0]*x[0] + x[1], (x[0]+x[1])*(x[0]+x[1])}
      // });
    };

    auto localView = basis.localView();
    auto const& indexSet = gridView.indexSet();
    for (auto const& e : elements(gridView))
    {
      std::cout << "Element[" << indexSet.index(e) << "]:" << std::endl;
      localView.bind(e);

      auto const& node = localView.tree();
      auto const& localFE = node.finiteElement();
      auto const& localIp = localFE.localInterpolation();

      auto local_f = [f,g=e.geometry()](auto const& x) { return f(g.global(x)); };

      std::vector<double> coeff;
      localIp.interpolate(local_f, coeff);
      for (std::size_t i = 0; i < coeff.size(); ++i)
        std::cout << "  c[" << localView.index(node.localIndex(i)) << "(" << i << ")] = " << coeff[i] << std::endl;
    }

  }

  // Check basis created using basis builder mechanism
  // {
  //   using namespace Functions::BasisFactory;
  //   auto basis = makeBasis(gridView, hhj<k>());
  //   test.subTest(checkBasis(basis, EnableNormalNormalContinuityCheck()));
  // }
}


int main (int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test;

  // Test with pure simplex grid
  // (Unfortunately there is no grid implementation available that only supports simplices.)
  std::cout<<"Testing Hellan-Hermann-Johnson basis in 2D with simplex grid\n";
  auto triangleGrid = Dune::StructuredGridFactory<UGGrid<2>>::createSimplexGrid({0.0,0.0},{1.0,1.0},{1u,1u});
  auto triangleGridView = triangleGrid->leafGridView();
  // testHellanHermannJohnsonBasis<0>(test, triangleGridView);
  testHellanHermannJohnsonBasis<1>(test, triangleGridView);
  // testHellanHermannJohnsonBasis<2>(test, triangleGridView);

  return test.exit();
}
