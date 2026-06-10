// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <array>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/vectorbasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;

namespace
{

template<class IndexMergingStrategy, class GridView>
Dune::TestSuite checkUniformVectorBasis(const GridView& gridView)
{
  Dune::TestSuite test;

  using SubPreBasis = Dune::Functions::LagrangePreBasis<GridView, 1>;
  using PreBasis = Dune::Functions::VectorPreBasis<IndexMergingStrategy, SubPreBasis>;
  using Basis = Dune::Functions::DefaultGlobalBasis<PreBasis>;

  std::vector<SubPreBasis> subPreBases;
  subPreBases.emplace_back(gridView);
  subPreBases.emplace_back(gridView);
  subPreBases.emplace_back(gridView);

  Basis basis(PreBasis{subPreBases});

  test.check(basis.preBasis().children() == 3, "vector basis children == 3");
  test.check(basis.preBasis().dimension() == 3 * basis.preBasis().subPreBasis(0).dimension(), "vector basis dimension");
  test.check(basis.preBasis().maxNodeSize() == 3 * basis.preBasis().subPreBasis(0).maxNodeSize(), "vector basis max node size");

  test.subTest(checkBasis(basis, EnableContinuityCheck()));

  return test;
}


template<class IndexMergingStrategy, class GridView>
Dune::TestSuite checkNonUniformVectorBasis(const GridView& gridView)
{
  Dune::TestSuite test;

  using SubPreBasis = Dune::Functions::LagrangePreBasis<GridView, -1>;
  using PreBasis = Dune::Functions::VectorPreBasis<IndexMergingStrategy, SubPreBasis>;
  using Basis = Dune::Functions::DefaultGlobalBasis<PreBasis>;

  std::vector<SubPreBasis> subPreBases;
  subPreBases.emplace_back(gridView, 1);
  subPreBases.emplace_back(gridView, 2);
  subPreBases.emplace_back(gridView, 3);

  Basis basis(PreBasis{subPreBases});

  test.check(basis.preBasis().children() == 3, "vector basis children == 3");
  test.check(basis.preBasis().dimension() ==
      basis.preBasis().subPreBasis(0).dimension()
    + basis.preBasis().subPreBasis(1).dimension()
    + basis.preBasis().subPreBasis(2).dimension(), "vector basis dimension");
  test.check(basis.preBasis().maxNodeSize() ==
      basis.preBasis().subPreBasis(0).maxNodeSize()
    + basis.preBasis().subPreBasis(1).maxNodeSize()
    + basis.preBasis().subPreBasis(2).maxNodeSize(), "vector basis max node size");

  test.subTest(checkBasis(basis, EnableContinuityCheck()));

  return test;
}

} // end anonymous namespace

int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;

  const int dim = 2;
  using Grid = Dune::YaspGrid<dim>;
  Dune::FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{4, 4}};
  Grid grid(l, elements);
  auto gridView = grid.leafGridView();

  using namespace Dune::Functions::BasisFactory;

  test.subTest(checkUniformVectorBasis<BlockedLexicographic>(gridView));
  test.subTest(checkUniformVectorBasis<FlatLexicographic>(gridView));

test.subTest(checkNonUniformVectorBasis<BlockedLexicographic>(gridView));
test.subTest(checkNonUniformVectorBasis<FlatLexicographic>(gridView));

  return test.exit();
}
