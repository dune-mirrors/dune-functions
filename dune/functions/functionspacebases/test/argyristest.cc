// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <dune/common/timer.hh>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/argyrisbasis.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/functionspacebases/test/enabledifferentiabilitycheck.hh>

using namespace Dune;
using namespace Dune::Functions;

// Hack: Disable test that has not been merged to master so far,
// by replacing it with a dummy.
template<int i=0>
class CheckLocalFiniteElementFlag {};

int main(int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test("argyris");

  using namespace Dune::Functions::BasisFactory;

  using Grid = UGGrid<2>;
  // using Grid = YaspGrid<2>;
  auto grid = StructuredGridFactory<Grid>::createSimplexGrid({0., 0.}, {1., 1.}, {{3, 3}});

  auto gridView = grid->leafGridView();
  // std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1)
  //           << " facettes and " << gridView.size(2) << " vertices" << std::endl;

  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, argyris());
  std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

  test.subTest(checkBasis(basis, EnableContinuityCheck(),
                          EnableDifferentiabilityCheck(),
                          CheckLocalFiniteElementFlag<2>()));

  auto lv = basis.localView();
  std::size_t repeat = 100000;
  Dune::Timer t;
  for (auto&& i : Dune::range(repeat))
    lv.bind(*std::begin(elements(gridView)));
  std::cout<<"Binding took "<<t.elapsed()<<"s\n";
  auto const& lb = lv.tree().finiteElement().localBasis();
  std::vector<Dune::FieldVector<double,1>> out;
  t.reset();
  for (auto&& i : Dune::range(repeat))
    lb.evaluateFunction({0.,0.},out);
  std::cout<<"Evaluations took "<<t.elapsed()<<"s\n";
  std::vector<Dune::FieldMatrix<double,1,2>> gradients;
  t.reset();
  for (auto&& i : Dune::range(repeat))
    lb.evaluateJacobian({0.,0.},gradients);
  std::cout<<"Jacobian Evaluations took "<<t.elapsed()<<"s\n";
  std::vector<Dune::FieldMatrix<double,2,2>> hessian;
  t.reset();
  for (auto&& i : Dune::range(repeat))
    lb.evaluateHessian({0.,0.},hessian);
  std::cout<<"Hessian Evaluations took "<<t.elapsed()<<"s\n";
  return test.exit();
}
