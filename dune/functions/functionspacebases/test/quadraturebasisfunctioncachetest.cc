// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <dune/common/test/testsuite.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/quadraturebasisfunctioncache.hh>

int main(int argc, char** argv)
{
  using namespace Dune;
  using namespace Dune::Functions;
  using namespace Dune::Functions::BasisFactory;

  MPIHelper::instance(argc,argv);
  TestSuite test;
  YaspGrid<1> firstGrid({1.0},{1});
  YaspGrid<1> secondGrid({2.0},{1});
  auto firstBasis = makeBasis(firstGrid.leafGridView(),lagrange<1>());
  auto secondBasis = makeBasis(secondGrid.leafGridView(),lagrange<1>());
  using LocalView = decltype(firstBasis)::LocalView;
  QuadratureBasisFunctionCache<
    typename LocalView::Tree,Derivatives::Value,Derivatives::Jacobian> cache;

  QuadratureRule<double,1> firstRule;
  firstRule.emplace_back(FieldVector<double,1>{0.25},1.0);
  QuadratureRule<double,1> secondRule;
  secondRule.emplace_back(FieldVector<double,1>{0.75},1.0);

  {
    auto localView = firstBasis.localView();
    localView.bind(*elements(firstBasis.gridView()).begin());
    cache.bind(localView);
    auto const& values = cache.get().evaluate(Derivatives::Value{},firstRule);
    test.check(std::abs(values[0][0][0]-0.75) < 1e-14);
  }

  auto localView = secondBasis.localView();
  localView.bind(*elements(secondBasis.gridView()).begin());
  cache.bind(localView);

  auto const& values = cache.get().evaluate(Derivatives::Value{},secondRule);
  test.check(std::abs(values[0][0][0]-0.25) < 1e-14,
    "equal-sized quadrature rules must not share precomputed values");

  auto const& jacobians = cache.get().evaluate(Derivatives::Jacobian{},secondRule);
  test.check(std::abs(jacobians[0][0][0][0]+0.5) < 1e-14,
    "cache must finalize with the currently bound local view");

  return test.exit();
}
