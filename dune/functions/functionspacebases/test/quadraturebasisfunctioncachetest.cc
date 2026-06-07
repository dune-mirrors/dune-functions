// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <dune/common/test/testsuite.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/quadraturebasisfunctioncache.hh>
#include <dune/functions/functionspacebases/raviartthomasbasis.hh>
#include <dune/functions/functionspacebases/transformed/derivative.hh>

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
    typename LocalView::Tree,Functions::Derivatives::Value,Functions::Derivatives::Jacobian> cache;

  auto const& firstRule = QuadratureRules<double,1>::rule(GeometryTypes::cube(1),1);
  auto const& secondRule = QuadratureRules<double,1>::rule(GeometryTypes::cube(1),3);

  {
    auto localView = firstBasis.localView();
    localView.bind(*elements(firstBasis.gridView()).begin());
    cache.bind(localView);
    auto const& values = cache.get().evaluate(Functions::Derivatives::Value{},firstRule);
    test.check(std::abs(values[0][0][0]-(1-firstRule[0].position()[0])) < 1e-14);
  }

  auto localView = secondBasis.localView();
  localView.bind(*elements(secondBasis.gridView()).begin());
  cache.bind(localView);

  auto const& values = cache.get().evaluate(Functions::Derivatives::Value{},secondRule);
  test.check(values.size() == secondRule.size() && values[0].size() == 2,
    "value cache must contain all basis functions");
  if (values.size() == secondRule.size() && values[0].size() == 2)
    test.check(std::abs(values[0][0][0]-(1-secondRule[0].position()[0])) < 1e-14,
      "quadrature rules with different keys must use distinct precomputations");

  auto const& jacobians = cache.get().evaluate(Functions::Derivatives::Jacobian{},secondRule);
  test.check(jacobians.size() == secondRule.size() && jacobians[0].size() == 2,
    "derivative caches must be initialized independently");
  if (jacobians.size() == secondRule.size() && jacobians[0].size() == 2)
    test.check(std::abs(jacobians[0][0][0][0]+0.5) < 1e-14,
      "cache must finalize with the currently bound local view");

  auto const& firstRuleValues = cache.get().evaluate(Functions::Derivatives::Value{},firstRule);
  test.check(firstRuleValues.size() == firstRule.size() && firstRuleValues[0].size() == 2);
  if (firstRuleValues.size() == firstRule.size() && firstRuleValues[0].size() == 2)
    test.check(std::abs(firstRuleValues[0][0][0]-(1-firstRule[0].position()[0])) < 1e-14,
      "precomputed buffers for different quadrature rules must coexist");


  // test mixed finite elements
  auto triGrid = Dune::StructuredGridFactory<UGGrid<2>>::createSimplexGrid({0.0,0.0},{1.0,1.0},{2,2});
  auto mixedBasis = makeBasis(triGrid->leafGridView(),composite(
    lagrange<1>(), raviartThomas<1>(), blockedLexicographic()));
  using MixedLocalView = decltype(mixedBasis)::LocalView;
  QuadratureBasisFunctionCache<
    typename MixedLocalView::Tree,Functions::Derivatives::Value,Functions::Derivatives::Jacobian,
    Functions::Derivatives::Divergence> mixedCache;

  {
    using namespace Dune::Indices;
    auto localView = mixedBasis.localView();
    localView.bind(*elements(mixedBasis.gridView()).begin());
    mixedCache.bind(localView);
    auto const& mixedRule = QuadratureRules<double,2>::rule(GeometryTypes::simplex(2),3);
    auto const& values = mixedCache.get(_0).evaluate(Functions::Derivatives::Value{},mixedRule);
    auto const& jacobians = mixedCache.get(_0).evaluate(Functions::Derivatives::Jacobian{},mixedRule);
    auto const& divs = mixedCache.get(_1).evaluate(Functions::Derivatives::Divergence{},mixedRule);
  }

  return test.exit();
}
