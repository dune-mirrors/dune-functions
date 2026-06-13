// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#include <doc/grids/gridfactory/hybridtestgrids.hh>

#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/quadraturebasisfunctioncache.hh>
#include <dune/functions/functionspacebases/raviartthomasbasis.hh>
#include <dune/functions/functionspacebases/transformed/derivative.hh>

namespace {

using namespace Dune;
using namespace Dune::Functions;
using namespace Dune::Functions::BasisFactory;

/**
 * \brief Check rule-key coexistence and rebinding after local-view destruction.
 */
void testRuleCachingAndRebinding(TestSuite& test)
{
  YaspGrid<1> firstGrid({1.0},{1});
  YaspGrid<1> secondGrid({2.0},{1});
  auto firstBasis = makeBasis(firstGrid.leafGridView(),lagrange<1>());
  auto secondBasis = makeBasis(secondGrid.leafGridView(),lagrange<1>());
  using LocalView = decltype(firstBasis)::LocalView;
  QuadratureBasisFunctionCache<
    typename LocalView::Tree,
    Functions::Derivatives::Value,Functions::Derivatives::Jacobian> cache;

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

  cache.clear();
  cache.bind(localView);
  auto const& valuesAfterClear =
    cache.get().evaluate(Functions::Derivatives::Value{},secondRule);
  test.check(valuesAfterClear.size() == secondRule.size());
  if (valuesAfterClear.size() == secondRule.size())
    test.check(std::abs(valuesAfterClear[0][0][0]
      -(1-secondRule[0].position()[0])) < 1e-14,
      "cache must remain usable after explicit invalidation");

  test.check(cache.get().size() == 2);
  test.check(cache.get().order() == 1);
  test.check(cache.get().basis().size() == cache.get().size());
}

/**
 * \brief Check that one tree cache supports simplex and cube elements.
 */
void testHybridGeometryTypes(TestSuite& test)
{
  auto grid = make2DHybridTestGrid<UGGrid<2>>();
  auto basis = makeBasis(grid->leafGridView(),lagrange<1>());
  auto localView = basis.localView();
  QuadratureBasisFunctionCache<
    typename decltype(localView)::Tree,
    Functions::Derivatives::Value,Functions::Derivatives::Jacobian> cache;

  bool testedSimplex = false;
  bool testedCube = false;
  for (auto const& element : elements(basis.gridView())) {
    localView.bind(element);
    cache.bind(localView);
    auto const& rule = QuadratureRules<double,2>::rule(element.type(),2);
    auto const& values = cache.get().evaluate(Functions::Derivatives::Value{},rule);
    auto const& jacobians = cache.get().evaluate(Functions::Derivatives::Jacobian{},rule);
    test.check(values.size() == rule.size());
    test.check(jacobians.size() == rule.size());
    test.check(values[0].size() == cache.get().size());
    test.check(jacobians[0].size() == cache.get().size());
    testedSimplex = testedSimplex || element.type().isSimplex();
    testedCube = testedCube || element.type().isCube();
  }
  test.check(testedSimplex && testedCube,
    "hybrid cache test must visit simplex and cube elements");
}

/**
 * \brief Check heterogeneous derivative tags and result range types.
 */
void testMixedDerivativeRanges(TestSuite& test)
{
  auto grid = StructuredGridFactory<UGGrid<2>>::createSimplexGrid(
    {0.0,0.0},{1.0,1.0},{2,2});
  auto basis = makeBasis(grid->leafGridView(),composite(
    lagrange<1>(),raviartThomas<1>(),blockedLexicographic()));
  auto localView = basis.localView();
  QuadratureBasisFunctionCache<
    typename decltype(localView)::Tree,
    Functions::Derivatives::Value,
    Functions::Derivatives::Jacobian,
    Functions::Derivatives::Divergence> cache;

  localView.bind(*elements(basis.gridView()).begin());
  cache.bind(localView);
  auto const& rule = QuadratureRules<double,2>::rule(GeometryTypes::simplex(2),3);

  using namespace Dune::Indices;
  auto& scalarCache = cache.get(_0);
  auto& piolaCache = cache.get(_1);
  auto const& scalarValues = scalarCache.evaluate(Functions::Derivatives::Value{},rule);
  auto const& scalarJacobians = scalarCache.evaluate(Functions::Derivatives::Jacobian{},rule);
  auto const& piolaValues = piolaCache.evaluate(Functions::Derivatives::Value{},rule);
  auto const& divergences = piolaCache.evaluate(Functions::Derivatives::Divergence{},rule);

  static_assert(std::same_as<
    std::remove_cvref_t<decltype(scalarValues[0][0])>,FieldVector<double,1>>);
  static_assert(std::same_as<
    std::remove_cvref_t<decltype(scalarJacobians[0][0])>,FieldMatrix<double,1,2>>);
  static_assert(std::same_as<
    std::remove_cvref_t<decltype(piolaValues[0][0])>,FieldVector<double,2>>);
  static_assert(std::same_as<
    std::remove_cvref_t<decltype(divergences[0][0])>,double>);

  test.check(scalarValues[0].size() == scalarCache.size());
  test.check(scalarJacobians[0].size() == scalarCache.size());
  test.check(piolaValues[0].size() == piolaCache.size());
  test.check(divergences[0].size() == piolaCache.size());
}

/**
 * \brief Check that Piola precomputations distinguish orientation variants.
 */
void testPiolaOrientationVariants(TestSuite& test)
{
  YaspGrid<2> grid({1.0,1.0},{2,2});
  auto basis = makeBasis(grid.leafGridView(),raviartThomas<1>());
  auto localView = basis.localView();
  QuadratureBasisFunctionCache<
    typename decltype(localView)::Tree,
    Functions::Derivatives::Value,
    Functions::Derivatives::Divergence> cache;
  auto const& rule =
    QuadratureRules<double,2>::rule(GeometryTypes::cube(2),3);

  for (auto const& element : elements(basis.gridView())) {
    localView.bind(element);
    cache.bind(localView);
    auto& leafCache = cache.get();
    auto const& cachedValues =
      leafCache.evaluate(Functions::Derivatives::Value{},rule);
    auto const& cachedDivergences =
      leafCache.evaluate(Functions::Derivatives::Divergence{},rule);

    using PhysicalBasis =
      typename decltype(localView)::Tree::FiniteElement::PhysicalBasis;
    std::vector<
      typename PhysicalBasis::template DerivativeRange<
        Functions::Derivatives::Value>> values;
    std::vector<
      typename PhysicalBasis::template DerivativeRange<
        Functions::Derivatives::Divergence>> divergences;

    for (std::size_t iq = 0; iq < rule.size(); ++iq) {
      auto const& position = rule[iq].position();
      leafCache.basis().evaluate(
        Functions::Derivatives::Value{},position,values);
      leafCache.basis().evaluate(
        Functions::Derivatives::Divergence{},position,divergences);

      for (std::size_t i = 0; i < values.size(); ++i) {
        test.check((values[i]-cachedValues[iq][i]).two_norm() < 1e-14,
          "cached Piola values must use the bound orientation variant");
        test.check(std::abs(divergences[i]-cachedDivergences[iq][i]) < 1e-14,
          "cached Piola divergences must use the bound orientation variant");
      }
    }
  }
}

} // namespace

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc,argv);
  Dune::TestSuite test;

  testRuleCachingAndRebinding(test);
  testHybridGeometryTypes(test);
  testMixedDerivativeRanges(test);
  testPiolaOrientationVariants(test);

  return test.exit();
}
