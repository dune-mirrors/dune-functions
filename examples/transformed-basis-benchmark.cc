// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/quadraturebasisfunctioncache.hh>
#include <dune/functions/functionspacebases/quadraturebasisfunctioncache2.hh>
#include <dune/functions/functionspacebases/raviartthomasbasis.hh>
#include <dune/functions/functionspacebases/transformed/derivative.hh>

namespace {

using namespace Dune;
using namespace Dune::Functions;
namespace FDerivatives = Dune::Functions::Derivatives;

struct Measurement
{
  double seconds;
  double checksum;
};

/**
 * \brief Add one quadrature-point contribution to an H(div) local matrix.
 *
 * The assembled form is the sum of the vector-valued mass product and the
 * divergence-divergence product.
 */
template<class Values, class Divergences>
void addQuadraturePoint(std::vector<double>& matrix,
                        Values const& values,
                        Divergences const& divergences,
                        double factor)
{
  auto const size = values.size();
  for (std::size_t i = 0; i < size; ++i)
    for (std::size_t j = 0; j < size; ++j)
      matrix[i*size+j] +=
        (values[i]*values[j] + divergences[i]*divergences[j])*factor;
}

/**
 * \brief Assemble with a callback evaluating one quadrature point at a time.
 */
template<class Basis, class LocalView, class QuadratureRule,
         class Values, class Divergences, class Evaluator>
double assemblePointwise(Basis const& basis,
                         LocalView& localView,
                         QuadratureRule const& quadrature,
                         Values& values,
                         Divergences& divergences,
                         Evaluator&& evaluate)
{
  std::vector<double> localMatrix;
  double checksum = 0;

  for (auto const& element : elements(basis.gridView())) {
    localView.bind(element);
    auto const geometry = element.geometry();
    auto const size = localView.tree().finiteElement().size();
    localMatrix.assign(size*size,0.0);

    for (std::size_t iq = 0; iq < quadrature.size(); ++iq) {
      auto const& position = quadrature[iq].position();
      evaluate(localView.tree().finiteElement(),position,values,divergences);
      auto const factor =
        quadrature[iq].weight()*geometry.integrationElement(position);
      addQuadraturePoint(localMatrix,values,divergences,factor);
    }

    checksum += std::accumulate(localMatrix.begin(),localMatrix.end(),0.0);
  }

  return checksum;
}

/**
 * \brief Assemble using cached reference evaluations and bound finalization.
 */
template<class Basis, class LocalView, class Cache, class QuadratureRule>
double assembleCached(Basis const& basis,
                      LocalView& localView,
                      Cache& cache,
                      QuadratureRule const& quadrature)
{
  std::vector<double> localMatrix;
  double checksum = 0;

  for (auto const& element : elements(basis.gridView())) {
    localView.bind(element);
    cache.bind(localView);
    auto const geometry = element.geometry();
    auto& leafCache = cache.get();
    auto const& values =
      leafCache.evaluate(FDerivatives::Value{},quadrature);
    auto const& divergences =
      leafCache.evaluate(FDerivatives::Divergence{},quadrature);
    localMatrix.assign(leafCache.size()*leafCache.size(),0.0);

    for (std::size_t iq = 0; iq < quadrature.size(); ++iq) {
      auto const& position = quadrature[iq].position();
      auto const factor =
        quadrature[iq].weight()*geometry.integrationElement(position);
      addQuadraturePoint(localMatrix,values[iq],divergences[iq],factor);
    }

    checksum += std::accumulate(localMatrix.begin(),localMatrix.end(),0.0);
  }

  return checksum;
}

/**
 * \brief Measure repeated assembly after one untimed warm-up pass.
 */
template<class F>
Measurement measure(std::size_t repetitions, F&& assemble)
{
  assemble();
  Timer timer;
  double checksum = 0;
  for (std::size_t repetition = 0; repetition < repetitions; ++repetition)
    checksum += assemble();
  return {timer.elapsed()/repetitions,checksum/repetitions};
}

/**
 * \brief Verify that a benchmark path assembles the same matrix entries.
 */
void checkChecksum(std::string const& name,
                   double checksum,
                   double reference)
{
  auto const tolerance =
    1e-11*std::max({1.0,std::abs(checksum),std::abs(reference)});
  if (std::abs(checksum-reference) > tolerance)
    DUNE_THROW(Exception,
      name << " checksum " << checksum
           << " differs from manual checksum " << reference);
}

void printMeasurement(std::string const& name,
                      Measurement const& measurement,
                      double referenceTime)
{
  std::cout << std::left << std::setw(30) << name
            << std::right << std::fixed << std::setprecision(3)
            << std::setw(10) << 1e3*measurement.seconds << " ms"
            << std::setw(10) << referenceTime/measurement.seconds << " x"
            << "  checksum " << std::scientific << std::setprecision(12)
            << measurement.checksum << std::defaultfloat << '\n';
}

} // namespace

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc,argv);

  try {
    std::size_t cellsPerDirection = 32;
    std::size_t repetitions = 5;
    int quadratureOrder = 4;
    if (argc > 1)
      cellsPerDirection = std::stoul(argv[1]);
    if (argc > 2)
      repetitions = std::stoul(argv[2]);
    if (argc > 3)
      quadratureOrder = std::stoi(argv[3]);

    constexpr int dim = 2;
    FieldVector<double,dim> upperRight(1.0);
    std::array<int,dim> cells;
    cells.fill(cellsPerDirection);
    YaspGrid<dim> grid(upperRight,cells);

    using namespace Dune::Functions::BasisFactory;
    auto basis = makeBasis(grid.leafGridView(),raviartThomas<1>());
    auto const& quadrature =
      QuadratureRules<double,dim>::rule(
        GeometryTypes::cube(dim),quadratureOrder);

    using LocalView = typename decltype(basis)::LocalView;
    using FiniteElement = typename LocalView::Tree::FiniteElement;
    using ReferenceBasis = typename FiniteElement::ReferenceBasis;
    using PhysicalBasis = typename FiniteElement::PhysicalBasis;
    using Value =
      typename PhysicalBasis::template DerivativeRange<FDerivatives::Value>;
    using Divergence =
      typename PhysicalBasis::template DerivativeRange<FDerivatives::Divergence>;

    std::vector<Value> values;
    std::vector<Divergence> divergences;

    auto manualLocalView = basis.localView();
    std::vector<typename ReferenceBasis::Traits::RangeType> referenceValues;
    std::vector<typename ReferenceBasis::Traits::JacobianType> referenceJacobians;
    auto manualAssembly = [&] {
      return assemblePointwise(
        basis,manualLocalView,quadrature,values,divergences,
        [&](auto const& finiteElement, auto const& position,
            auto& physicalValues, auto& physicalDivergences) {
          auto const& referenceBasis = finiteElement.referenceBasis();
          referenceBasis.evaluateFunction(position,referenceValues);
          referenceBasis.evaluateJacobian(position,referenceJacobians);

          physicalValues.resize(referenceValues.size());
          physicalDivergences.resize(referenceJacobians.size());
          auto const geometry = manualLocalView.element().geometry();
          auto const jacobian = geometry.jacobian(position);
          auto const integrationElement = geometry.integrationElement(position);
          for (std::size_t i = 0; i < referenceValues.size(); ++i) {
            jacobian.mv(referenceValues[i],physicalValues[i]);
            physicalValues[i] /= integrationElement;
            physicalDivergences[i] = 0;
            for (std::size_t j = 0; j < dim; ++j)
              physicalDivergences[i] += referenceJacobians[i][j][j];
            physicalDivergences[i] /= integrationElement;
          }
        });
    };

    auto evaluateLocalView = basis.localView();
    auto evaluateAssembly = [&] {
      return assemblePointwise(
        basis,evaluateLocalView,quadrature,values,divergences,
        [](auto const& finiteElement, auto const& position,
          auto& physicalValues, auto& physicalDivergences) {
          finiteElement.physicalBasis().evaluate(
            FDerivatives::Value{},position,physicalValues);
          finiteElement.physicalBasis().evaluate(
            FDerivatives::Divergence{},position,physicalDivergences);
        });
    };

    auto splitLocalView = basis.localView();
    typename PhysicalBasis::template PrecomputeBuffer<FDerivatives::Value>
      valuePrecomputeBuffer;
    typename PhysicalBasis::template PrecomputeBuffer<FDerivatives::Divergence>
      divergencePrecomputeBuffer;
    auto splitAssembly = [&] {
      return assemblePointwise(
        basis,splitLocalView,quadrature,values,divergences,
        [&](auto const& finiteElement, auto const& position,
            auto& physicalValues, auto& physicalDivergences) {
          auto const& physicalBasis = finiteElement.physicalBasis();
          physicalBasis.precompute(
            FDerivatives::Value{},position,valuePrecomputeBuffer);
          physicalBasis.finalize(
            FDerivatives::Value{},position,valuePrecomputeBuffer,physicalValues);
          physicalBasis.precompute(
            FDerivatives::Divergence{},position,divergencePrecomputeBuffer);
          physicalBasis.finalize(
            FDerivatives::Divergence{},position,
            divergencePrecomputeBuffer,physicalDivergences);
        });
    };

    auto cacheLocalView = basis.localView();
    QuadratureBasisFunctionCache<
      typename LocalView::Tree,
      FDerivatives::Value,
      FDerivatives::Divergence> cache;
    auto cachedAssembly = [&] {
      return assembleCached(
        basis,cacheLocalView,cache,quadrature);
    };

    QuadratureBasisFunctionCache2<
      typename LocalView::Tree,
      FDerivatives::Value,
      FDerivatives::Divergence> cache2;
    auto cachedAssembly2 = [&] {
      return assembleCached(
        basis,cacheLocalView,cache2,quadrature);
    };

    auto const manual = measure(repetitions,manualAssembly);
    auto const evaluate = measure(repetitions,evaluateAssembly);
    auto const split = measure(repetitions,splitAssembly);
    auto const cached = measure(repetitions,cachedAssembly);
    auto const cached2 = measure(repetitions,cachedAssembly2);

    checkChecksum("evaluate",evaluate.checksum,manual.checksum);
    checkChecksum("precompute/finalize",split.checksum,manual.checksum);
    checkChecksum("quadrature cache",cached.checksum,manual.checksum);
    checkChecksum("quadrature cache2",cached2.checksum,manual.checksum);

    std::cout << "RT1 H(div) mass + div-div assembly\n"
              << "grid: " << cellsPerDirection << " x "
              << cellsPerDirection << ", elements: "
              << basis.gridView().size(0)
              << ", quadrature points: " << quadrature.size()
              << ", repetitions: " << repetitions << "\n\n";
#ifndef NDEBUG
    std::cout << "Warning: assertions are enabled; use a release build for "
                 "representative timings.\n\n";
#endif
    std::cout << std::left << std::setw(30) << "path"
              << std::right << std::setw(15) << "time/pass"
              << std::setw(14) << "relative"
              << '\n';
    printMeasurement("manual reference + Piola",manual,manual.seconds);
    printMeasurement("evaluate()",evaluate,manual.seconds);
    printMeasurement("precompute() + finalize()",split,manual.seconds);
    printMeasurement("quadrature cache",cached,manual.seconds);
    printMeasurement("quadrature cache2",cached2,manual.seconds);
  }
  catch (Dune::Exception const& exception) {
    std::cerr << exception << '\n';
    return 1;
  }
  catch (std::exception const& exception) {
    std::cerr << exception.what() << '\n';
    return 1;
  }
}
