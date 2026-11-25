// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <memory>
#include <numeric>
#include <random>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/uggrid.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lobattobasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

using namespace Dune;
using namespace Dune::Functions;


// Create a structured grid of the unit cube with n^dim cells, where each cell is either a
// cube element or split into simplices (Kuhn triangulation), depending on `isSimplexCell`.
// To obtain non-trivial orientations of edges and faces, the vertices are inserted in random
// order and the vertex lists of the elements are transformed by random orientation-preserving
// symmetries of the reference element.
template<class Grid, class IsSimplexCell>
std::unique_ptr<Grid> createShuffledGrid (int n, IsSimplexCell&& isSimplexCell, unsigned int seed = 42)
{
  constexpr int dim = Grid::dimension;
  std::mt19937 rng(seed);

  // lattice vertices, inserted in random order
  const int nv = Dune::power(n+1, dim);
  std::vector<unsigned int> insertionIndex(nv);
  std::iota(insertionIndex.begin(), insertionIndex.end(), 0);
  std::shuffle(insertionIndex.begin(), insertionIndex.end(), rng);
  std::vector<FieldVector<double,dim>> positions(nv);
  for (int v = 0; v < nv; ++v)
    for (int d = 0, r = v; d < dim; ++d, r /= (n+1))
      positions[v][d] = double(r % (n+1)) / n;

  GridFactory<Grid> factory;
  std::vector<int> vertexAt(nv);
  for (int k = 0; k < nv; ++k)
    vertexAt[insertionIndex[k]] = k;
  for (int k = 0; k < nv; ++k)
    factory.insertVertex(positions[insertionIndex[k]]);

  auto vertexIndex = [&](std::array<int,dim> c) {
    int v = 0;
    for (int d = dim-1; d >= 0; --d)
      v = v*(n+1) + c[d];
    return (unsigned int)(vertexAt[v]);
  };

  // orientation of a simplex given by its vertex positions
  auto simplexOrientation = [&](const std::vector<unsigned int>& vertices) {
    FieldMatrix<double,dim,dim> J;
    for (int j = 0; j < dim; ++j)
      for (int d = 0; d < dim; ++d)
        J[j][d] = positions[insertionIndex[vertices[j+1]]][d] - positions[insertionIndex[vertices[0]]][d];
    return J.determinant();
  };

  // all orientation-preserving symmetries of the reference cube as maps of corner indices
  std::vector<std::vector<int>> cubeSymmetries;
  {
    std::array<int,dim> perm;
    std::iota(perm.begin(), perm.end(), 0);
    do {
      int inversions = 0;
      for (int i = 0; i < dim; ++i)
        for (int j = i+1; j < dim; ++j)
          inversions += perm[i] > perm[j];
      for (int flips = 0; flips < (1 << dim); ++flips) {
        if ((inversions + std::popcount((unsigned int)flips)) % 2 != 0)
          continue;
        std::vector<int> map(1 << dim);
        for (int k = 0; k < (1 << dim); ++k) {
          int image = 0;
          for (int d = 0; d < dim; ++d)
            image |= ((((k >> perm[d]) & 1) ^ ((flips >> d) & 1)) << d);
          map[k] = image;
        }
        cubeSymmetries.push_back(map);
      }
    } while (std::next_permutation(perm.begin(), perm.end()));
  }

  std::array<int,dim> cell{};
  for (int cellIndex = 0; cellIndex < Dune::power(n, dim); ++cellIndex) {
    for (int d = 0, r = cellIndex; d < dim; ++d, r /= n)
      cell[d] = r % n;

    auto corner = [&](int k) {
      std::array<int,dim> c = cell;
      for (int d = 0; d < dim; ++d)
        c[d] += (k >> d) & 1;
      return vertexIndex(c);
    };

    if (isSimplexCell(cell)) {
      // Kuhn triangulation of the cell
      std::array<int,dim> perm;
      std::iota(perm.begin(), perm.end(), 0);
      do {
        std::vector<unsigned int> vertices{corner(0)};
        int k = 0;
        for (int j = 0; j < dim; ++j) {
          k |= (1 << perm[j]);
          vertices.push_back(corner(k));
        }
        if (simplexOrientation(vertices) < 0)
          std::swap(vertices[0], vertices[1]);
        // random even permutation of the vertices
        std::shuffle(vertices.begin(), vertices.end(), rng);
        if (simplexOrientation(vertices) < 0)
          std::swap(vertices[0], vertices[1]);
        factory.insertElement(GeometryTypes::simplex(dim), vertices);
      } while (std::next_permutation(perm.begin(), perm.end()));
    }
    else {
      const auto& symmetry = cubeSymmetries[rng() % cubeSymmetries.size()];
      std::vector<unsigned int> vertices(1 << dim);
      for (int k = 0; k < (1 << dim); ++k)
        vertices[k] = corner(symmetry[k]);
      factory.insertElement(GeometryTypes::cube(dim), vertices);
    }
  }

  return factory.createGrid();
}


// Check that the interpolation of a polynomial of order p is exact
template<class Basis, class F>
void checkPolynomialInterpolation (TestSuite& testSuite, const Basis& basis, F f, std::string name)
{
  TestSuite t(name);
  std::vector<double> coefficients;
  interpolate(basis, coefficients, f);

  auto fh = makeDiscreteGlobalBasisFunction<double>(basis, coefficients);
  auto localFh = localFunction(fh);
  double maxError = 0;
  for (const auto& e : elements(basis.gridView())) {
    localFh.bind(e);
    const auto geometry = e.geometry();
    for (const auto& qp : QuadratureRules<double,Basis::GridView::dimension>::rule(e.type(), 4))
      maxError = std::max(maxError, std::abs(localFh(qp.position()) - f(geometry.global(qp.position()))));
  }
  t.check(maxError < 1e-8) << "interpolation error " << maxError;
  testSuite.subTest(t);
}


template<class GridView>
void testUniformOrder (TestSuite& testSuite, const GridView& gridView, unsigned int maxOrder, std::string name)
{
  constexpr int dim = GridView::dimension;
  for (unsigned int p = 1; p <= maxOrder; ++p) {
    auto basis = makeBasis(gridView, BasisFactory::lobatto(p));
    testSuite.subTest(checkBasis(basis, EnableContinuityCheck()));

    // polynomial of total degree p
    auto f = [p](const FieldVector<double,dim>& x) {
      double s = 0.3;
      for (int d = 0; d < dim; ++d)
        s += (d+1) * x[d];
      return std::pow(s, p);
    };
    checkPolynomialInterpolation(testSuite, basis, f, name + " interpolation p=" + std::to_string(p));
  }
}


template<class GridView>
void testAdaptiveOrder (TestSuite& testSuite, const GridView& gridView, unsigned int maxOrder, std::string name)
{
  constexpr int dim = GridView::dimension;
  using Element = typename GridView::template Codim<0>::Entity;
  MultipleCodimMultipleGeomTypeMapper<GridView> elementMapper(gridView, mcmgElementLayout());

  // pseudo-random element orders that can be changed after the basis is created
  std::vector<unsigned int> elementOrders(elementMapper.size());
  for (std::size_t i = 0; i < elementOrders.size(); ++i)
    elementOrders[i] = 1 + (i*7 + 3) % maxOrder;
  auto elementOrder = [&](const Element& e) { return elementOrders[elementMapper.index(e)]; };

  // constant order: same dimension as the uniform basis
  {
    TestSuite t(name + " adaptive with constant order");
    for (unsigned int p = 1; p <= maxOrder; ++p) {
      auto uniform = makeBasis(gridView, BasisFactory::lobatto(p));
      auto adaptive = makeBasis(gridView, BasisFactory::lobatto([p](const auto&) { return p; }));
      t.check(uniform.dimension() == adaptive.dimension())
        << "dimension of adaptive basis " << adaptive.dimension() << " != " << uniform.dimension();
      testSuite.subTest(checkBasis(adaptive, EnableContinuityCheck()));
    }
    testSuite.subTest(t);
  }

  // variable order with minimum rule
  {
    auto basis = makeBasis(gridView, BasisFactory::lobatto(elementOrder));
    testSuite.subTest(checkBasis(basis, EnableContinuityCheck()));

    // the minimum rule is satisfied
    TestSuite t(name + " minimum rule");
    for (const auto& e : elements(gridView)) {
      const auto& re = referenceElement(e);
      for (int c = 1; c < dim; ++c)
        for (int i = 0; i < re.size(c); ++i)
          t.check(basis.preBasis().subEntityOrder(e, i, c) <= elementOrder(e), "minimum rule violated");
    }
    testSuite.subTest(t);

    // linear functions are reproduced
    auto f = [](const FieldVector<double,dim>& x) { return 0.5 + x[0] - 2*x[dim-1]; };
    checkPolynomialInterpolation(testSuite, basis, f, name + " adaptive interpolation");

    // p-refinement: increase the order on some elements and update the basis
    const auto oldDimension = basis.dimension();
    for (std::size_t i = 0; i < elementOrders.size(); i += 3)
      elementOrders[i] = std::min(elementOrders[i] + 1, maxOrder + 1);
    basis.update(gridView);
    TestSuite tr(name + " p-refinement");
    tr.check(basis.dimension() >= oldDimension, "dimension decreased after p-refinement");
    testSuite.subTest(tr);
    testSuite.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  // override of the sub-entity orders
  {
    auto subEntityOrder = [&](const Element& e, int i, int c) {
      return std::max(1u, elementOrder(e) - (unsigned int)(i % 2));
    };
    auto basis = makeBasis(gridView, BasisFactory::lobatto(elementOrder, subEntityOrder));
    testSuite.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  // anisotropic orders (cube elements only)
  bool cubeGrid = true;
  for (auto type : gridView.indexSet().types(0))
    cubeGrid = cubeGrid && type.isCube();
  if (cubeGrid) {
    auto anisotropicOrder = [&](const Element& e) {
      std::array<unsigned int,dim> p;
      for (int k = 0; k < dim; ++k)
        p[k] = 1 + (elementMapper.index(e) + k) % maxOrder;
      return p;
    };
    auto basis = makeBasis(gridView, BasisFactory::lobatto(anisotropicOrder));
    testSuite.subTest(checkBasis(basis, EnableContinuityCheck()));
  }
}


int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);
  TestSuite testSuite;

  {
    auto grid = createShuffledGrid<UGGrid<2>>(3, [](auto) { return true; });
    testUniformOrder(testSuite, grid->leafGridView(), 5, "triangles");
    testAdaptiveOrder(testSuite, grid->leafGridView(), 4, "triangles");
  }
  {
    auto grid = createShuffledGrid<UGGrid<2>>(3, [](auto) { return false; });
    testUniformOrder(testSuite, grid->leafGridView(), 5, "quadrilaterals");
    testAdaptiveOrder(testSuite, grid->leafGridView(), 4, "quadrilaterals");
  }
  {
    auto grid = createShuffledGrid<UGGrid<2>>(4, [](auto c) { return c[0] < 2; });
    testUniformOrder(testSuite, grid->leafGridView(), 4, "mixed 2d");
    testAdaptiveOrder(testSuite, grid->leafGridView(), 4, "mixed 2d");
  }
  {
    auto grid = createShuffledGrid<UGGrid<3>>(2, [](auto) { return true; });
    testUniformOrder(testSuite, grid->leafGridView(), 4, "tetrahedra");
    testAdaptiveOrder(testSuite, grid->leafGridView(), 4, "tetrahedra");
  }
  {
    auto grid = createShuffledGrid<UGGrid<3>>(2, [](auto) { return false; });
    testUniformOrder(testSuite, grid->leafGridView(), 4, "hexahedra");
    testAdaptiveOrder(testSuite, grid->leafGridView(), 3, "hexahedra");
  }

  return testSuite.exit();
}
