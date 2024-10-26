// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/functions/functionspacebases/multidomain.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

int main(int argc, char** argv)
{
  using namespace Dune::Functions::BasisFactory;
  using namespace Dune::Functions::MultiDomain;

  Dune::MPIHelper::instance(argc, argv);

  /* ---- setup a mesh, partition it and define subdomains ---- */

  constexpr int dim = 2;
  Dune::YaspGrid<dim> grid{ {1.0,1.0}, {3,1} };
  const auto & gridView = grid.leafGridView();

  std::vector<int> partitions(grid.size(0), 0);
  auto const& indexSet = grid.leafIndexSet();
  for (auto const& e : elements(grid.leafGridView())) {
    auto center = e.geometry().center();
    partitions[indexSet.index(e)] =
      center[0] < 0.45 ? 0 :
      (center[0] > 0.55 ? 1 : 2);
  }
  // partitions 0,2,1 ... 0: domain 0, 1: domain 1, 2: overlap
  // -> subdomains domain0 = {0,2}, domain1 = {1,2}, interface = {2}
  auto domainInfo = createPartitionedDomainInfo(std::move(partitions), {{0,2}, {1,2}, {2}});

  /*      _______
    dom0: |_|_|X|
          _______
    dom1: |X|_|_|
   */
  // CutCellDomainInfo domainInfo(levelsets, {{0,2}, {1,2}, {2}});


  /* ---- run tests on a single subdomain or in a multidomain setup ---- */

  Dune::TestSuite test;

  {
    auto stokesPreBasis = composite(power<dim>(lagrange<2>()), lagrange<1>());
    auto darcyPreBasis = lagrange<1>();
    auto multiDomainPreBasis =
      multiDomainComposite(domainInfo,
        restrict(stokesPreBasis, subdomain(0)),
        restrict(darcyPreBasis, subdomain(1)));

    auto basis = makeBasis(gridView,multiDomainPreBasis);

    using namespace Dune::Indices;
    const auto& restrictedGridView0 = basis.preBasis().subPreBasis(_0).gridView();
    const auto& restrictedGridView1 = basis.preBasis().subPreBasis(_1).gridView();

    // 3x1 mesh
    // subdomains : 2x1 mesh -> 2 cells, 7 edges, 6 vertices
    test.require(restrictedGridView0.size(0) == 2, "check # of subdomain 0 cells");
    test.require(restrictedGridView1.size(0) == 2, "check # of subdomain 1 cells");
    test.require(restrictedGridView0.size(1) == 7, "check # of subdomain 0 edges");
    test.require(restrictedGridView1.size(1) == 7, "check # of subdomain 1 edges");
    test.require(restrictedGridView0.size(2) == 6, "check # of subdomain 0 vertices");
    test.require(restrictedGridView1.size(2) == 6, "check # of subdomain 1 vertices");
  }

  {
    auto stokesPreBasis = composite(power<dim>(lagrange<2>()), lagrange<1>());
    auto darcyPreBasis = lagrange<1>();
    auto multiDomainPreBasis =
      multiDomainComposite(domainInfo,
        restrict(stokesPreBasis, subdomain(0)),
        restrict(darcyPreBasis, subdomain(1)));

    auto basis = makeBasis(gridView,multiDomainPreBasis);

    auto checkGlobalSize = [](const auto & basis)
    {
      // TH(dom0)   :
      //         v  :  30 DOFs ( 2 per cell, edge, vertex)
      //         p  :   6 DOFs
      // Lag(dom1)  :   6 DOFs
      // total      :  42 DOFs
      int sz = basis.dimension();
      return sz == 42;
    };

    auto checkLocalSize = [](const auto & basis)
    {
      // Stokes(dom0): 2 x (4 x vertex + 4 x edge + cell) + 4 x vertex -> 22 DOFs
      // Darcy(dom1):  4 x vertex -> 4 DOFs
      std::array<int,3> refSizes = { 22, 26, 4 };
      std::array<int,3> sizes;
      auto & gv = basis.gridView(); // -> wrapped original gridView
      auto localView = basis.localView();
      int i = 0;
      for (auto && e : elements(gv))
      {
        localView.bind(e); // internally: domainInfo.bind(e), then child.bind(e)
        sizes[i++] = localView.size();
      };
      return std::ranges::equal(refSizes, sizes);
    };

    test.subTest(checkBasis(basis));
    test.require(checkGlobalSize(basis), "check global size of restricted Stokes-Darcy basis");
    test.require(checkLocalSize(basis), "check local sizes of restricted Stokes-Darcy basis");
  }

  return test.exit();
}
