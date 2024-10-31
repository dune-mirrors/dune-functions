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

int main(int argc, char** argv)
{
  using namespace Dune::Functions::BasisFactory;
  using namespace Dune::Functions::MultiDomain;

  Dune::MPIHelper::instance(argc, argv);

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
  // CutCellDomainInfo domainInfo(levelsets, {{0,2}, {1,2}, {2}});


  auto stokesPreBasis = restrict(composite(power<dim>(lagrange<2>()), lagrange<1>()), subdomain(0));
  // we definitely
  auto darcyPreBasis = restrict(lagrange<1>(), subdomain(1));
  auto multiDomainPreBasis = multiDomainComposite(domainInfo, stokesPreBasis, darcyPreBasis);
  auto fooDomainPreBasis = multiDomainComposite(domainInfo, darcyPreBasis);
  // auto basis = makeBasis(gridView,fooDomainPreBasis);
  auto basis = makeBasis(gridView,multiDomainPreBasis);

  basis.gridView(); // -> original gridView
  auto localView = basis.localView();
  int i = 0;
  for (auto && e : elements(gridView))
  {
    localView.bind(e); // -> domainInfo.bind(e), dann child.bind(e)
    std::cout << "cell " << i++
              << " has " << localView.size() << " DOFs" << std::endl;
  }

  using namespace Dune::Indices;
  std::cout << "DOM0: " << basis.preBasis().subPreBasis(_0).gridView().size(0) << " cells\n";
  std::cout << "DOM1: " << basis.preBasis().subPreBasis(_0).gridView().size(0) << " cells\n";
  std::cout << "DOM0: " << basis.preBasis().subPreBasis(_0).gridView().size(1) << " edges\n";
  std::cout << "DOM1: " << basis.preBasis().subPreBasis(_0).gridView().size(1) << " edges\n";
  std::cout << "DOM0: " << basis.preBasis().subPreBasis(_0).gridView().size(2) << " vertices\n";
  std::cout << "DOM1: " << basis.preBasis().subPreBasis(_0).gridView().size(2) << " vertices\n";

  // after setting a subdomain id for all elements, the number of DOFs has changed
  // 3x1 mesh
  // subdomains : 2x1 mesh -> 2 cells, 7 edges, 6 vertices
  // TH(dom0)   :
  //         v  :  30 DOFs ( 2 per cell, edge, vertex)
  //         p  :   6 DOFs
  // Lag(dom1)  :   6 DOFs
  // total      :  42 DOFs
  int sz = basis.dimension();
  std::cout << "Functionsspace has " << sz << " global DOFs" << std::endl;
}
