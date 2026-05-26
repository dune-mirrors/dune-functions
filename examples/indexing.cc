// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <iostream>

#include <dune/common/indices.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/utility/gridinfo.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/printgrid.hh>
// #include <dune/vtk/datacollectors/lagrangedatacollector.hh>
// #include <dune/vtk/vtkwriter.hh>

template <class GV>
void printGridView (const GV& gridView)
{
  auto const& indexSet = gridView.indexSet();
  for (auto const& e : elements(gridView))
  {
    auto refElem = referenceElement(e);
    std::cout << "Element " << indexSet.index(e) << ":" << std::endl;
    std::cout << "  Edges (" << e.subEntities(1) << "):" << std::endl;
    for (unsigned int i = 0; i < e.subEntities(1); ++i) {
      std::cout << "    Edge " << indexSet.subIndex(e,i,1) << ":" << std::endl;
      auto ei = e.template subEntity<1>(i);
      std::cout << "      Center: " << ei.geometry().center() << std::endl;
      std::cout << "      Vertices:";
      for (int ii = 0; ii < refElem.size(i,1,2); ++ii)
        std::cout << " " << indexSet.subIndex(e,refElem.subEntity(i,1,ii,2),2);
      std::cout << std::endl;
    }
    std::cout << "  Vertices (" << e.subEntities(2) << "):" << std::endl;
    for (unsigned int i = 0; i < e.subEntities(2); ++i) {
      std::cout << "    Vertex " << indexSet.subIndex(e,i,2) << ":" << std::endl;
      auto ei = e.template subEntity<2>(i);
      std::cout << "      Center: " << ei.geometry().center() << std::endl;
    }
  }
}

template <class Basis>
void printBasis (const Basis& basis, std::string name)
{
  using namespace Dune::Indices;
  std::cout << std::endl << name << std::endl << "----------------------" << std::endl;

  auto localView = basis.localView();
  auto const& indexSet = basis.gridView().indexSet();
  for (auto const& e : elements(basis.gridView()))
  {
    std::cout << "Element " << indexSet.index(e) << ":" << std::endl;

    localView.bind(e);
    auto const& tree = localView.tree();
    auto const& vNode = tree.child(_0);
    auto const& pNode = tree.child(_1);

    std::cout << "  Velocity:" << std::endl;
    for (std::size_t i = 0; i < vNode.degree(); ++i)
    {
      std::cout << "    V[" << i << "]:" << std::endl;
      auto const& viNode = vNode.child(i);
      for (std::size_t j = 0; j < viNode.size(); ++j) {
        auto localIndex = viNode.localIndex(j);
        auto globalIndex = localView.index(localIndex);
        std::cout << "      " << j << " => " << localIndex << " => " << globalIndex << std::endl;
      }
    }

    std::cout << "  Pressure:" << std::endl;
    for (std::size_t j = 0; j < pNode.size(); ++j) {
      auto localIndex = pNode.localIndex(j);
      auto globalIndex = localView.index(localIndex);
      std::cout << "      " << j << " => " << localIndex << " => " << globalIndex << std::endl;
    }
  }
}


int main (int argc, char *argv[])
{
  auto& helper = Dune::MPIHelper::instance(argc, argv);

  using GridType = Dune::UGGrid<2>;
  auto gridPtr = Dune::StructuredGridFactory<GridType>::createSimplexGrid({0.0,0.0},{1.0,1.0},{1u,1u});

  printGridView(gridPtr->leafGridView());

  using namespace Dune::Functions::BasisFactory;
  auto basis1 = makeBasis(gridPtr->leafGridView(),
    composite(
      power<2>(
        lagrange<2>(),
        blockedLexicographic()
      ),
      lagrange<1>(),
      blockedLexicographic()
      )
    );
  printBasis(basis1, "BLBL");


  auto basis2 = makeBasis(gridPtr->leafGridView(),
    composite(
      power<2>(
        lagrange<2>(),
        blockedInterleaved()
      ),
      lagrange<1>(),
      blockedLexicographic()
      )
    );
  printBasis(basis2, "BIBL");


  auto basis3 = makeBasis(gridPtr->leafGridView(),
    composite(
      power<2>(
        lagrange<2>(),
        flatInterleaved()
      ),
      lagrange<1>(),
      flatLexicographic()
      )
    );
  printBasis(basis3, "FIFL");

}