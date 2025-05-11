// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dune/grid/multidomaingrid.hh>

#include <dune/istl/bvector.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/transformedgridviewprebasis.hh>
#include <dune/functions/functionspacebases/transformedgridviewbasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>

int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;

  // Create a 2D grid
  std::shared_ptr grid2d = [](){
    const int dim = 2;
    using Grid = Dune::UGGrid<dim>;

    Dune::GridFactory<Grid> factory;
    factory.insertVertex({0,0});
    factory.insertVertex({0,1});
    factory.insertVertex({1,0});
    factory.insertVertex({1,1});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {0,1,2});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {1,2,3});

    auto grid = factory.createGrid();
    grid->globalRefine(4);
    return grid;
  }();

  // Create a multi-domain grid with 2 sub-domains
  auto mdgrid2d = []<class Grid>(const std::shared_ptr<Grid>& grid) {
    // lets have only 2 sub-domains
    using MDTraits = Dune::mdgrid::FewSubDomainsTraits<Grid::dimension, 3>;
    using MDGrid = Dune::mdgrid::MultiDomainGrid<Grid, MDTraits>;
    auto mdgrid = std::make_unique<MDGrid>(grid);

    mdgrid->startSubDomainMarking();
    for (const auto &cell :
        elements(mdgrid->leafGridView(), Dune::Partitions::interior)) {
      if (cell.geometry().center()[0] < 0.5)
        mdgrid->addToSubDomain(0, cell);
      else
        mdgrid->addToSubDomain(1, cell);
    }

    mdgrid->preUpdateSubDomains();
    mdgrid->updateSubDomains();
    mdgrid->postUpdateSubDomains();
    return mdgrid;
  }(std::move(grid2d));

  using namespace Dune::Functions::BasisFactory;

  // Defnie extension of pre-basis factory
  // -> pose pre-basis on a sub-domain grid and transform it to a pre-basis posed on a multi-domain grid
  static constexpr auto restrictPreBasisFactory = [](auto&& pre_basis_factory, std::size_t subDomain)
  {
    return [=]<class MultiDomainGridView>(const MultiDomainGridView& multiDomainGridView){
      using SubDomainGridView = typename MultiDomainGridView::Grid::SubDomainGrid::LeafGridView;
      return Dune::Functions::TransformedGridViewPreBasis{
        pre_basis_factory(multiDomainGridView.grid().subDomain(subDomain).leafGridView()),
        multiDomainGridView,
        [subDomain](const MultiDomainGridView& multiDomainGridView) { return multiDomainGridView.grid().subDomain(subDomain).leafGridView(); },
        [](const auto& subDomainGridView, const MultiDomainGridView& multiDomainGridView, const auto& entity)
        {
          if (subDomainGridView.contains(entity))
            return subDomainGridView.grid().subDomainEntity(entity);
          else
            return typename SubDomainGridView::template Codim<0>::Entity{};
        }
      };
    };
  };

  // Define restriction of basis factory
  // -> transform basis posed on a multi-domain grid to basis posed on a sub-domain grid
  static constexpr auto restrictBasis = []<class Basis>(const Basis& basis, std::size_t subDomain)
  {
    auto gridView = basis.gridView().grid().subDomain(subDomain).leafGridView();
    return Dune::Functions::TransformedGridViewBasis{
      basis,
      basis.gridView().grid().subDomain(subDomain).leafGridView(),
      [](const auto& subDomainGridView){
        return subDomainGridView.grid().multiDomainGrid().leafGridView();
      },
      []<class MultiDomainGridView>(const MultiDomainGridView& multiDomainGridView, const auto& subDomainGridView, const auto& entity)
        -> const MultiDomainGridView::template Codim<0>::Entity&
      {
        return multiDomainGridView.grid().multiDomainEntity(entity);
      }
    };
  };

  auto basis = makeBasis(mdgrid2d->leafGridView(),
    composite(
      restrictPreBasisFactory(lagrange<1>(), 0),
      restrictPreBasisFactory(lagrange<2>(), 1),
      lagrange<1>(),
      flatLexicographic()
    )
  );

  Dune::BlockVector<double> x;
  auto xBackend = Dune::Functions::istlVectorBackend(x);
  xBackend.resize(basis);

  auto writeBasisDiscreteFunction = [&](std::string name, const auto& basis, auto f){
    interpolate(basis, x, f);

    auto xFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(basis, x);

    Dune::SubsamplingVTKWriter vtkWriter(basis.gridView(), Dune::refinementLevels(2));

    vtkWriter.addVertexData(
      xFunction,
      Dune::VTK::FieldInfo("x", Dune::VTK::FieldInfo::Type::scalar, 1));

    vtkWriter.write(name);
  };

  using namespace Dune::Indices;
  writeBasisDiscreteFunction("sub-domain-0", restrictBasis(subspaceBasis(basis, Dune::TypeTree::treePath(_0)), 0),   [=](auto x){ return std::sin(10*x.two_norm()); });
  writeBasisDiscreteFunction("sub-domain-1", restrictBasis(subspaceBasis(basis, Dune::TypeTree::treePath(_1)), 1),   [=](auto x){ return std::sin(10*x[0]) + std::cos(10*x[1]); });
  writeBasisDiscreteFunction("multi-domain-0",               subspaceBasis(basis, Dune::TypeTree::treePath(_0)    ), [=](auto x){ return std::sin(10*x.two_norm()); });
  writeBasisDiscreteFunction("multi-domain-1",               subspaceBasis(basis, Dune::TypeTree::treePath(_1)    ), [=](auto x){ return std::sin(10*x[0]) + std::cos(10*x[1]); });
  writeBasisDiscreteFunction("multi-domain-2",               subspaceBasis(basis, Dune::TypeTree::treePath(_2)    ), [=](auto x){ return 2 + x.two_norm(); });
}
catch (Dune::Exception& e)
{
  std::cout << e.what() << std::endl;
}
