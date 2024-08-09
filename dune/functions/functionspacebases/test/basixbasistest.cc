// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/functions/functionspacebases/basixbasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

using namespace Dune;
using namespace Dune::Functions;

template <class Grid>
struct PerturbedGrid
{
  using ctype = typename Grid::ctype ;
  static constexpr int dimension = Grid::dimension;
  static constexpr int dimensionworld = Grid::dimensionworld;
};

template <class Grid>
class Dune::GridFactory<PerturbedGrid<Grid>>
    : public GridFactory<Grid>
{
  using Super = GridFactory<Grid>;
  static constexpr int dim = Grid::dimension;
  static constexpr int dow = Grid::dimensionworld;
  using ctype = typename Grid::ctype;

public:
  template <class... Args>
  explicit GridFactory (Args&&... args)
    : Super(std::forward<Args>(args)...)
    , rd_()
    , g_(rd_())
  {}

  void insertElement(const Dune::GeometryType& type,
                     const std::vector<unsigned int>& vertices) override
  {
    std::vector<unsigned int> new_vertices(vertices);
    std::shuffle(new_vertices.begin(), new_vertices.end(), g_);
    Super::insertElement(type,new_vertices);
  }
  using Super::insertElement;

  void insertBoundarySegment(const std::vector<unsigned int>& vertices)
  {
    std::vector<unsigned int> new_vertices(vertices);
    std::shuffle(new_vertices.begin(), new_vertices.end(), g_);
    Super::insertBoundarySegment(new_vertices);
  }
  using Super::insertBoundarySegment;

private:
  std::random_device rd_;
  std::mt19937 g_;
};

int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test;

  using namespace Dune::Functions::BasisFactory;

  { // dim = 2
    const int dim = 2;
    using Grid = Dune::UGGrid<dim>;
    using GridFactory = Dune::GridFactory<Grid>;
    using Factory = Dune::StructuredGridFactory<Grid>;
    // using GridFactory = Dune::GridFactory<PerturbedGrid<Grid>>;
    // using Factory = Dune::StructuredGridFactory<PerturbedGrid<Grid>>;

    { // simplex grid
      auto gridFactory = GridFactory{};
      Factory::createSimplexGrid(gridFactory, {0.0,0.0}, {1.0,1.0}, {2,2});
      auto grid = gridFactory.createGrid();
      //grid->globalRefine(1);
      auto gridView = grid->leafGridView();

      for (int degree = 1; degree < 4; ++degree)
      {
        std::cout << "triangle Lagrange (deg=" << degree << "):" << std::endl;
        auto basis_lag = makeBasis(gridView, basix_lagrange(degree));
        test.subTest(checkBasis(basis_lag, EnableContinuityCheck()));

        std::cout << "triangle Lagrange-dg (deg=" << degree << "):" << std::endl;
        auto basis_lagdg = makeBasis(gridView, basix_lagrangedg(degree));
        test.subTest(checkBasis(basis_lagdg));

        std::cout << "triangle Lobatto (deg=" << degree << "):" << std::endl;
        auto basis_lob = makeBasis(gridView, basix_lagrange(degree,::basix::element::lagrange_variant::gll_warped));
        test.subTest(checkBasis(basis_lob, EnableContinuityCheck()));

        std::cout << "triangle Nedelec (deg=" << degree << "):" << std::endl;
        auto basis_ned = makeBasis(gridView, basix_nedelec(degree));
        test.subTest(checkBasis(basis_ned, EnableTangentialContinuityCheck()));

        std::cout << "triangle Raviart-Thomas (deg=" << degree << "):" << std::endl;
        auto basis_rt = makeBasis(gridView, basix_rt(degree));
        test.subTest(checkBasis(basis_rt, EnableNormalContinuityCheck()));

        std::cout << "triangle Brezzi-Douglas-Marini (deg=" << degree << "):" << std::endl;
        auto basis_bdm = makeBasis(gridView, basix_bdm(degree));
        test.subTest(checkBasis(basis_bdm, EnableNormalContinuityCheck()));

        if (degree == 1) {
          std::cout << "triangle Crouzeix-Raviart (deg=" << degree << "):" << std::endl;
          auto basis_cr = makeBasis(gridView, basix_cr());
          test.subTest(checkBasis(basis_cr, EnableCenterContinuityCheck()));
        }

#if 0 // not yet implemented
        if (degree < 3) {
          std::cout << "triangle Regge (deg=" << degree << "):" << std::endl;
          auto basis_regge = makeBasis(gridView, basix_regge(degree));
          test.subTest(checkBasis(basis_regge));
          // TODO: Found a constant zero basis function for degree >= 3
        }

        if (degree < 3) {
          std::cout << "triangle Hellan-Herrmann-Johnson (deg=" << degree << "):" << std::endl;
          auto basis_hhj = makeBasis(gridView, basix_hhj(degree));
          test.subTest(checkBasis(basis_hhj));
          // TODO: Found a constant zero basis function for degree >= 3
        }
#endif
      }
    }

    { // cube grid
      auto gridFactory = GridFactory{};
      Factory::createCubeGrid(gridFactory, {0.0,0.0}, {1.0,1.0}, {2,2});
      auto grid = gridFactory.createGrid();
      //grid->globalRefine(1);
      auto gridView = grid->leafGridView();

      for (int degree = 1; degree < 4; ++degree)
      {
        std::cout << "quadrilateral (deg=" << degree << "):" << std::endl;
        auto basis_lag= makeBasis(gridView, basix_lagrange(degree));
        test.subTest(checkBasis(basis_lag, EnableContinuityCheck()));

        std::cout << "quadrilateral Lagrange-dg (deg=" << degree << "):" << std::endl;
        auto basis_lagdg = makeBasis(gridView, basix_lagrangedg(degree));
        test.subTest(checkBasis(basis_lagdg));

        std::cout << "quadrilateral Lobatto (deg=" << degree << "):" << std::endl;
        auto basis_lob = makeBasis(gridView, basix_lagrange(degree,::basix::element::lagrange_variant::gll_warped));
        test.subTest(checkBasis(basis_lob, EnableContinuityCheck()));
      }
    }
  }


  { // dim = 3
    const int dim = 3;
    using Grid = Dune::UGGrid<dim>;
    using GridFactory = Dune::GridFactory<Grid>;
    using Factory = Dune::StructuredGridFactory<Grid>;
    // using GridFactory = Dune::GridFactory<PerturbedGrid<Grid>>;
    // using Factory = Dune::StructuredGridFactory<PerturbedGrid<Grid>>;

    { // simplex grid
      auto gridFactory = GridFactory{};
      Factory::createSimplexGrid(gridFactory, {0.0,0.0,0.0}, {1.0,1.0,1.0}, {2,2,2});
      auto grid = gridFactory.createGrid();
      grid->globalRefine(1);
      auto gridView = grid->leafGridView();

      for (int degree = 1; degree < 5; ++degree)
      {
        std::cout << "tetrahedron (deg=" << degree << "):" << std::endl;
        auto basis_tet = makeBasis(gridView, basix_lagrange(degree));
        test.subTest(checkBasis(basis_tet, EnableContinuityCheck()));
      }
    }

    { // cube grid
      auto gridFactory = GridFactory{};
      Factory::createCubeGrid(gridFactory, {0.0,0.0,0.0}, {1.0,1.0,1.0}, {3,3,3});
      auto grid = gridFactory.createGrid();
      // grid->globalRefine(1);  // Crashes UGGrid
      auto gridView = grid->leafGridView();

      for (int degree = 1; degree < 5; ++degree)
      {
        std::cout << "hexahedron (deg=" << degree << "):" << std::endl;
        auto basis_hex = makeBasis(gridView, basix_lagrange(degree));
        test.subTest(checkBasis(basis_hex, EnableContinuityCheck()));
      }
    }
  }

  return test.exit();
}
