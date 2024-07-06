// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/basixbasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <basix/e-lagrange.h>

using namespace Dune;
using namespace Dune::Functions;

int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test;

  using namespace Dune::Functions::BasisFactory;

  { // dim = 2
    const int dim = 2;
    using Grid = Dune::UGGrid<dim>;
    using Factory = Dune::StructuredGridFactory<Grid>;

    { // simplex grid
      auto grid = Factory::createSimplexGrid({0.0,0.0}, {1.0,1.0}, {2,2});
      grid->globalRefine(2);
      auto gridView = grid->leafGridView();

      for (int degree = 1; degree < 4; ++degree)
      {
        std::cout << "triangle (deg=" << degree << "):" << std::endl;
        auto lagrange_tri = ::basix::element::create_lagrange<double>(::basix::cell::type::triangle, degree,
          ::basix::element::lagrange_variant::equispaced, false);
        auto basis_tri = makeBasis(gridView, BasisFactory::basix<RangeClass::scalar>(std::move(lagrange_tri)));
        test.subTest(checkBasis(basis_tri, EnableContinuityCheck()));
      }
    }

    { // cube grid
      auto grid = Factory::createCubeGrid({0.0,0.0}, {1.0,1.0}, {2,2});
      grid->globalRefine(2);
      auto gridView = grid->leafGridView();

      for (int degree = 1; degree < 4; ++degree)
      {
        std::cout << "quadrilateral (deg=" << degree << "):" << std::endl;
        auto lagrange_quad = ::basix::element::create_lagrange<double>(::basix::cell::type::quadrilateral, degree,
          ::basix::element::lagrange_variant::equispaced, false);
        auto basis_quad= makeBasis(gridView, BasisFactory::basix<RangeClass::scalar>(std::move(lagrange_quad)));
        test.subTest(checkBasis(basis_quad, EnableContinuityCheck()));
      }
    }
  }

  return test.exit();
}
