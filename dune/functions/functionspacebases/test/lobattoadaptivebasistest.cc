// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/lobattoadaptivebasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

using namespace Dune;
using namespace Dune::Functions;

template<int dim>
void test (Dune::TestSuite& testSuite)
{
  using Grid = Dune::YaspGrid<dim>;

  Dune::FieldVector<double,dim> lower; lower = 0.0;
  Dune::FieldVector<double,dim> upper; upper = 1.0;
  auto elems = Dune::filledArray<dim,unsigned int>(2);
  elems[0] = 2;

  auto gridPtr = StructuredGridFactory<Grid>::createCubeGrid(lower, upper,elems);
  auto gridView = gridPtr->leafGridView();

  using namespace Dune::Functions::BasisFactory;

  auto basis = makeBasis(gridView, lobattoAdaptive(LobattoOrders<dim>{GeometryTypes::cube(dim)}));
  for (unsigned int p = 1; p < 6; ++p) {
    LobattoOrders<dim> o{GeometryTypes::cube(dim), p};
    for (auto const& e : elements(gridView))
      basis.preBasis().setOrders(e, o);
    basis.preBasis().enforceMinimumRule();
    basis.preBasis().initializeIndices();

    testSuite.subTest(checkBasis(basis, EnableContinuityCheck()));
  }
}


int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite testSuite;
  test<1>(testSuite);
  test<2>(testSuite);
  // test<3>(testSuite);

  return testSuite.exit();
}