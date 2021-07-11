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

#include <dune/functions/functionspacebases/lobattobasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

using namespace Dune;
using namespace Dune::Functions;


template<int dim>
void test (Dune::TestSuite& testSuite)
{
  using Grid = Dune::UGGrid<dim>;

  Dune::FieldVector<double,dim> lower; lower = 0.0;
  Dune::FieldVector<double,dim> upper; upper = 0.0;
  auto elems = Dune::filledArray<unsigned int,dim>(10);

  auto gridPtr = StructuredGridFactory<Grid>::createCubeGrid(lower, upper,elems);
  auto gridView = gridPtr->leafGridView();

  auto basis1 = makeBasis(gridView, lobatto(1));
  auto basis2 = makeBasis(gridView, lobatto(LobattoOrder<dim>{2}));

  testSuite.subTest(checkBasis(basis1, EnableContinuityCheck()));
  testSuite.subTest(checkBasis(basis2, EnableContinuityCheck()));
}


int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite testSuite;
  test<1>(testSuite);
  test<2>(testSuite);
  test<3>(testSuite);

  return test.exit();
}
