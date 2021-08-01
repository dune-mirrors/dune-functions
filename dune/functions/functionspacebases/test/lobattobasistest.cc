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

#include <dune/functions/functionspacebases/lobattobasis.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

using namespace Dune;
using namespace Dune::Functions;

template <class Basis>
void print (Basis const& basis)
{
  std::cout << "basis.dimension = " << basis.dimension() << std::endl;

  auto localView = basis.localView();
  for (auto const& e : elements(basis.gridView()))
  {
    std::cout << "element " << basis.gridView().indexSet().index(e) << " {" << std::endl;
    localView.bind(e);

    auto const& node = localView.tree();
    auto const& localFE = node.finiteElement();

    for (std::size_t i = 0; i < localFE.size(); ++i)
      std::cout << "  " << i << ": " << localFE.localCoefficients().localKey(i) << " => "
                << localView.index(node.localIndex(i)) << std::endl;

    std::cout << "}" << std::endl;
  }
}

template <class GridView>
void printGridView (GridView const& gridView)
{
  const int dim = GridView::dimension;

  std::cout << "grid = {" << std::endl;
  for (int c = 0; c <= dim; ++c)
    for (auto&& t : gridView.indexSet().types(c))
      std::cout << "  size[" << t << "] = " << gridView.size(t) << std::endl;
  auto const& indexSet = gridView.indexSet();
  for (auto const& e : elements(gridView))
  {
    std::cout << "  " << indexSet.index(e) << ": " << e.geometry().center() << std::endl;

    auto refElem = referenceElement(e);
    for (int c = 0; c <= dim; ++c) {
      std::cout << "    codim " << c << ": [";
      for (int i = 0; i < refElem.size(c); ++i) {
        std::cout << " " << indexSet.subIndex(e,i,c);
      }
      std::cout << " ]" << std::endl;
    }
  }
  std::cout << "}" << std::endl;
}


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
  // printGridView(gridView);

  using namespace Dune::Functions::BasisFactory;

  for (unsigned int p = 1; p < 6; ++p) {
    auto basis = makeBasis(gridView, lobatto(p));
    testSuite.subTest(checkBasis(basis, EnableContinuityCheck()));
  }

  for (unsigned int p = 1; p < 6; ++p) {
    auto basis = makeBasis(gridView, lobatto(LobattoOrders<dim>{std::uint8_t(p)}));
    testSuite.subTest(checkBasis(basis, EnableContinuityCheck()));
  }
}


int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite testSuite;
  test<1>(testSuite);
  test<2>(testSuite);
  test<3>(testSuite);

  return testSuite.exit();
}
