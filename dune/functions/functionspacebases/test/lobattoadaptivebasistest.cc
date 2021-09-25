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

//! Enforce that the polynomial degree on entities shared by two elements is the minimum of the
//! polynomial degrees of the element interiors.
template <class GridView, class Orders>
void enforceMinimumRule (GridView const& gridView, Orders& orders)
{
  const int dim = GridView::dimension;
  if constexpr(dim > 1) {
    if (orders[1].empty()) // no facet orders assigned
      return;

    const auto& indexSet = gridView.indexSet();
    for (auto const& e : elements(gridView)) {
      const std::uint8_t inside = orders[0][indexSet.index(e)];

      auto refElem = referenceElement(e);
      for (auto const& is : intersections(gridView, e))
      {
        if (is.neighbor()) {
          const std::uint8_t outside = orders[0][indexSet.index(is.outside())];
          const std::uint8_t min_p = std::min(inside, outside);

          const int i = is.indexInOutside();
          const int c = 1;
          const auto idx = indexSet.subIndex(e,i,c);
          for (int k = 0; k < dim-1; ++k)
            orders[c][idx] = std::min(min_p, orders[c][idx]);

          if constexpr(dim > 2) {
            if (orders[2].empty()) // no sub-facet orders assigned
              continue;

            for (int cc = 2; cc < dim; ++cc) {
              for (int ii = 0; ii < refElem.size(i,c,cc); ++ii) {
                const auto sub_idx = indexSet.subIndex(e, refElem.subEntity(i,c,ii,cc), cc);
                orders[cc][sub_idx] = std::min(orders[c][idx], orders[cc][sub_idx]);
              }
            }
          }
        }
      }
    }
  }
}

template<int dim>
void test (Dune::TestSuite& testSuite)
{
  std::cout << "dim=" << dim << std::endl;
  using Grid = Dune::YaspGrid<dim>;

  Dune::FieldVector<double,dim> lower; lower = 0.0;
  Dune::FieldVector<double,dim> upper; upper = 1.0;
  auto elems = Dune::filledArray<dim,unsigned int>(2);
  elems[0] = 2;

  auto gridPtr = StructuredGridFactory<Grid>::createCubeGrid(lower,upper,elems);
  auto gridView = gridPtr->leafGridView();

  using namespace Dune::Functions::BasisFactory;

  LobattoEntityOrders orders{gridView};

  auto basis = makeBasis(gridView, lobatto(orders));
  for (unsigned int p = 1; p < 6; ++p) {
    orders.clear();
    for (auto const& e : elements(gridView))
      orders.set(e, p); // set interior order to p

    basis.preBasis().initializeIndices();
    basis.preBasis().debug();
    std::cout << "  p=" << p << " => basis.dimension=" << basis.dimension() << std::endl;

    testSuite.subTest(checkBasis(basis, EnableContinuityCheck()));
  }
}


void test2d (Dune::TestSuite& testSuite)
{
  std::cout << "flipped element in 2d" << std::endl;

  const int dim = 2;
  using Grid = Dune::UGGrid<dim>;

  Dune::GridFactory<Grid> factory;
  factory.insertVertex({0.0,0.0});
  factory.insertVertex({1.0,0.0});
  factory.insertVertex({2.0,0.0});
  factory.insertVertex({0.0,1.0});
  factory.insertVertex({1.0,1.0});
  factory.insertVertex({2.0,1.0});

  factory.insertElement(GeometryTypes::cube(2), {0u,1u,3u,4u});
  factory.insertElement(GeometryTypes::cube(2), {4u,5u,1u,2u}); // flipped in y-direction

  auto gridPtr = factory.createGrid();
  auto gridView = gridPtr->leafGridView();

  using namespace Dune::Functions::BasisFactory;

  LobattoEntityOrders orders{gridView};
  auto basis = makeBasis(gridView, lobatto(orders));

  for (unsigned int pb = 1; pb < 6; ++pb) {
    for (auto const& e : elements(gridView))
      orders.set(e, pb);

    for (unsigned int pe = 1; pe <= pb; ++pe) {
      for (auto const& e : elements(gridView)) {
        auto refElem = referenceElement(e);
        for (int i = 0; i < refElem.size(1); ++i)
          orders.set(e.subEntity<1>(i), pe);
      }

      basis.preBasis().initializeIndices();
      std::cout << "  pb=" << pb << ", pe=" << pe
                << " => basis.dimension=" << basis.dimension() << std::endl;
      testSuite.subTest(checkBasis(basis, EnableContinuityCheck()));
    }
  }
}


int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite testSuite;
  test<1>(testSuite);
  test<2>(testSuite);
  test<3>(testSuite);

  test2d(testSuite);

  return testSuite.exit();
}
