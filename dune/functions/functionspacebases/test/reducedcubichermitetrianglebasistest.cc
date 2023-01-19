// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/reducedcubichermitetrianglebasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>

using namespace Dune;

int main(int argc, char *argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test;

  const int dim = 2;

  using Grid = UGGrid<dim>;
  using GridView = Grid::LeafGridView;
  FieldVector<double, dim> l(1);
  std::shared_ptr<Grid> grid = StructuredGridFactory<Grid>::createSimplexGrid({0.0, 0.0}, l, {{10, 10}});
  auto gridView = grid->leafGridView();

  // check ReducedCubicHermiteTriangleBasis created 'manually'
  {
    Functions::ReducedCubicHermiteTriangleBasis<GridView> basis(gridView);
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    test.subTest(checkBasis(basis, EnableVertexJacobianContinuityCheck()));
  }

  // check ReducedCubicHermiteTriangleBasis created using basis builder mechanism
  {
    using namespace Functions::BasisFactory;
    auto basis = makeBasis(gridView, reducedCubicHermiteTriangle());
    test.subTest(checkBasis(basis, EnableContinuityCheck()));
    test.subTest(checkBasis(basis, EnableVertexJacobianContinuityCheck()));
  }

  return test.exit();
}