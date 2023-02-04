// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/reducedcubichermitetrianglebasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>

using namespace Dune;

// Test function: The natural embedding of the 2d grid into R^3.
// Is there no way to do this shorter, e.g., with two lambdas?
class IdentityGridEmbedding
{
public:
  FieldVector<double,3> operator() (const FieldVector<double,2>& x) const
  {
    return {x[0], x[1], 0.0};
  }

  friend auto derivative(const IdentityGridEmbedding& p)
  {
    return [](const FieldVector<double,2>& x) { return FieldMatrix<double,3,2>({{1,0}, {0,1}, {0,0}}); };
  }
};

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
    test.subTest(checkBasis(basis, EnableContinuityCheck(), EnableVertexJacobianContinuityCheck()));
  }

  // check ReducedCubicHermiteTriangleBasis created using basis builder mechanism
  {
    using namespace Functions::BasisFactory;
    auto basis = makeBasis(gridView, reducedCubicHermiteTriangle());
    test.subTest(checkBasis(basis, EnableContinuityCheck(), EnableVertexJacobianContinuityCheck()));

    /**
     * @brief TODO: Test with PowerBasis leads to Segmentation fault!
     */
    auto powerBasis = makeBasis(gridView,
                                power<3>(
                                reducedCubicHermiteTriangle(),
                                blockedInterleaved()));

    // For debugging:
    auto localView = powerBasis.localView();
    for (const auto& e : elements(powerBasis.gridView()))
    {
      localView.bind(e);

      // This call works...
      std::cout << localView.tree().child(0).finiteElement().element().type() << std::endl;

      // ... this call crashes.
      std::cout << localView.tree().child(0).finiteElement().localBasis().element().type() << std::endl;
    }

    test.subTest(checkBasis(powerBasis, EnableContinuityCheck(), EnableVertexJacobianContinuityCheck()));

    // Test whether we can call 'interpolate'
    std::vector<FieldVector<double,3> > x;
    IdentityGridEmbedding identityGridEmbedding;
    Functions::interpolate(powerBasis, x, identityGridEmbedding);
  }

  return test.exit();
}
