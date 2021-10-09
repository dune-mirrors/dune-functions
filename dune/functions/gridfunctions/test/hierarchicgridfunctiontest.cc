// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/istl/bvector.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/hierarchicgridfunction.hh>

#include <dune/functions/gridfunctions/test/gridfunctiontest.hh>

using namespace Dune;
using namespace Dune::Functions;
using namespace Dune::Functions::Test;

int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  TestSuite suite;

  // Generate grid for testing
  Dune::YaspGrid<2> grid({1.0,1.0}, {10,10});
  grid.globalRefine(2);
  auto gridView1 = grid.levelGridView(0);
  auto gridView2 = grid.leafGridView();

  using namespace Functions::BasisBuilder;
  auto basis1 = makeBasis(gridView1, lagrange<2>());
  auto basis2 = makeBasis(gridView2, lagrange<2>());

  using Range = FieldVector<double,1>;
  Dune::BlockVector<Range> coeff1, coeff2;

  {
    // Inner test function f is a polynomial of degree 2.
    auto f = [](const auto& x){
      Range y;
      for (typename Range::size_type i = 0; i < y.size(); ++i)
        y[i] = (x[i]+i)*x[i];
      return y;
    };

    // Interpolate f wrt basis.
    interpolate(basis1, coeff1, f);
    auto gf1 = Dune::Functions::makeDiscreteGlobalBasisFunction<Range>(basis1, coeff1);

    // Compute integral (order 4 is sufficient).
    double integral1 = integrateGridViewFunction(gridView1, gf1, 4);

    HierarchicGridFunction hgf{gridView2, gf1};
    interpolate(basis2, coeff2, hgf);
    auto gf2 = Dune::Functions::makeDiscreteGlobalBasisFunction<Range>(basis2, coeff2);

    // Compute integral (order 4 is sufficient).
    double integral2 = integrateGridViewFunction(gridView2, gf2, 4);

    std::cout << "integral1 = " << integral1 << std::endl;
    std::cout << "integral2 = " << integral2 << std::endl;
    std::cout << "diff = " << std::abs(integral1 - integral2) << std::endl;
    suite.check(std::abs(integral1 - integral2) < 4*std::numeric_limits<double>::epsilon());

    // suite.check(
    //     checkGridViewFunction(gridView2, hgf, integral, 4),
    //     "Check if HierarchicGridFunction has correct integral");
  }


  return suite.exit();

} catch ( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
  return 1;
}
