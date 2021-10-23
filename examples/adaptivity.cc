// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/istl/bvector.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/persistentgridview.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/coarsenedgridfunction.hh>
#include <dune/functions/gridfunctions/refinedgridfunction.hh>

#include <dune/functions/gridfunctions/test/gridfunctiontest.hh>

using namespace Dune;
using namespace Dune::Functions;
using namespace Dune::Functions::Test;

int main (int argc, char* argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  TestSuite suite;

  // Generate grid for testing
  Dune::YaspGrid<2> grid({1.0,1.0}, {1,1});
  auto gridView = grid.leafGridView();

  using namespace Functions::BasisBuilder;
  auto basis1 = makeBasis(gridView, lagrange<2>());
  auto basis2 = makeBasis(gridView, lagrange<2>());
  using Basis2 = decltype(basis2);

  using Range = FieldVector<double,1>;
  Dune::BlockVector<Range> coeff1, coeff2;

  // Inner test function f is a polynomial of degree 2.
  auto f = [](const auto& x){
    Range y;
    for (typename Range::size_type i = 0; i < y.size(); ++i)
      y[i] = (x[i]+i)*x[i];
    return y;
  };

  // Interpolate f wrt basis.
  interpolate(basis1, coeff1, f);

  auto oldGridView = Dune::Functions::Experimental::persistentGridView(gridView);
  auto oldBasis = makeBasis(oldGridView, lagrange<2>());
  auto gf1 = Dune::Functions::makeDiscreteGlobalBasisFunction<Range>(oldBasis, coeff1);

  grid.globalRefine(2);

  auto newGridView = grid.leafGridView();
  basis2.update(newGridView);

  coeff2.resize(sizeInfo(basis2));
  auto gf2 = Dune::Functions::makeDiscreteGlobalBasisFunction<Range>(basis2, coeff2);

  // Obtain a local view of f
  auto lf1 = localFunction(gf1);
  auto lf2 = localFunction(gf2);

  RefinedGridFunction r_gf1{newGridView, gf1};
  auto r_lf1 = localFunction(r_gf1);

  CoarsenedGridFunction c_gf1{newGridView, gf1};
  auto c_lf1 = localFunction(c_gf1);

  using Tree = typename Basis2::LocalView::Tree;
  using NTRE = Dune::Functions::HierarchicNodeToRangeMap;
  NTRE ntre{};

  using BV = Dune::Functions::Imp::AllTrueBitSetVector;
  BV bitVector{};
  // auto bitVector = istlVectorBackend(bv);

  auto localView = basis2.localView();
  for (const auto& e : elements(gridView))
  {
    localView.bind(e);

    if (e.isNew()) {
      c_lf1.bind(e);
      Dune::Functions::Imp::LocalInterpolateVisitor<Basis2, Tree, NTRE, decltype(coeff2), decltype(c_lf1), decltype(bitVector)> localInterpolateVisitor(basis2, coeff2, bitVector, c_lf1, localView, ntre);
      Dune::TypeTree::applyToTree(localView.tree(),localInterpolateVisitor);
    } else {
      r_lf1.bind(e);
      Dune::Functions::Imp::LocalInterpolateVisitor<Basis2, Tree, NTRE, decltype(coeff2), decltype(r_lf1), decltype(bitVector)> localInterpolateVisitor(basis2, coeff2, bitVector, r_lf1, localView, ntre);
      Dune::TypeTree::applyToTree(localView.tree(),localInterpolateVisitor);
    }
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
