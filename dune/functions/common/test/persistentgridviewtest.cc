// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <array>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/common/persistentgridview.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

using namespace Dune;
using namespace Dune::Functions;

int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  auto test = Dune::TestSuite();

  // Generate grid for testing
  using Grid = Dune::YaspGrid<2>;
  auto grid = Grid({1., 1.},{10, 10});
  auto gridView = grid.leafGridView();

  // Create a basis for testing
  using Dune::Functions::BasisFactory::lagrange;
  using Dune::Functions::BasisFactory::makeBasis;
  auto persistentBasis = makeBasis(Dune::Functions::Experimental::persistentGridView(gridView), lagrange<3>());

  auto accessFromBasis = [&](){
    auto localView = persistentBasis.localView();
    for(const auto& e : Dune::elements(gridView))
      localView.bind(e);
  };

  // Check if filling the PersistentGridView lazily
  // in mutable state works
  test.requireNoThrow(accessFromBasis)
    << "Filling PersistentGridView in mutable state failed.";

  persistentBasis.gridView().indexSet().setMutable(false);

  // Check if accessing the PersistentGridView in non-mutable state works
  test.requireNoThrow(accessFromBasis)
    << "Accessing PersistentGridView indices in non-mutable state failed.";

  grid.globalRefine(1);

  // Check if accessing the PersistentGridView in non-mutable state works
  // after grid refinement
  {
    auto accessFromCoarseGridView = [&](){
      // Create coarse grid view because you cannot iterate over
      // a PersistentGridView
      auto coarseGridView = grid.levelGridView(grid.maxLevel()-1);
      auto localView = persistentBasis.localView();
      for(const auto& e : Dune::elements(coarseGridView))
        localView.bind(e);
    };

    test.requireNoThrow(accessFromCoarseGridView)
      << "Accessing PersistentGridView indices in non-mutable state failed after grid refinement.";
  }

  // Check if the coarse elements are known to the PersistentGridView
  {
    bool containsCoarseElements = true;
    auto pgv = persistentBasis.gridView();

    auto coarseGridView = grid.levelGridView(grid.maxLevel()-1);
    for (const auto& e : Dune::elements(coarseGridView)) {
      containsCoarseElements = containsCoarseElements && pgv.contains(e);
    }
    test.require(containsCoarseElements, "Check contains()")
      << "Some elements that are part of the PersistentGridView are not contained.";
  }

  // Reversely, check if the new elements are unknown to the PersistentGridView
  {
    bool containsNoFineElements = true;
    auto pgv = persistentBasis.gridView();

    for (const auto& e : Dune::elements(grid.leafGridView())) {
      containsNoFineElements = containsNoFineElements && not pgv.contains(e);
    }
    test.require(containsNoFineElements, "Check contains()")
      << "Some elements that are part of the PersistentGridView are not contained.";
  }

  return test.exit();
}
