// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <array>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dune/grid/yaspgrid.hh>

#include "stokes.hh"

#define BLOCKEDBASIS 2

// { using_namespace_dune_begin }
using namespace Dune;
using namespace Dune::Functions::Examples::Stokes;
// { using_namespace_dune_end }

// { main_begin }
int main (int argc, char *argv[]) try
{
  // Set up MPI, if available
  MPIHelper::instance(argc, argv);
  // { mpi_setup_end }

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  // { grid_setup_begin }
  const int dim = 2;
  using GridType = YaspGrid<dim>;
  FieldVector<double,dim> upperRight = {1, 1};
  std::array<int,dim> elements = {{4, 4}};
  GridType grid(upperRight,elements);

  using GridView = typename GridType::LeafGridView;
  GridView gridView = grid.leafGridView();
  // { grid_setup_end }

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////

  // { function_space_basis_begin }
  using namespace Functions::BasisFactory;

  constexpr std::size_t p = 1; // pressure order for Taylor-Hood

#if BLOCKEDBASIS
  auto taylorHoodBasis = makeBasis(
          gridView,
          composite(
            power<dim>(
              lagrange<p+1>(),
              blockedInterleaved()),
            lagrange<p>()
          ));
  run("taylor-hood",taylorHoodBasis,ContainerLayout<BLOCKED>{});
#else
  auto taylorHoodBasis = makeBasis(
          gridView,
          composite(
            power<dim>(
              lagrange<p+1>(),
              flatInterleaved()),
            lagrange<p>()
          ));
  run("taylor-hood",taylorHoodBasis,ContainerLayout<FLAT>{});
#endif
  // { function_space_basis_end }
 }
// Error handling
 catch (Exception& e) {
    std::cout << e.what() << std::endl;
 }
