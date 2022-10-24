#pragma once

#include <array>
#include <tuple>
#include <type_traits>
#include <vector>

#include <dune/common/indices.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/functions/common/utility.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

namespace Dune {
namespace Functions {

template <int dim>
class BasisFactories
{
  using GridType = YaspGrid<dim>;

public:
  static const std::size_t num_bases = 11;
  static const std::size_t num_false_bases = 1;

  BasisFactories()
    : grid_({1.0, 1.0}, {2, 2})
  {}

  GridType& grid() { return grid_; }
  auto gridView() const { return grid_.leafGridView(); }

  auto basis(index_constant<0>) const // Root: blockedLexicographic, Velocity: flatLexicographic
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<dim>(
          lagrange<2>(),
          flatLexicographic()),
        lagrange<1>(),
        blockedLexicographic()
      ));
  }

  auto prefixes(index_constant<0>) const
  {
    using namespace Dune::Indices;
    using namespace Dune::TypeTree;
    return std::tuple{
      hybridTreePath(_0), hybridTreePath(_1), hybridTreePath(_0, 0u)
    };
  }


  auto basis(index_constant<1>) const // Root: blockedLexicographic, Velocity: flatInterleaved
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<dim>(
          lagrange<2>(),
          flatInterleaved()),
        lagrange<1>(),
        blockedLexicographic()
      ));
  }

  auto prefixes(index_constant<1>) const
  {
    using namespace Dune::Indices;
    using namespace Dune::TypeTree;
    return std::tuple{
      hybridTreePath(_0), hybridTreePath(_1), hybridTreePath(_0, 0u)
    };
  }


  auto basis(index_constant<2>) const // Root: blockedLexicographic, Velocity: blockedLexicographic
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<dim>(
          lagrange<2>(),
          blockedLexicographic()),
        lagrange<1>(),
        blockedLexicographic()
      ));
  }

  auto prefixes(index_constant<2>) const
  {
    using namespace Dune::Indices;
    using namespace Dune::TypeTree;
    return std::tuple{
      hybridTreePath(_0), hybridTreePath(_1), hybridTreePath(_0, 0u), hybridTreePath(_0, 1u),
      hybridTreePath(_0, 0u, 0u), hybridTreePath(_0, 1u, 0u)
    };
  }


  auto basis(index_constant<3>) const // Root: blockedLexicographic, Velocity: blockedInterleaved
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<dim>(
          lagrange<2>(),
          blockedInterleaved()),
        lagrange<1>(),
        blockedLexicographic()
      ));
  }

  auto prefixes(index_constant<3>) const
  {
    using namespace Dune::Indices;
    using namespace Dune::TypeTree;
    return std::tuple{
      hybridTreePath(_0), hybridTreePath(_1), hybridTreePath(_0, 0u), hybridTreePath(_0, 1u),
      hybridTreePath(_0, 0u, 0u), hybridTreePath(_0, 0u, 1u)
    };
  }


  auto basis(index_constant<4>) const // Root: flatLexicographic, Velocity/Pressure: flatLexicographic
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<dim>(
          lagrange<2>(),
          flatLexicographic()),
        power<1>(
          lagrange<1>(),
          flatLexicographic()),
        flatLexicographic()
      ));
  }

  auto prefixes(index_constant<4>) const
  {
    using namespace Dune::Indices;
    using namespace Dune::TypeTree;
    return std::tuple{
      hybridTreePath(0u)
    };
  }


  auto basis(index_constant<5>) const // Root: blockedLexicographic, Velocity/Pressure: blockedLexicographic
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<dim>(
          lagrange<2>(),
          blockedLexicographic()),
        power<1>(
          lagrange<1>(),
          blockedLexicographic()),
        blockedLexicographic()
      ));
  }

  auto prefixes(index_constant<5>) const
  {
    using namespace Dune::Indices;
    using namespace Dune::TypeTree;
    return std::tuple{
      hybridTreePath(_0), hybridTreePath(_1), hybridTreePath(_0, 0u), hybridTreePath(_0, 1u),
      hybridTreePath(_1, 0u), hybridTreePath(_0, 0u, 0u), hybridTreePath(_0, 0u, 1u),
      hybridTreePath(_1, 0u, 0u)
    };
  }


  auto basis(index_constant<6>) const // Root: blockedLexicographic, Velocity/Pressure: on the same level
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        lagrange<2>(),
        lagrange<2>(),
        lagrange<1>(),
        blockedLexicographic()
      ));
  }

  auto prefixes(index_constant<6>) const
  {
    using namespace Dune::Indices;
    using namespace Dune::TypeTree;
    return std::tuple{
      hybridTreePath(_0), hybridTreePath(_1), hybridTreePath(_2), hybridTreePath(_0, 0u),
      hybridTreePath(_1, 0u), hybridTreePath(_2, 0u)
    };
  }


  auto basis(index_constant<7>) const
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<1>(
          power<dim>(
            lagrange<2>(),
            blockedLexicographic() ),
          blockedLexicographic() ),
        lagrange<1>(),
        blockedLexicographic()
      ));
  }

  auto prefixes(index_constant<7>) const
  {
    using namespace Dune::Indices;
    using namespace Dune::TypeTree;
    return std::tuple{
      hybridTreePath(_0), hybridTreePath(_1), hybridTreePath(_0, 0u), hybridTreePath(_0, 0u, 0u),
      hybridTreePath(_0, 0u, 0u, 0u), hybridTreePath(_1, 0u)
    };
  }


  auto basis(index_constant<8>) const
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<1>(
          power<dim>(
            lagrange<2>(),
            blockedLexicographic() ),
          blockedLexicographic() ),
        power<1>(
          lagrange<1>(),
          blockedLexicographic() ),
        blockedLexicographic()
      ));
  }

  auto prefixes(index_constant<8>) const
  {
    using namespace Dune::Indices;
    using namespace Dune::TypeTree;
    return std::tuple{
      hybridTreePath(_0), hybridTreePath(_1), hybridTreePath(_0, 0u), hybridTreePath(_0, 0u, 0u),
      hybridTreePath(_0, 0u, 0u, 0u), hybridTreePath(_1, 0u), hybridTreePath(_1, 0u, 0u)
    };
  }


  auto basis(index_constant<9>) const
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<2>(lagrange<2>()),  // Cahn-Hilliard equation
        composite(                  // Stokes equation
          power<dim>(
            lagrange<2>()
          ),
          lagrange<1>()
        )
    ));
  }

  auto prefixes(index_constant<9>) const
  {
    using namespace Dune::Indices;
    using namespace Dune::TypeTree;
    return std::tuple{
      hybridTreePath(_0), hybridTreePath(_1), hybridTreePath(_0, 0u), hybridTreePath(_0, 0u, 0u),
      hybridTreePath(_1, _0), hybridTreePath(_1, _0, 0u), hybridTreePath(_1, _0, 0u, 0u),
      hybridTreePath(_1, _1), hybridTreePath(_1, _1, 0u)
    };
  }


  auto basis(index_constant<10>) const
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        lagrange<2>(),              // Diffusion equation
        power<2>(                   // Cahn-Hilliard equation
          lagrange<2>(),
          blockedLexicographic()),
        composite(                  // Stokes equation
          power<dim>(
            lagrange<2>(),
            blockedLexicographic()
          ),
          lagrange<1>(),
          blockedLexicographic()
        ),
        blockedLexicographic()
    ));
  }

  auto prefixes(index_constant<10>) const
  {
    using namespace Dune::Indices;
    using namespace Dune::TypeTree;
    return std::tuple{
      hybridTreePath(_0), hybridTreePath(_1), hybridTreePath(_2), hybridTreePath(_0, 0u),
      hybridTreePath(_1, 0u), hybridTreePath(_1, 0u, 0u),
      hybridTreePath(_2, _0), hybridTreePath(_2, _0, 0u), hybridTreePath(_2, _0, 0u, 0u),
      hybridTreePath(_2, _1), hybridTreePath(_2, _1, 0u)
    };
  }


  auto false_basis(index_constant<0>) const
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<dim>(
          lagrange<2>(),
          blockedLexicographic()),
        power<1>(
          lagrange<1>(),
          blockedLexicographic()),
        flatLexicographic()
      ));
  }

private:
  GridType grid_;
};

}} // end namespace Dune::Functions
