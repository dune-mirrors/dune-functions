#pragma once

#include <array>
#include <tuple>
#include <type_traits>
#include <vector>

#include <dune/common/indices.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/functions/common/utility.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

namespace Dune {

struct StdTraits
{
  template <class T>
  using DynamicVector = std::vector<T>;

  template <class B, std::size_t N>
  using PowerVector = std::array<B,N>;

  template <class... Rows>
  using CompositeVector = std::tuple<Rows...>;

  template <class T, std::size_t N>
  using LeafBlockVector = std::vector<FieldVector<T,int(N)>>;
};

namespace Functions {

template <int dim>
class BasisFactories
{
  using GridType = YaspGrid<dim>;

public:
  static const std::size_t K = 1;
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
          lagrange<K+1>(),
          flatLexicographic()),
        lagrange<K>(),
        blockedLexicographic()
      ));
  }

  template <class T>
  auto vector(index_constant<0>) const
  {
    return StdTraits::CompositeVector<
      StdTraits::DynamicVector<T>,
      StdTraits::DynamicVector<T>
      >{};
  }

  auto basis(index_constant<1>) const // Root: blockedLexicographic, Velocity: flatInterleaved
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<dim>(
          lagrange<K+1>(),
          flatInterleaved()),
        lagrange<K>(),
        blockedLexicographic()
      ));
  }

  template <class T>
  auto vector(index_constant<1>) const
  {
    return StdTraits::CompositeVector<
      StdTraits::DynamicVector<T>,
      StdTraits::DynamicVector<T>
      >{};
  }

  auto basis(index_constant<2>) const // Root: blockedLexicographic, Velocity: blockedLexicographic
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<dim>(
          lagrange<K+1>(),
          blockedLexicographic()),
        lagrange<K>(),
        blockedLexicographic()
      ));
  }

  template <class T>
  auto vector(index_constant<2>) const
  {
    using Vector = StdTraits::CompositeVector<
      StdTraits::PowerVector<StdTraits::DynamicVector<T>, dim>,
      StdTraits::DynamicVector<T>
      >;
    return Vector{};
  }

  auto basis(index_constant<3>) const // Root: blockedLexicographic, Velocity: blockedInterleaved
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<dim>(
          lagrange<K+1>(),
          blockedInterleaved()),
        lagrange<K>(),
        blockedLexicographic()
      ));
  }

  template <class T>
  auto vector(index_constant<3>) const
  {
    using Vector = StdTraits::CompositeVector<
      StdTraits::LeafBlockVector<T, dim>,
      StdTraits::DynamicVector<T>
      >;
    return Vector{};
  }

  auto basis(index_constant<4>) const // Root: flatLexicographic, Velocity/Pressure: flatLexicographic
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<dim>(
          lagrange<K+1>(),
          flatLexicographic()),
        power<1>(
          lagrange<K>(),
          flatLexicographic()),
        flatLexicographic()
      ));
  }

  template <class T>
  auto vector(index_constant<4>) const
  {
    using Vector = StdTraits::DynamicVector<T>;
    return Vector{};
  }

  auto basis(index_constant<5>) const // Root: blockedLexicographic, Velocity/Pressure: blockedLexicographic
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<dim>(
          lagrange<K+1>(),
          blockedLexicographic()),
        power<1>(
          lagrange<K>(),
          blockedLexicographic()),
        blockedLexicographic()
      ));
  }

  template <class T>
  auto vector(index_constant<5>) const
  {
    using Vector = StdTraits::CompositeVector<
      StdTraits::PowerVector<StdTraits::DynamicVector<T>, dim>,
      StdTraits::PowerVector<StdTraits::DynamicVector<T>, 1>
      >;
    return Vector{};
  }

  auto basis(index_constant<6>) const // Root: blockedLexicographic, Velocity/Pressure: on the same level
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        lagrange<K+1>(),
        lagrange<K+1>(),
        lagrange<K>(),
        blockedLexicographic()
      ));
  }

  template <class T>
  auto vector(index_constant<6>) const
  {
    using Vector = StdTraits::CompositeVector<
      StdTraits::DynamicVector<T>,
      StdTraits::DynamicVector<T>,
      StdTraits::DynamicVector<T>
      >;
    return Vector{};
  }

  auto basis(index_constant<7>) const
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<1>(
          power<dim>(
            lagrange<K+1>(),
            blockedLexicographic() ),
          blockedLexicographic() ),
        lagrange<K>(),
        blockedLexicographic()
      ));
  }

  template <class T>
  auto vector(index_constant<7>) const
  {
    using Vector = StdTraits::CompositeVector<
      StdTraits::PowerVector<StdTraits::PowerVector<StdTraits::DynamicVector<T>, dim>, 1>,
      StdTraits::DynamicVector<T>
      >;
    return Vector{};
  }

  auto basis(index_constant<8>) const
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<1>(
          power<dim>(
            lagrange<K+1>(),
            blockedLexicographic() ),
          blockedLexicographic() ),
        power<1>(
          lagrange<K>(),
          blockedLexicographic() ),
        blockedLexicographic()
      ));
  }

  template <class T>
  auto vector(index_constant<8>) const
  {
    using Vector = StdTraits::CompositeVector<
      StdTraits::PowerVector<StdTraits::PowerVector<StdTraits::DynamicVector<T>, dim>, 1>,
      StdTraits::PowerVector<StdTraits::DynamicVector<T>, 1>
      >;
    return Vector{};
  }

  auto basis(index_constant<9>) const
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<2>(lagrange<K+1>()),  // Cahn-Hilliard equation
        composite(                  // Stokes equation
          power<dim>(
            lagrange<K+1>()
          ),
          lagrange<K>()
        )
    ));
  }

  template <class T>
  auto vector(index_constant<9>) const
  {
    using Vector = StdTraits::CompositeVector<
      StdTraits::LeafBlockVector<T, 2>,
      StdTraits::CompositeVector<
        StdTraits::LeafBlockVector<T, dim>,
        StdTraits::DynamicVector<T>
        >
      >;
    return Vector{};
  }

  auto basis(index_constant<10>) const
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        lagrange<K+1>(),            // Diffusion equation
        power<2>(                   // Cahn-Hilliard equation
          lagrange<K+1>(),
          blockedLexicographic()),
        composite(                  // Stokes equation
          power<dim>(
            lagrange<K+1>(),
            blockedLexicographic()
          ),
          lagrange<K>(),
          blockedLexicographic()
        ),
        blockedLexicographic()
    ));
  }

  template <class T>
  auto vector(index_constant<10>) const
  {
    using Vector = StdTraits::CompositeVector<
      StdTraits::DynamicVector<T>,
      StdTraits::PowerVector<StdTraits::DynamicVector<T>, 2>,
      StdTraits::CompositeVector<
        StdTraits::PowerVector<StdTraits::DynamicVector<T>, dim>,
        StdTraits::DynamicVector<T>
        >
      >;
    return Vector{};
  }

  auto false_basis(index_constant<0>) const
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(gridView(),
      composite(
        power<dim>(
          lagrange<K+1>(),
          blockedLexicographic()),
        power<1>(
          lagrange<K>(),
          blockedLexicographic()),
        flatLexicographic()
      ));
  }

private:
  GridType grid_;
};

}} // end namespace Dune::Functions
