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
  using Traits = StdTraits;

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

  template <class T>
  auto vector(index_constant<0>) const
  {
    return Traits::CompositeVector<
      Traits::DynamicVector<T>,
      Traits::DynamicVector<T>
      >{};
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

  template <class T>
  auto vector(index_constant<1>) const
  {
    return Traits::CompositeVector<
      Traits::DynamicVector<T>,
      Traits::DynamicVector<T>
      >{};
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

  template <class T>
  auto vector(index_constant<2>) const
  {
    using Vector = Traits::CompositeVector<
      Traits::PowerVector<Traits::DynamicVector<T>, dim>,
      Traits::DynamicVector<T>
      >;
    return Vector{};
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

  template <class T>
  auto vector(index_constant<3>) const
  {
    using Vector = Traits::CompositeVector<
      Traits::LeafBlockVector<T, dim>,
      Traits::DynamicVector<T>
      >;
    return Vector{};
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

  template <class T>
  auto vector(index_constant<4>) const
  {
    using Vector = Traits::DynamicVector<T>;
    return Vector{};
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

  template <class T>
  auto vector(index_constant<5>) const
  {
    using Vector = Traits::CompositeVector<
      Traits::PowerVector<Traits::DynamicVector<T>, dim>,
      Traits::PowerVector<Traits::DynamicVector<T>, 1>
      >;
    return Vector{};
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

  template <class T>
  auto vector(index_constant<6>) const
  {
    using Vector = Traits::CompositeVector<
      Traits::DynamicVector<T>,
      Traits::DynamicVector<T>,
      Traits::DynamicVector<T>
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
            lagrange<2>(),
            blockedLexicographic() ),
          blockedLexicographic() ),
        lagrange<1>(),
        blockedLexicographic()
      ));
  }

  template <class T>
  auto vector(index_constant<7>) const
  {
    using Vector = Traits::CompositeVector<
      Traits::PowerVector<Traits::PowerVector<Traits::DynamicVector<T>, dim>, 1>,
      Traits::DynamicVector<T>
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
            lagrange<2>(),
            blockedLexicographic() ),
          blockedLexicographic() ),
        power<1>(
          lagrange<1>(),
          blockedLexicographic() ),
        blockedLexicographic()
      ));
  }

  template <class T>
  auto vector(index_constant<8>) const
  {
    using Vector = Traits::CompositeVector<
      Traits::PowerVector<Traits::PowerVector<Traits::DynamicVector<T>, dim>, 1>,
      Traits::PowerVector<Traits::DynamicVector<T>, 1>
      >;
    return Vector{};
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

  template <class T>
  auto vector(index_constant<9>) const
  {
    using Vector = Traits::CompositeVector<
      Traits::LeafBlockVector<T, 2>,
      Traits::CompositeVector<
        Traits::LeafBlockVector<T, dim>,
        Traits::DynamicVector<T>
        >
      >;
    return Vector{};
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

  template <class T>
  auto vector(index_constant<10>) const
  {
    using Vector = Traits::CompositeVector<
      Traits::DynamicVector<T>,
      Traits::PowerVector<Traits::DynamicVector<T>, 2>,
      Traits::CompositeVector<
        Traits::PowerVector<Traits::DynamicVector<T>, dim>,
        Traits::DynamicVector<T>
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
