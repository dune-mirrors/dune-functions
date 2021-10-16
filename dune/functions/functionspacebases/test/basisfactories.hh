#pragma once

#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/functions/common/utility.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

namespace Dune {

struct ISTLTraits
{
  template <class T>
  using DynamicVector = BlockVector<T>;

  template <class B, std::size_t N>
  using PowerVector = BlockVector<B>;

  template <bool same, class... Rows>
  struct CompositeVectorImpl;

  template <class... Rows>
  struct CompositeVectorImpl<false, Rows...>
  {
    using type = MultiTypeBlockVector<Rows...>;
  };

  template <class Row0, class... Rows>
  struct CompositeVectorImpl<true, Row0, Rows...>
  {
    using type = PowerVector<Row0, (sizeof...(Rows)+1)>;
  };

  // if all components are the same, using a PowerVector otherwise a MultiTypeVector
  template <class... Rows>
  using CompositeVector
    = typename CompositeVectorImpl<Dune::Functions::isAllSame<Rows...>::value, Rows...>::type;

  template <class T, std::size_t N>
  using LeafBlockVector = BlockVector<FieldVector<T,int(N)>>;
};


namespace Functions {

template <int dim>
class BasisFactories
{
  using GridType = YaspGrid<dim>;

public:
  static const std::size_t K = 1;
  static const std::size_t num_bases = 10;
  static const std::size_t num_false_bases = 1;

  BasisFactories()
    : grid_({1.0, 1.0}, {2, 2})
  {}

  GridType& grid() { return grid_; }
  auto gridView() const { return grid_.leafGridView(); }

  auto basis(index_constant<0>) const // Root: blockedLexicographic, Velocity: flatLexicographic
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(
      gridView(),
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
    return ISTLTraits::CompositeVector<
      ISTLTraits::DynamicVector<T>,
      ISTLTraits::DynamicVector<T>
      >{};
  }

  auto basis(index_constant<1>) const // Root: blockedLexicographic, Velocity: flatInterleaved
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(
      gridView(),
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
    return ISTLTraits::CompositeVector<
      ISTLTraits::DynamicVector<T>,
      ISTLTraits::DynamicVector<T>
      >{};
  }

  auto basis(index_constant<2>) const // Root: blockedLexicographic, Velocity: blockedLexicographic
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(
      gridView(),
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
    using Vector = ISTLTraits::CompositeVector<
      ISTLTraits::PowerVector<ISTLTraits::DynamicVector<T>, dim>,
      ISTLTraits::DynamicVector<T>
      >;
    return Vector{};
  }

  auto basis(index_constant<3>) const // Root: blockedLexicographic, Velocity: blockedInterleaved
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(
      gridView(),
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
    using Vector = ISTLTraits::CompositeVector<
      ISTLTraits::LeafBlockVector<T, dim>,
      ISTLTraits::DynamicVector<T>
      >;
    return Vector{};
  }

  auto basis(index_constant<4>) const // Root: flatLexicographic, Velocity/Pressure: flatLexicographic
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(
      gridView(),
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
    using Vector = ISTLTraits::DynamicVector<T>;
    return Vector{};
  }

  auto basis(index_constant<5>) const // Root: blockedLexicographic, Velocity/Pressure: blockedLexicographic
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(
      gridView(),
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
    using Vector = ISTLTraits::CompositeVector<
      ISTLTraits::PowerVector<ISTLTraits::DynamicVector<T>, dim>,
      ISTLTraits::PowerVector<ISTLTraits::DynamicVector<T>, 1>
      >;
    return Vector{};
  }

  auto basis(index_constant<6>) const // Root: blockedLexicographic, Velocity/Pressure: on the same level
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(
      gridView(),
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
    using Vector = ISTLTraits::CompositeVector<
      ISTLTraits::DynamicVector<T>,
      ISTLTraits::DynamicVector<T>,
      ISTLTraits::DynamicVector<T>
      >;
    return Vector{};
  }

  auto basis(index_constant<7>) const
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(
      gridView(),
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
    using Vector = ISTLTraits::CompositeVector<
      ISTLTraits::PowerVector<ISTLTraits::PowerVector<ISTLTraits::DynamicVector<T>, dim>, 1>,
      ISTLTraits::DynamicVector<T>
      >;
    return Vector{};
  }

  auto basis(index_constant<8>) const
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(
      gridView(),
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
    using Vector = ISTLTraits::CompositeVector<
      ISTLTraits::PowerVector<ISTLTraits::PowerVector<ISTLTraits::DynamicVector<T>, dim>, 1>,
      ISTLTraits::PowerVector<ISTLTraits::DynamicVector<T>, 1>
      >;
    return Vector{};
  }

  auto basis(index_constant<9>) const // Root: blockedLexicographic, Velocity: blockedInterleaved
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(
      gridView(),
      composite(
        power<2>(lagrange<K+1>()),  // cahn-hilliard basis
        composite(                  // stokes basis
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
    using Vector = ISTLTraits::CompositeVector<
      ISTLTraits::LeafBlockVector<T, 2>,
      ISTLTraits::CompositeVector<
        ISTLTraits::LeafBlockVector<T, dim>,
        ISTLTraits::DynamicVector<T>
        >
      >;
    return Vector{};
  }

  auto false_basis(index_constant<0>) const // Root: flatLexicographic, Velocity/Pressure: blockedLexicographic
  {
    using namespace Dune::Functions::BasisFactory;
    return makeBasis(
      gridView(),
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
