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
};

namespace Functions {
namespace BasisFactory {

/**
 * \brief Collection of basis factories with additional info to be used in tests.
 *
 * This collection of factories to generate global bases with various index merging
 * strategies. These factories represent different ways to encode a taylor-hood
 * basis, i.e., `dim` velocity components and a pressure component.
 *
 * In addition to the basis factories, the class provides methods for constructing
 * vector containers compatible with the bases. Therefore, we use the `StdTraits`
 * type definition for different kinds of vectors.
 **/
template <int dim>
class BasisFactories
{
  using Traits = StdTraits;

public:
  static const std::size_t num_bases = 16;
  static const std::size_t num_false_bases = 1;

  // Root: blockedLexicographic, Velocity: flatLexicographic
  static auto basis(index_constant<0>)
  {
    return composite(
        power<dim>(
          lagrange<2>(),
          flatLexicographic()),
        lagrange<1>(),
        blockedLexicographic()
      );
  }

  template <class T>
  static auto vector(index_constant<0>)
  {
    return Traits::CompositeVector<
      Traits::DynamicVector<T>,
      Traits::DynamicVector<T>
      >{};
  }


  // Root: blockedLexicographic, Velocity: flatInterleaved
  static auto basis(index_constant<1>)
  {
    return composite(
        power<dim>(
          lagrange<2>(),
          flatInterleaved()),
        lagrange<1>(),
        blockedLexicographic()
      );
  }

  template <class T>
  static auto vector(index_constant<1>)
  {
    return Traits::CompositeVector<
      Traits::DynamicVector<T>,
      Traits::DynamicVector<T>
      >{};
  }


  // Root: blockedLexicographic, Velocity: blockedLexicographic
  static auto basis(index_constant<2>)
  {
    return composite(
        power<dim>(
          lagrange<2>(),
          blockedLexicographic()),
        lagrange<1>(),
        blockedLexicographic()
      );
  }

  template <class T>
  static auto vector(index_constant<2>)
  {
    using Vector = Traits::CompositeVector<
      Traits::PowerVector<Traits::DynamicVector<T>, dim>,
      Traits::DynamicVector<T>
      >;
    return Vector{};
  }


  // Root: blockedLexicographic, Velocity: blockedInterleaved
  static auto basis(index_constant<3>)
  {
    return composite(
        power<dim>(
          lagrange<2>(),
          blockedInterleaved()),
        lagrange<1>(),
        blockedLexicographic()
      );
  }

  template <class T>
  static auto vector(index_constant<3>)
  {
    using Vector = Traits::CompositeVector<
      Traits::DynamicVector<Dune::FieldVector<T, dim>>,
      Traits::DynamicVector<T>
      >;
    return Vector{};
  }


  // Root: flatLexicographic, Velocity/Pressure: flatLexicographic
  static auto basis(index_constant<4>)
  {
    return composite(
        power<dim>(
          lagrange<2>(),
          flatLexicographic()),
        power<1>(
          lagrange<1>(),
          flatLexicographic()),
        flatLexicographic()
      );
  }

  template <class T>
  static auto vector(index_constant<4>)
  {
    using Vector = Traits::DynamicVector<T>;
    return Vector{};
  }


  // Root: blockedLexicographic, Velocity/Pressure: blockedLexicographic
  static auto basis(index_constant<5>)
  {
    return composite(
        power<dim>(
          lagrange<2>(),
          blockedLexicographic()),
        power<1>(
          lagrange<1>(),
          blockedLexicographic()),
        blockedLexicographic()
      );
  }

  template <class T>
  static auto vector(index_constant<5>)
  {
    using Vector = Traits::CompositeVector<
      Traits::PowerVector<Traits::DynamicVector<T>, dim>,
      Traits::PowerVector<Traits::DynamicVector<T>, 1>
      >;
    return Vector{};
  }


  // Root: blockedLexicographic, Velocity/Pressure: on the same level
  static auto basis(index_constant<6>)
  {
    return composite(
        lagrange<2>(),
        lagrange<2>(),
        lagrange<1>(),
        blockedLexicographic()
      );
  }

  template <class T>
  static auto vector(index_constant<6>)
  {
    using Vector = Traits::CompositeVector<
      Traits::DynamicVector<T>,
      Traits::DynamicVector<T>,
      Traits::DynamicVector<T>
      >;
    return Vector{};
  }


  static auto basis(index_constant<7>)
  {
    return composite(
        power<1>(
          power<dim>(
            lagrange<2>(),
            blockedLexicographic() ),
          blockedLexicographic() ),
        lagrange<1>(),
        blockedLexicographic()
      );
  }

  template <class T>
  static auto vector(index_constant<7>)
  {
    using Vector = Traits::CompositeVector<
      Traits::PowerVector<Traits::PowerVector<Traits::DynamicVector<T>, dim>, 1>,
      Traits::DynamicVector<T>
      >;
    return Vector{};
  }


  static auto basis(index_constant<8>)
  {
    return composite(
        power<1>(
          power<dim>(
            lagrange<2>(),
            blockedLexicographic() ),
          blockedLexicographic() ),
        power<1>(
          lagrange<1>(),
          blockedLexicographic() ),
        blockedLexicographic()
      );
  }

  template <class T>
  static auto vector(index_constant<8>)
  {
    using Vector = Traits::CompositeVector<
      Traits::PowerVector<Traits::PowerVector<Traits::DynamicVector<T>, dim>, 1>,
      Traits::PowerVector<Traits::DynamicVector<T>, 1>
      >;
    return Vector{};
  }


  static auto basis(index_constant<9>)
  {
    return composite(
        power<2>(lagrange<2>()),  // Cahn-Hilliard equation
        composite(                // Stokes equation
          power<dim>(
            lagrange<2>()
          ),
          lagrange<1>()
        )
    );
  }

  template <class T>
  static auto vector(index_constant<9>)
  {
    using Vector = Traits::CompositeVector<
      Traits::DynamicVector<FieldVector<T, 2>>,
      Traits::CompositeVector<
        Traits::DynamicVector<FieldVector<T, dim>>,
        Traits::DynamicVector<T>
        >
      >;
    return Vector{};
  }


  static auto basis(index_constant<10>)
  {
    return composite(
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
    );
  }

  template <class T>
  static auto vector(index_constant<10>)
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


  static auto basis(index_constant<11>)
  {
    return power<2>(
        power<2>(lagrange<2>(), blockedLexicographic()),
        flatInterleaved()
        );
  }

  template <class T>
  static auto vector(index_constant<11>)
  {
    using Vector = Traits::PowerVector<Traits::DynamicVector<T>, 4>;
    return Vector{};
  }


  static auto basis(index_constant<12>)
  {
    return power<2>(
        power<2>(lagrange<2>(), blockedLexicographic()),
        flatLexicographic()
        );
  }

  template <class T>
  static auto vector(index_constant<12>)
  {
    using Vector = Traits::PowerVector<Traits::DynamicVector<T>, 4>;
    return Vector{};
  }


  static auto basis(index_constant<13>)
  {
    return power<2>(
        composite(
          power<1>(power<1>(lagrange<1>())),
          power<2>(lagrange<1>()),
          power<3>(lagrange<1>())
        ),
        flatLexicographic()
      );
  }

  template <class T>
  static auto vector(index_constant<13>)
  {
    using Vector = Traits::CompositeVector<
      Traits::DynamicVector<Traits::PowerVector<Dune::FieldVector<T,1>,1>>,
      Traits::DynamicVector<Dune::FieldVector<T,2>>,
      Traits::DynamicVector<Dune::FieldVector<T,3>>,
      Traits::DynamicVector<Traits::PowerVector<Dune::FieldVector<T,1>,1>>,
      Traits::DynamicVector<Dune::FieldVector<T,2>>,
      Traits::DynamicVector<Dune::FieldVector<T,3>>
      >;
    return Vector{};
  }


  static auto basis(index_constant<14>)
  {
    return power<2>(
        composite(
          power<1>(power<1>(lagrange<1>())),
          power<2>(lagrange<1>()),
          power<3>(lagrange<1>())
        ),
        flatInterleaved()
      );
  }

  template <class T>
  static auto vector(index_constant<14>)
  {
    using Vector = Traits::CompositeVector<
      Traits::DynamicVector<Traits::PowerVector<Dune::FieldVector<T,1>,1>>,
      Traits::DynamicVector<Traits::PowerVector<Dune::FieldVector<T,1>,1>>,
      Traits::DynamicVector<Dune::FieldVector<T,2>>,
      Traits::DynamicVector<Dune::FieldVector<T,2>>,
      Traits::DynamicVector<Dune::FieldVector<T,3>>,
      Traits::DynamicVector<Dune::FieldVector<T,3>>
      >;
    return Vector{};
  }


  static auto basis(index_constant<15>)
  {
    return composite(
        composite(
          power<2>(lagrange<1>()),
          power<1>(power<2>(lagrange<1>()))
        ),
        composite(
          power<1>(power<1>(lagrange<1>())),
          power<2>(lagrange<1>()),
          power<3>(lagrange<1>())
        ),
        flatLexicographic()
      );
  }

  template <class T>
  static auto vector(index_constant<15>)
  {
    using Vector = Traits::CompositeVector<
      Traits::DynamicVector<Dune::FieldVector<T,2>>,
      Traits::DynamicVector<Traits::PowerVector<Dune::FieldVector<T,1>,2>>,
      Traits::DynamicVector<Traits::PowerVector<Dune::FieldVector<T,1>,1>>,
      Traits::DynamicVector<Dune::FieldVector<T,2>>,
      Traits::DynamicVector<Dune::FieldVector<T,3>>
      >;
    return Vector{};
  }


  static auto false_basis(index_constant<0>)
  {
    return composite(
        power<dim>(
          lagrange<2>(),
          blockedLexicographic()),
        power<1>(
          lagrange<1>(),
          blockedLexicographic()),
        flatLexicographic()
      );
  }
};

}}} // end namespace Dune::Functions::BasisFactory
