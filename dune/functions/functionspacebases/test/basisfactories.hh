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
 * This collection of factories can be used to generate global bases with various
 * index merging strategies. Some factories represent different ways to encode a
 * taylor-hood basis, i.e., a combination of velocity components and a pressure component.
 * Other bases are tailored to trigger complicated cases for the index-tree / size
 * construction.
 *
 * The basis factories are numbered from `[0, BasisFactory::size)` and can be accessed
 * using the static methods `basis(index)` where `index` is an `index_constant` in
 * that range.
 *
 * In addition to the basis factories, the class provides methods for constructing
 * vector containers compatible with the bases. Therefore, we use the `StdTraits`
 * type definition for different kinds of vectors. The vectors can be accessed using
 * the static method `vector<T>(index)` where `T` is the field-type to be used in the
 * vector container.
 **/
class BasisFactories
{
  using Traits = StdTraits;

public:
  static constexpr std::size_t size = 16;

  static auto basis(index_constant<0>)
  {
    return composite(
        power<3>(
          lagrange<2>(),
          flatLexicographic()),
        lagrange<1>(),
        blockedLexicographic()
      );
  }

  template <class T>
  static auto vector(index_constant<0>)
  {
    return Traits::PowerVector<Traits::DynamicVector<T>,2>{};
  }


  static auto basis(index_constant<1>)
  {
    return composite(
        power<3>(
          lagrange<2>(),
          flatInterleaved()),
        lagrange<1>(),
        blockedLexicographic()
      );
  }

  template <class T>
  static auto vector(index_constant<1>)
  {
    return Traits::PowerVector<Traits::DynamicVector<T>,2>{};
  }


  static auto basis(index_constant<2>)
  {
    return composite(
        power<3>(
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
      Traits::PowerVector<Traits::DynamicVector<T>, 3>,
      Traits::DynamicVector<T>
      >;
    return Vector{};
  }


  static auto basis(index_constant<3>)
  {
    return composite(
        power<3>(
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
      Traits::DynamicVector<Dune::FieldVector<T, 3>>,
      Traits::DynamicVector<T>
      >;
    return Vector{};
  }


  static auto basis(index_constant<4>)
  {
    return composite(
        power<3>(
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


  static auto basis(index_constant<5>)
  {
    return composite(
        power<3>(
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
      Traits::PowerVector<Traits::DynamicVector<T>, 3>,
      Traits::PowerVector<Traits::DynamicVector<T>, 1>
      >;
    return Vector{};
  }


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
    using Vector = Traits::PowerVector<Traits::DynamicVector<T>, 3>;
    return Vector{};
  }


  static auto basis(index_constant<7>)
  {
    return composite(
        power<1>(
          power<3>(
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
      Traits::PowerVector<Traits::PowerVector<Traits::DynamicVector<T>, 3>, 1>,
      Traits::DynamicVector<T>
      >;
    return Vector{};
  }


  static auto basis(index_constant<8>)
  {
    return composite(
        power<1>(
          power<3>(
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
      Traits::PowerVector<Traits::PowerVector<Traits::DynamicVector<T>, 3>, 1>,
      Traits::PowerVector<Traits::DynamicVector<T>, 1>
      >;
    return Vector{};
  }


  static auto basis(index_constant<9>)
  {
    return composite(
        power<2>(lagrange<2>()),
        composite(
          power<3>(
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
        Traits::DynamicVector<FieldVector<T, 3>>,
        Traits::DynamicVector<T>
        >
      >;
    return Vector{};
  }


  static auto basis(index_constant<10>)
  {
    return composite(
        lagrange<2>(),
        power<2>(
          lagrange<2>(),
          blockedLexicographic()),
        composite(
          power<3>(
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
        Traits::PowerVector<Traits::DynamicVector<T>, 3>,
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
};

}}} // end namespace Dune::Functions::BasisFactory
