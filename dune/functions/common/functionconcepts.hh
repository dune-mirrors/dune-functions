// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_COMMON_FUNCTIONCONCEPT_HH
#define DUNE_FUNCTIONS_COMMON_FUNCTIONCONCEPT_HH

#include <concepts>

#include <dune/common/typelist.hh>
#include <dune/common/concept.hh>
#include <dune/functions/common/signature.hh>
#include <dune/functions/gridfunctions/localderivativetraits.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>

namespace Dune {
namespace Functions {
namespace Concept {



// Callable concept ############################################################
/**
 * \brief Concept objects that can be called with given argument list
 *
 * \ingroup FunctionConcepts
 *
 * \tparam Args Argument list for function call
 */
template<class F, class... Args>
concept Callable = requires(F f, Args... args) {
  f(args...);
};


// Function concept ############################################################
/**
 * \brief Concept for a function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 */
template<class F, class Range, class Domain>
concept Function = requires(F f, Domain arg) {
  { f(arg) } -> std::convertible_to<Range>;
};

/// Check if F models the Function concept with given signature \ingroup FunctionConcepts
template<class F, class Range, class Domain>
constexpr auto isFunction()
{ return std::bool_constant<Concept::Function<F, Range, Domain> >{}; }

/// Check if f models the Function concept with given signature \ingroup FunctionConcepts
template<class F, class Range, class Domain>
constexpr auto isFunction(const F& f, SignatureTag<Range(Domain)>)
{
  return std::bool_constant<Function<F,Range,Domain>>{};
}



// DifferentiableFunction concept ##############################################
/**
 * \brief Concept for a differentiable function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * The derivative range is derived from the provided \p DerivativeTraits
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 * \tparam DerivativeTraits Traits class for computation of derivative range
 */
template<class F, class Range, class Domain, template<class> class DerivativeTraits = DefaultDerivativeTraits>
concept DifferentiableFunction = Function<F,Range,Domain>
&& requires (F f) {
  { derivative(f) } -> Function<typename DerivativeTraits<Range(Domain)>::Range, Domain>;
};

/// Check if F models the DifferentiableFunction concept with given signature \ingroup FunctionConcepts
template <class F, class Range, class Domain,
          template<class> class DerivativeTraits = DefaultDerivativeTraits>
constexpr auto isDifferentiableFunction()
{ return std::bool_constant<Concept::DifferentiableFunction<F, Range, Domain, DerivativeTraits> >{}; }

/// Check if f models the DifferentiableFunction concept with given signature \ingroup FunctionConcepts
template<class F, class Range, class Domain>
constexpr auto isDifferentiableFunction(const F& f, SignatureTag<Range(Domain)>)
{
  return std::bool_constant<DifferentiableFunction<F,Range,Domain>>{};
}


// LocalFunction concept ##############################################
/**
 * \brief Concept for a local function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 * \tparam LocalContext The local context this function is defined on
 */
template<class F, class Range, class Domain, class LocalContext>
concept LocalFunction = Function<F,Range,Domain>
&& requires(F& f, LocalContext ctx) {
  f.bind(ctx);
  f.unbind();
} && requires(const F& f) {
  { f.bound() } -> std::convertible_to<bool>;
  { f. localContext() } -> std::convertible_to<LocalContext>;
};

/// Check if F models the LocalFunction concept with given signature and local context \ingroup FunctionConcepts
template<class F, class Range, class Domain, class LocalContext>
constexpr auto isLocalFunction()
{ return std::bool_constant<Concept::LocalFunction<F, Range, Domain, LocalContext> >{}; }



// DifferentiableLocalFunction concept ##############################################
/**
 * \brief Concept for a differentiable local function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * The derivative range is derived from the provided \p DerivativeTraits
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 * \tparam LocalContext The local context this function is defined on
 * \tparam DerivativeTraits Traits class for computation of derivative range
 */
template <class F, class Range, class Domain, class LocalContext,
          template<class> class DerivativeTraits = DefaultDerivativeTraits>
concept DifferentiableLocalFunction =
   DifferentiableFunction<F, Range, Domain, DerivativeTraits>
&& LocalFunction<F, Range, Domain, LocalContext>
&& requires(F& f, LocalContext ctx) {
    f.bind(ctx);
    f.unbind();
} && requires(const F& f) {
    { f.localContext() } -> std::convertible_to<LocalContext>;
};

/// Check if F models the DifferentiableLocalFunction concept with given signature and local context \ingroup FunctionConcepts
template <class F, class Range, class Domain, class LocalContext,
          template<class> class DerivativeTraits = DefaultDerivativeTraits>
constexpr auto isDifferentiableLocalFunction()
{ return std::bool_constant<Concept::DifferentiableLocalFunction<F, Range, Domain, LocalContext, DerivativeTraits> >{}; }



// EntitySet concept ##############################################
/**
 * \brief Concept for an entity set for a \ref Concept::GridFunction<Range(Domain), EntitySet, DerivativeTraits>
 *
 * \ingroup FunctionConcepts
 *
 * This describes the set of entities on which a grid function
 * can be localized.
 *
 */
template<class E>
concept EntitySet = requires(E entitySet) {
  typename E::Element;
  typename E::LocalCoordinate;
  typename E::GlobalCoordinate;
};

/// Check if F models the GridFunction concept with given signature and entity set \ingroup FunctionConcepts
template<class E>
constexpr auto isEntitySet()
{ return std::bool_constant<Concept::EntitySet<E>>{}; }

// GridFunction concept ##############################################
/**
 * \brief Concept for a grid function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 * \tparam EntitySet Set of entities on which the function can be localized
 */
template<class F, class Range, class Domain, class ES>
concept GridFunction = Function<F,Range,Domain>
&& requires(const F& f) {
  { localFunction(f) } -> LocalFunction<Range, typename ES::LocalCoordinate, typename ES::Element>;
  { f.entitySet() } -> std::convertible_to<ES>;

  requires EntitySet<ES>;
  requires std::convertible_to<typename ES::GlobalCoordinate, Domain>;
};
/// Check if F models the GridFunction concept with given signature and entity set \ingroup FunctionConcepts
template<class F, class Range, class Domain, class EntitySet>
constexpr auto isGridFunction()
{ return std::bool_constant<Concept::GridFunction<F, Range, Domain, EntitySet> >{}; }


namespace Imp {

  template <class ES, template<class> class DerivativeTraits>
  struct LocalDerivativeTraits
  {
    template <class R>
    using Traits = typename Dune::Functions::LocalDerivativeTraits<ES, DerivativeTraits>::template Traits<R>;
  };

} // end namespace Imp

// DifferentiableGridFunction concept ##############################################
/**
 * \brief Concept for a differentiable grid function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * The derivative range is derived from the provided \p DerivativeTraits.
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 * \tparam EntitySet Set of entities on which the function can be localized
 * \tparam DerivativeTraits Traits class for computation of derivative range
 */
template <class F, class Range, class Domain, class ES,
          template<class> class DerivativeTraits = DefaultDerivativeTraits>
concept DifferentiableGridFunction =
   DifferentiableFunction<F, Range, Domain, DerivativeTraits>
&& GridFunction<F, Range, Domain, ES>
&& requires(const F& f) {
  { localFunction(f) } -> DifferentiableLocalFunction<
          Range, typename ES::LocalCoordinate, typename ES::Element,
          Imp::LocalDerivativeTraits<ES, DerivativeTraits>::template Traits>;
};

/// Check if F models the DifferentiableGridFunction concept with given signature and entity set \ingroup FunctionConcepts
template <class F, class Range, class Domain, class EntitySet,
          template<class> class DerivativeTraits = DefaultDerivativeTraits>
constexpr auto isDifferentiableGridFunction()
{ return std::bool_constant<Concept::DifferentiableGridFunction<F, Range, Domain, EntitySet, DerivativeTraits> >{}; }


// GridViewFunction concept ##############################################
/**
 * \brief Concept for a grid view function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * This exactly the \ref Concept::GridFunction<Range(Domain), EntitySet>
 * concept with a \ref GridViewEntitySet as \p EntitySet.
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 * \tparam GridView GridView on which the function can be localized
 */
template<class F, class Range, class Domain, class GridView>
concept GridViewFunction = GridFunction<F, Range, Domain, GridViewEntitySet<GridView,0> >;

/// Check if F models the GridViewFunction concept with given signature \ingroup FunctionConcepts
template<class F, class Range, class Domain, class GridView>
constexpr auto isGridViewFunction()
{ return std::bool_constant<Concept::GridViewFunction<F, Range, Domain, GridView> >{}; }


// DifferentiableGridViewFunction concept ##############################################
/**
 * \brief Concept for a differentiable grid view function mapping \p Domain to \p Range
 *
 * \ingroup FunctionConcepts
 *
 * This exactly the \ref Concept::GridFunction<Range(Domain), EntitySet, DerivativeTraits>
 * concept with a \ref GridViewEntitySet as \p EntitySet.
 *
 * \tparam Domain Domain type
 * \tparam Range Range type
 * \tparam GridView GridView on which the function can be localized
 * \tparam DerivativeTraits Traits class for computation of derivative range
 */
template <class F, class Range, class Domain, class GridView,
          template<class> class DerivativeTraits = DefaultDerivativeTraits>
concept DifferentiableGridViewFunction
  = DifferentiableGridFunction<F, Range, Domain, GridViewEntitySet<GridView,0>, DerivativeTraits>;

/// Check if F models the DifferentiableGridViewFunction concept with given signature \ingroup FunctionConcepts
template <class F, class Range, class Domain, class GridView,
          template<class> class DerivativeTraits = DefaultDerivativeTraits>
constexpr bool isDifferentiableGridViewFunction()
{ return std::bool_constant<Concept::DifferentiableGridViewFunction<F, Range, Domain, GridView, DerivativeTraits> >{}; }


}}} // namespace Dune::Functions::Concept

#endif // DUNE_FUNCTIONS_COMMON_FUNCTIONCONCEPT_HH
