// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_CONCEPTS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_CONCEPTS_HH

#include <concepts>
#include <cstddef>
#include <vector>

#include <dune/geometry/type.hh>

namespace Dune::Functions::Concept {

/**
 * \brief Minimal bind context for local finite-element transformations.
 *
 * A bind context is the geometry-dependent state made available to local
 * finite-element transformations.  The minimal common denominator is an
 * element geometry and its GeometryType.  Transformations with additional
 * needs, such as face orientations or mesh-size data, should refine this
 * concept with their own requirements.
 */
template<class Context>
concept LocalFiniteElementBindContext =
  requires(Context const& context) {
    context.geometry();
    { context.type() } -> std::convertible_to<GeometryType>;
  };

/**
 * \brief Staged transformation for a selected derivative quantity.
 *
 * The transformation owns the element-dependent part of the evaluation
 * protocol for the derivative tag \p Derivative.  It maps reference-element
 * precomputation buffers to the public transformed range type.
 *
 * This concept is intentionally checked per derivative tag.  For example, an
 * H(div) transformation may support values, Jacobians, and divergences without
 * claiming support for Hessians.
 */
template<class Transformation, class LocalBasis, class Context, class Derivative>
concept LocalBasisTransformation =
  LocalFiniteElementBindContext<Context> &&
  requires(Transformation transformation,
           Transformation const& constTransformation,
           LocalBasis const& localBasis,
           Context const& context,
           typename LocalBasis::Traits::DomainType const& x,
           Derivative derivative,
           typename Transformation::template Traits<LocalBasis,Context>::template PrecomputeBuffer<Derivative>& precomputed,
           typename Transformation::template Traits<LocalBasis,Context>::template PrecomputeBuffer<Derivative> const& constPrecomputed,
           std::vector<typename Transformation::template Traits<LocalBasis,Context>::template DerivativeRange<Derivative>>& out)
  {
    typename Transformation::template Traits<LocalBasis,Context>;
    typename Transformation::template Traits<LocalBasis,Context>::template PrecomputeBuffer<Derivative>;
    typename Transformation::template Traits<LocalBasis,Context>::template DerivativeRange<Derivative>;

    transformation.bind(context);
    constTransformation.precompute(derivative, localBasis, x, precomputed);
    constTransformation.finalize(derivative, localBasis, x, constPrecomputed, out);
  };

/**
 * \brief Transformation with an inverse value map for local interpolation.
 *
 * A transformed local interpolation receives global-valued functions.  The
 * transformation policy has to pull those functions back to the local-valued
 * callable expected by the reference local interpolation.
 */
template<class Transformation, class Context, class Function>
concept LocalInterpolationTransformation =
  LocalFiniteElementBindContext<Context> &&
  requires(Transformation transformation,
           Transformation const& constTransformation,
           Context const& context,
           Function const& f)
  {
    transformation.bind(context);
    constTransformation.localFunctionPullback(f);
  };

/**
 * \brief Transformed local basis with split precompute/finalize evaluation.
 *
 * The concept describes the experimental local-basis interface used by
 * transformed finite elements.  As for transformations, it is checked for one
 * derivative tag at a time.
 */
template<class Basis, class Derivative>
concept StagedTransformedLocalBasis =
  requires(Basis const& basis,
           typename Basis::Domain const& x,
           Derivative derivative,
           typename Basis::template PrecomputeBuffer<Derivative>& precomputed,
           typename Basis::template PrecomputeBuffer<Derivative> const& constPrecomputed,
           std::vector<typename Basis::template DerivativeRange<Derivative>>& out)
  {
    typename Basis::Domain;
    typename Basis::Range;
    typename Basis::template PrecomputeBuffer<Derivative>;
    typename Basis::template DerivativeRange<Derivative>;

    { basis.size() } -> std::convertible_to<std::size_t>;
    { basis.order() } -> std::convertible_to<int>;

    basis.precompute(derivative, x, precomputed);
    basis.finalize(derivative, x, constPrecomputed, out);
    basis.evaluate(derivative, x, out);
  };

} // end namespace Dune::Functions::Concept

#endif
