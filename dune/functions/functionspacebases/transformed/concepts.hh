// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_CONCEPTS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_CONCEPTS_HH

#include <concepts>
#include <cstddef>
#include <utility>
#include <vector>

#include <dune/geometry/type.hh>

namespace Dune::Functions::Concept {

/**
 * \brief Structural interface of a reference local finite element.
 */
template<class LocalFiniteElement>
concept ReferenceLocalFiniteElement =
  requires(LocalFiniteElement const& finiteElement) {
    typename LocalFiniteElement::Traits;
    typename LocalFiniteElement::Traits::LocalBasisType;
    typename LocalFiniteElement::Traits::LocalCoefficientsType;
    typename LocalFiniteElement::Traits::LocalInterpolationType;

    { finiteElement.localBasis() }
      -> std::same_as<typename LocalFiniteElement::Traits::LocalBasisType const&>;
    { finiteElement.localCoefficients() }
      -> std::same_as<typename LocalFiniteElement::Traits::LocalCoefficientsType const&>;
    { finiteElement.localInterpolation() }
      -> std::same_as<typename LocalFiniteElement::Traits::LocalInterpolationType const&>;
  };

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
 * \brief Structural requirements shared by all local-basis transformations.
 */
template<class Transformation, class LocalBasis, class Context>
concept LocalBasisTransformationPolicy =
  LocalFiniteElementBindContext<Context> &&
  requires(Transformation transformation, Context const& context) {
    typename Transformation::template Traits<LocalBasis,Context>;
    transformation.bind(context);
  };

/**
 * \brief Structural requirements shared by interpolation transformations.
 */
template<class Transformation, class Context>
concept InterpolationTransformationPolicy =
  LocalFiniteElementBindContext<Context> &&
  requires(Transformation transformation, Context const& context) {
    transformation.bind(context);
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
  LocalBasisTransformationPolicy<Transformation,LocalBasis,Context> &&
  requires(Transformation transformation,
           LocalBasis const& localBasis,
           Context const& context,
           typename LocalBasis::Traits::DomainType const& x,
           Derivative derivative,
           typename Transformation::template Traits<LocalBasis,Context>::template PrecomputeBuffer<Derivative>& precomputed,
           std::vector<typename Transformation::template Traits<LocalBasis,Context>::template DerivativeRange<Derivative>>& out)
  {
    typename Transformation::template Traits<LocalBasis,Context>;
    typename Transformation::template Traits<LocalBasis,Context>::template PrecomputeBuffer<Derivative>;
    typename Transformation::template Traits<LocalBasis,Context>::template DerivativeRange<Derivative>;

    transformation.bind(context);
    std::as_const(transformation).precompute(derivative, localBasis, x, precomputed);
    std::as_const(transformation).finalize(derivative, localBasis, x, std::as_const(precomputed), out);
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
  InterpolationTransformationPolicy<Transformation,Context> &&
  requires(Transformation transformation,
           Context const& context,
           Function const& f)
  {
    transformation.bind(context);
    std::as_const(transformation).localFunctionPullback(f);
  };

/**
 * \brief One typed stage in a local-basis transformation pipeline.
 *
 * A stage may change the range type.  Its output therefore depends on the
 * derivative, local basis, bind context, and input range.
 */
template<class Stage, class Derivative, class LocalBasis, class Context, class InputRange>
concept TransformationStage =
  LocalFiniteElementBindContext<Context> &&
  requires(Stage stage,
           Context const& context,
           Derivative derivative,
           LocalBasis const& localBasis,
           typename LocalBasis::Traits::DomainType const& x,
           std::vector<InputRange> const& in,
           std::vector<typename Stage::template OutputRange<
             Derivative,LocalBasis,Context,InputRange>>& out)
  {
    typename Stage::template OutputRange<Derivative,LocalBasis,Context,InputRange>;
    stage.bind(context);
    std::as_const(stage).transform(derivative,localBasis,x,in,out);
  };

/**
 * \brief Transformed local basis with split precompute/finalize evaluation.
 *
 * The concept describes the local-basis interface used by transformed finite elements.
 * As for transformations, it is checked for one derivative tag at a time.
 */
template<class Basis, class Derivative>
concept TransformedLocalBasis =
  requires(Basis const& basis,
           typename Basis::Domain const& x,
           Derivative derivative,
           typename Basis::template PrecomputeBuffer<Derivative>& precomputed,
           std::vector<typename Basis::template DerivativeRange<Derivative>>& out)
  {
    typename Basis::Domain;
    typename Basis::Range;
    typename Basis::template PrecomputeBuffer<Derivative>;
    typename Basis::template DerivativeRange<Derivative>;

    { basis.size() } -> std::convertible_to<std::size_t>;
    { basis.order() } -> std::convertible_to<int>;

    basis.precompute(derivative, x, precomputed);
    basis.finalize(derivative, x, std::as_const(precomputed), out);
    basis.evaluate(derivative, x, out);
  };

} // end namespace Dune::Functions::Concept

#endif
