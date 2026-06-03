// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_LOCALFINITEELEMENT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_LOCALFINITEELEMENT_HH

#include <array>
#include <cassert>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/geometry/type.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>

#include <dune/functions/functionspacebases/transformed/concepts.hh>
#include <dune/functions/functionspacebases/transformed/derivative.hh>

namespace Dune::Functions {

namespace TransformedLocalFiniteElementLocalBasis {
  struct Reference {};
  struct Physical {};
}

/**
 * \brief Local-basis adapter driven by a transformation policy.
 *
 * This adapter wraps a reference local basis and applies a transformation policy
 * to evaluate basis functions on the physical element.  The evaluation is split
 * into a precompute/finalize protocol, allowing element-dependent data to be reused
 * across multiple evaluation points.
 *
 * \tparam LocalBasis Type of the reference local basis to transform (e.g., Lagrange
 *   basis on the reference element).
 * \tparam Context Type of the element context providing geometry information
 *   (must satisfy LocalFiniteElementBindContext).
 * \tparam Transformation Type of the transformation policy implementing the
 *   precompute/finalize protocol (must satisfy LocalBasisTransformation for the
 *   desired derivatives).
 *
 * The class is designed to work with different transformation policies, such as
 * covariant Piola transforms for H(curl) spaces or contravariant Piola transforms
 * for H(div) spaces.
 * The transformation policy owns the element-dependent logic (e.g., Jacobian
 * computations).
 */
template<class LocalBasis, class Context, class Transformation>
class TransformedLocalBasis
{
    using DerivativeTraits = typename Transformation::template Traits<LocalBasis,Context>;

  public:
    //! Buffer type used by precompute() for the quantity selected by D.
    template<class D>
    using PrecomputeBuffer = typename DerivativeTraits::template PrecomputeBuffer<D>;

    //! Public output range type for the quantity selected by D.
    template<class D>
    using DerivativeRange = typename DerivativeTraits::template DerivativeRange<D>;

    //! Reference-element coordinate type.
    using Domain = typename LocalBasis::Traits::DomainType;

    //! Output range type for Derivatives::Value.
    using Range = DerivativeRange<Derivatives::Value>;

    //! LocalBasis traits for compatibility with the dune-localfunctions interface.
    using Traits = typename LocalBasis::Traits;

    /**
     * \brief Default constructor creating an unbound basis.
     *
     * The basis needs to be bound to a LocalBasis and a Context before use.
     */
    TransformedLocalBasis() = default;

    /**
     * \brief Constructor with transformation policy.
     *
     * \param transformation The transformation policy to use for basis evaluation.
     */
    explicit TransformedLocalBasis(Transformation transformation)
      : transformation_(std::move(transformation))
    {}

    /**
     * \brief Bind to a reference local basis.
     *
     * \param localBasis The reference local basis to transform.
     */
    void bind(LocalBasis const& localBasis)
    {
      localBasis_ = &localBasis;
    }

    /**
     * \brief Bind to an element context.
     *
     * This forwards the context to the transformation policy, which can extract
     * any element-dependent data it needs (e.g., the geometry for Piola transforms).
     *
     * \param context The element context containing geometry information.
     */
    void bind(Context const& context)
    {
      transformation_.bind(context);
    }

    /**
     * \brief Return the number of shape functions.
     *
     * \return Number of shape functions in the underlying local basis.
     */
    std::size_t size() const
    {
      assert(!!localBasis_);
      return localBasis_->size();
    }

    /**
     * \brief Return the polynomial order of the basis.
     *
     * \return Maximum polynomial degree of shape functions in the underlying local basis.
     */
    int order() const
    {
      assert(!!localBasis_);
      return localBasis_->order();
    }

    /**
     * \brief Precompute element-dependent data for a derivative.
     *
     * This is the first stage of the split evaluation protocol.  The precomputed
     * data is stored in \p out and later passed to finalize().
     *
     * \tparam D The derivative tag (e.g., Derivatives::Value, Derivatives::Jacobian).
     * \param derivative Tag selecting which derivative to compute.
     * \param x Point in the reference element where to evaluate.
     * \param[out] out Buffer to store precomputed data.
     */
    template<class D>
      requires Concept::LocalBasisTransformation<Transformation,LocalBasis,Context,D>
    void precompute(D derivative,
                    Domain const& x,
                    PrecomputeBuffer<D>& out) const
    {
      assert(!!localBasis_);
      transformation_.precompute(derivative, *localBasis_, x, out);
    }

    /**
     * \brief Finalize evaluation using precomputed data.
     *
     * This is the second stage of the split evaluation protocol.  It uses the
     * precomputed buffer from precompute() to produce the final transformed values.
     *
     * \tparam D The derivative tag (e.g., Derivatives::Value, Derivatives::Jacobian).
     * \param derivative Tag selecting which derivative to compute.
     * \param x Point in the reference element where to evaluate.
     * \param in Precomputed data buffer from precompute().
     * \param[out] out Vector to store the evaluated derivative for each shape function.
     */
    template<class D>
      requires Concept::LocalBasisTransformation<Transformation,LocalBasis,Context,D>
    void finalize(D derivative,
                  Domain const& x,
                  PrecomputeBuffer<D> const& in,
                  std::vector<DerivativeRange<D>>& out) const
    {
      assert(!!localBasis_);
      transformation_.finalize(derivative, *localBasis_, x, in, out);
    }

    /**
     * \brief Evaluate a derivative at a point (combined precompute and finalize).
     *
     * This is a convenience method that performs both precomputation and finalization
     * in one call.  For performance-critical code, prefer the split protocol.
     *
     * \tparam D The derivative tag (e.g., Derivatives::Value, Derivatives::Jacobian).
     * \param derivative Tag selecting which derivative to compute.
     * \param x Point in the reference element where to evaluate.
     * \param[out] out Vector to store the evaluated derivative for each shape function.
     */
    template<class D>
      requires Concept::LocalBasisTransformation<Transformation,LocalBasis,Context,D>
    void evaluate(D derivative,
                  Domain const& x,
                  std::vector<DerivativeRange<D>>& out) const
    {
      PrecomputeBuffer<D> precomputed;
      precompute(derivative, x, precomputed);
      finalize(derivative, x, precomputed, out);
    }

    //! Evaluate all shape functions using the classic LocalBasis interface.
    void evaluateFunction(Domain const& x, std::vector<DerivativeRange<Derivatives::Value>>& out) const
      requires (Concept::LocalBasisTransformation<Transformation,LocalBasis,Context,Derivatives::Value>)
    {
      evaluate(Derivatives::Value{}, x, out);
    }

    //! Evaluate all shape function Jacobians using the classic LocalBasis interface.
    void evaluateJacobian(Domain const& x, std::vector<DerivativeRange<Derivatives::Jacobian>>& out) const
      requires (Concept::LocalBasisTransformation<Transformation,LocalBasis,Context,Derivatives::Jacobian>)
    {
      evaluate(Derivatives::Jacobian{}, x, out);
    }

    /**
     * \brief Access the transformation policy.
     *
     * \return Const reference to the transformation policy.
     */
    Transformation const& transformation() const
    {
      return transformation_;
    }

  private:
    LocalBasis const* localBasis_ = nullptr;
    Transformation transformation_;
};

/**
 * \brief Local interpolation adapter driven by a transformation policy.
 *
 * This adapter wraps a reference local interpolation and applies a transformation
 * policy to handle interpolation of global-valued functions.  The transformation
 * policy provides a pullback mechanism that converts global-valued functions into
 * the local-valued functions expected by the reference interpolation.
 *
 * \tparam LocalInterpolation Type of the reference local interpolation to use.
 * \tparam Context Type of the element context providing geometry information
 *   (must satisfy LocalFiniteElementBindContext).
 * \tparam Transformation Type of the transformation policy providing the pullback
 *   (must satisfy LocalInterpolationTransformation).
 *
 * For example, a Piola transformation would pull back a physical vector field to
 * the reference element before passing it to the reference interpolation.
 */
template<class LocalInterpolation, class Context, class Transformation>
class TransformedLocalInterpolation
{
  public:
    /**
     * \brief Default constructor creating an unbound interpolation.
     */
    TransformedLocalInterpolation() = default;

    /**
     * \brief Constructor with transformation policy.
     *
     * \param transformation The transformation policy providing the pullback.
     */
    explicit TransformedLocalInterpolation(Transformation transformation)
      : transformation_(std::move(transformation))
    {}

    /**
     * \brief Bind to a reference local interpolation.
     *
     * \param localInterpolation The reference local interpolation to use.
     */
    void bind(LocalInterpolation const& localInterpolation)
    {
      localInterpolation_ = &localInterpolation;
    }

    /**
     * \brief Bind to an element context.
     *
     * \param context The element context containing geometry information.
     */
    void bind(Context const& context)
    {
      transformation_.bind(context);
    }

    /**
     * \brief Interpolate a global-valued function onto the transformed basis.
     *
     * The transformation policy pulls back the global-valued function \p f to
     * a local-valued function that the reference interpolation can handle.
     *
     * \tparam F Type of the function to interpolate.
     * \tparam C Coefficient type for the interpolation result.
     * \param f The function to interpolate.
     * \param[out] out Vector to store the interpolation coefficients.
     */
    template<class F, class C>
      requires Concept::LocalInterpolationTransformation<Transformation,Context,F>
    void interpolate(F const& f, std::vector<C>& out) const
    {
      assert(!!localInterpolation_);
      auto localValuedFunction = transformation_.localFunctionPullback(f);
      localInterpolation_->interpolate(localValuedFunction, out);
    }

    /**
     * \brief Access the transformation policy.
     *
     * \return Const reference to the transformation policy.
     */
    Transformation const& transformation() const
    {
      return transformation_;
    }

  private:
    LocalInterpolation const* localInterpolation_ = nullptr;
    Transformation transformation_;
};

/**
 * \brief Transformed local finite-element adapter.
 *
 * This class wraps a reference local finite element and exposes a transformed
 * local basis and interpolation.  It forwards coefficient access to the reference
 * local finite element, while basis evaluation and interpolation are handled
 * through transformation policies.
 *
 * The class supports two modes controlled by the LocalBasisTag template parameter:
 * - Physical (default): localBasis() returns the transformed physical basis, and
 *   localInterpolation() returns the transformed interpolation.
 * - Reference: localBasis() returns the reference basis, and localInterpolation()
 *   returns the reference interpolation.
 *
 * \tparam LocalFiniteElement Type of the reference local finite element to wrap.
 * \tparam Context Type of the element context providing geometry information
 *   (must satisfy LocalFiniteElementBindContext).
 * \tparam Transformation Type of the transformation policy for basis and interpolation
 *   (must satisfy LocalBasisTransformation and LocalInterpolationTransformation
 *   as needed).
 * \tparam LocalBasisTag Tag to select whether localBasis() returns the reference or
 *   physical basis. Use TransformedLocalFiniteElementLocalBasis::Reference or
 *   ::Physical.
 *
 * Example usage:
 * \code
 * using LFE = ...; // reference local finite element type
 * using Context = ...; // element context type
 * using PiolaTransform = ...; // Piola transformation policy
 *
 * TransformedLocalFiniteElement<LFE, Context, PiolaTransform> tlfe(piola);
 * tlfe.bind(lfe, context);
 * auto const& basis = tlfe.localBasis(); // returns transformed physical basis
 * \endcode
 */
template<class LocalFiniteElement,
         class Context,
         class Transformation,
         class LocalBasisTag = TransformedLocalFiniteElementLocalBasis::Physical>
class TransformedLocalFiniteElement
{
    using ReferenceLocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;
    using ReferenceLocalInterpolation = typename LocalFiniteElement::Traits::LocalInterpolationType;
    using ReferenceLocalCoefficients = typename LocalFiniteElement::Traits::LocalCoefficientsType;

  public:
    //! Type of the reference local basis.
    using ReferenceBasis = ReferenceLocalBasis;

    //! Type of the transformed physical local basis.
    using PhysicalBasis = TransformedLocalBasis<ReferenceLocalBasis,Context,Transformation>;

    //! Type exposed by localBasis() for compatibility with existing bases.
    using LocalBasis = std::conditional_t<
      std::is_same_v<LocalBasisTag,TransformedLocalFiniteElementLocalBasis::Reference>,
      ReferenceBasis,
      PhysicalBasis>;

    //! Type of the transformed physical local interpolation.
    using PhysicalLocalInterpolation = TransformedLocalInterpolation<ReferenceLocalInterpolation,Context,Transformation>;

    //! Type exposed by localInterpolation() for compatibility with existing bases.
    using LocalInterpolation = std::conditional_t<
      std::is_same_v<LocalBasisTag,TransformedLocalFiniteElementLocalBasis::Reference>,
      ReferenceLocalInterpolation,
      PhysicalLocalInterpolation>;

    //! Export number types, dimensions, etc.
    using Traits = LocalFiniteElementTraits<LocalBasis,ReferenceLocalCoefficients,LocalInterpolation>;

    /**
     * \brief Default constructor creating an unbound finite element.
     */
    TransformedLocalFiniteElement() = default;

    /**
     * \brief Constructor with transformation policy.
     *
     * \param transformation The transformation policy to use for basis and interpolation.
     */
    explicit TransformedLocalFiniteElement(Transformation transformation)
      : basis_(std::move(transformation))
      , physicalInterpolation_(basis_.transformation())
    {}

    /**
     * \brief Bind to a reference local finite element.
     *
     * This binds the internal basis and interpolation to the corresponding
     * components of the reference local finite element.
     *
     * \param lfe The reference local finite element to wrap.
     */
    void bind(LocalFiniteElement const& lfe)
    {
      lfe_ = &lfe;
      basis_.bind(lfe.localBasis());
      physicalInterpolation_.bind(lfe.localInterpolation());
    }

    /**
     * \brief Bind to an element context.
     *
     * This forwards the context to the transformation policy.
     *
     * \param context The element context containing geometry information.
     */
    void bind(Context const& context)
    {
      basis_.bind(context);
      physicalInterpolation_.bind(context);
      type_ = context.type();
    }

    /**
     * \brief Convenience method to bind both reference LFE and context at once.
     *
     * \param lfe The reference local finite element to wrap.
     * \param context The element context containing geometry information.
     */
    void bind(LocalFiniteElement const& lfe, Context const& context)
    {
      bind(lfe);
      bind(context);
    }

    /**
     * \brief Access the local basis.
     *
     * Depending on the LocalBasisTag template parameter, this returns either
     * the reference basis (for Reference tag) or the transformed physical basis
     * (for Physical tag, the default).
     *
     * \return The local basis (reference or physical).
     */
    LocalBasis const& localBasis() const
    {
      if constexpr (std::is_same_v<LocalBasisTag,TransformedLocalFiniteElementLocalBasis::Reference>)
        return referenceBasis();
      else
        return physicalBasis();
    }

    /**
     * \brief Access the transformed physical basis.
     *
     * \return Const reference to the transformed physical local basis.
     */
    PhysicalBasis const& physicalBasis() const
    {
      return basis_;
    }

    /**
     * \brief Access the reference basis.
     *
     * \return Const reference to the underlying reference local basis.
     */
    ReferenceBasis const& referenceBasis() const
    {
      assert(!!lfe_);
      return lfe_->localBasis();
    }

    /**
     * \brief Access the local coefficients.
     *
     * These are forwarded from the underlying reference local finite element.
     *
     * \return Const reference to the local coefficients.
     */
    ReferenceLocalCoefficients const& localCoefficients() const
    {
      assert(!!lfe_);
      return lfe_->localCoefficients();
    }

    /**
     * \brief Access the local interpolation.
     *
     * Depending on the LocalBasisTag template parameter, this returns either
     * the reference interpolation (for Reference tag) or the transformed physical
     * interpolation (for Physical tag, the default).
     *
     * \return The local interpolation (reference or physical).
     */
    LocalInterpolation const& localInterpolation() const
    {
      if constexpr (std::is_same_v<LocalBasisTag,TransformedLocalFiniteElementLocalBasis::Reference>) {
        assert(!!lfe_);
        return lfe_->localInterpolation();
      }
      else
        return physicalInterpolation_;
    }

    /**
     * \brief Return the number of shape functions.
     *
     * \return Number of shape functions.
     */
    std::size_t size() const
    {
      return basis_.size();
    }

    /**
     * \brief Return the geometry type of the bound element.
     *
     * \return The GeometryType of the element from the bound context.
     */
    GeometryType type() const
    {
      return type_;
    }

  private:
    LocalFiniteElement const* lfe_ = nullptr;
    PhysicalBasis basis_;
    PhysicalLocalInterpolation physicalInterpolation_;
    GeometryType type_ = {};
};

} // end namespace Dune::Functions

#endif
