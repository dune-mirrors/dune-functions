// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_LOCALFINITEELEMENT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_LOCALFINITEELEMENT_HH

#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

#include <dune/geometry/type.hh>

#include <dune/functions/functionspacebases/transformed/concepts.hh>
#include <dune/functions/functionspacebases/transformed/derivative.hh>

namespace Dune::Functions {

/**
 * \brief Draft local-basis adapter driven by a transformation policy.
 *
 * The adapter factors transformed basis evaluation into a reference local
 * basis and a transformation object.  The transformation object owns the
 * element-dependent logic and exposes the staged precompute/finalize protocol
 * checked by Concept::LocalBasisTransformation.
 *
 * This class is intentionally small and experimental.  It is meant to make the
 * intended structure visible while the transformation policy interface is still
 * being refined.
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

    TransformedLocalBasis() = default;

    explicit TransformedLocalBasis(Transformation transformation)
      : transformation_(std::move(transformation))
    {}

    void bind(LocalBasis const& localBasis)
    {
      localBasis_ = &localBasis;
    }

    void bind(Context const& context)
    {
      transformation_.bind(context);
    }

    std::size_t size() const
    {
      assert(!!localBasis_);
      return localBasis_->size();
    }

    int order() const
    {
      assert(!!localBasis_);
      return localBasis_->order();
    }

    template<class D>
      requires Concept::LocalBasisTransformation<Transformation,LocalBasis,Context,D>
    void precompute(D derivative,
                    Domain const& x,
                    PrecomputeBuffer<D>& out) const
    {
      assert(!!localBasis_);
      transformation_.precompute(derivative, *localBasis_, x, out);
    }

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

    Transformation const& transformation() const
    {
      return transformation_;
    }

  private:
    LocalBasis const* localBasis_ = nullptr;
    Transformation transformation_;
};

/**
 * \brief Draft transformed local finite-element adapter.
 *
 * The class wraps a reference local finite element and exposes a transformed
 * local basis.  Coefficients are forwarded to the reference finite element.
 * Interpolation is deliberately not modeled yet, because inverse value
 * transformations and transformed dual functionals need a separate policy
 * hook.
 */
template<class LocalFiniteElement, class Context, class Transformation>
class TransformedLocalFiniteElement
{
    using ReferenceLocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;

  public:
    //! Type of the transformed local basis exposed by localBasis() and basis().
    using Basis = TransformedLocalBasis<ReferenceLocalBasis,Context,Transformation>;

    TransformedLocalFiniteElement() = default;

    explicit TransformedLocalFiniteElement(Transformation transformation)
      : basis_(std::move(transformation))
    {}

    void bind(LocalFiniteElement const& lfe)
    {
      lfe_ = &lfe;
      basis_.bind(lfe.localBasis());
    }

    void bind(Context const& context)
    {
      basis_.bind(context);
      type_ = context.type();
    }

    void bind(LocalFiniteElement const& lfe, Context const& context)
    {
      bind(lfe);
      bind(context);
    }

    Basis const& localBasis() const
    {
      return basis_;
    }

    Basis const& basis() const
    {
      return basis_;
    }

    auto const& localCoefficients() const
    {
      assert(!!lfe_);
      return lfe_->localCoefficients();
    }

    std::size_t size() const
    {
      return basis_.size();
    }

    GeometryType type() const
    {
      return type_;
    }

  private:
    LocalFiniteElement const* lfe_ = nullptr;
    Basis basis_;
    GeometryType type_ = {};
};

} // end namespace Dune::Functions

#endif
