// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_LOCALFINITEELEMENT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_LOCALFINITEELEMENT_HH

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

    //! LocalBasis traits for compatibility with the dune-localfunctions interface.
    using Traits = typename LocalBasis::Traits;

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

    //! Evaluate all shape functions using the classic LocalBasis interface.
    void evaluateFunction(Domain const& x, std::vector<typename Traits::RangeType>& out) const
      requires Concept::LocalBasisTransformation<Transformation,LocalBasis,Context,Derivatives::Value>
    {
      evaluate(Derivatives::Value{}, x, out);
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
 * \brief Local interpolation adapter driven by the same transformation policy.
 *
 * The transformation policy provides the pullback used to convert a
 * global-valued function into the local-valued function expected by the
 * reference local interpolation.
 */
template<class LocalInterpolation, class Context, class Transformation>
class TransformedLocalInterpolation
{
  public:
    TransformedLocalInterpolation() = default;

    explicit TransformedLocalInterpolation(Transformation transformation)
      : transformation_(std::move(transformation))
    {}

    void bind(LocalInterpolation const& localInterpolation)
    {
      localInterpolation_ = &localInterpolation;
    }

    void bind(Context const& context)
    {
      transformation_.bind(context);
    }

    template<class F, class C>
      requires Concept::LocalInterpolationTransformation<Transformation,Context,F>
    void interpolate(F const& f, std::vector<C>& out) const
    {
      assert(!!localInterpolation_);
      auto localValuedFunction = transformation_.localFunctionPullback(f);
      localInterpolation_->interpolate(localValuedFunction, out);
    }

    Transformation const& transformation() const
    {
      return transformation_;
    }

  private:
    LocalInterpolation const* localInterpolation_ = nullptr;
    Transformation transformation_;
};

/**
 * \brief Draft transformed local finite-element adapter.
 *
 * The class wraps a reference local finite element and exposes a transformed
 * local basis.  Coefficients are forwarded to the reference finite element,
 * and interpolation is adapted through the same transformation policy.
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
    using LocalBasis = std::conditional_t<std::is_same_v<LocalBasisTag,TransformedLocalFiniteElementLocalBasis::Reference>,
                                          ReferenceBasis,
                                          PhysicalBasis>;

    //! Type of the transformed physical local interpolation.
    using PhysicalLocalInterpolation = TransformedLocalInterpolation<ReferenceLocalInterpolation,Context,Transformation>;

    //! Type exposed by localInterpolation() for compatibility with existing bases.
    using LocalInterpolation = std::conditional_t<std::is_same_v<LocalBasisTag,TransformedLocalFiniteElementLocalBasis::Reference>,
                                                  ReferenceLocalInterpolation,
                                                  PhysicalLocalInterpolation>;

    //! Export number types, dimensions, etc.
    using Traits = LocalFiniteElementTraits<LocalBasis,ReferenceLocalCoefficients,LocalInterpolation>;

    TransformedLocalFiniteElement() = default;

    explicit TransformedLocalFiniteElement(Transformation transformation)
      : basis_(std::move(transformation))
      , physicalInterpolation_(basis_.transformation())
    {}

    void bind(LocalFiniteElement const& lfe)
    {
      lfe_ = &lfe;
      basis_.bind(lfe.localBasis());
      physicalInterpolation_.bind(lfe.localInterpolation());
    }

    void bind(Context const& context)
    {
      basis_.bind(context);
      physicalInterpolation_.bind(context);
      type_ = context.type();
    }

    void bind(LocalFiniteElement const& lfe, Context const& context)
    {
      bind(lfe);
      bind(context);
    }

    LocalBasis const& localBasis() const
    {
      if constexpr (std::is_same_v<LocalBasisTag,TransformedLocalFiniteElementLocalBasis::Reference>)
        return referenceBasis();
      else
        return physicalBasis();
    }

    PhysicalBasis const& physicalBasis() const
    {
      return basis_;
    }

    ReferenceBasis const& referenceBasis() const
    {
      assert(!!lfe_);
      return lfe_->localBasis();
    }

    ReferenceBasis const& referenceLocalBasis() const
    {
      return referenceBasis();
    }

    auto const& localCoefficients() const
    {
      assert(!!lfe_);
      return lfe_->localCoefficients();
    }

    LocalInterpolation const& localInterpolation() const
    {
      if constexpr (std::is_same_v<LocalBasisTag,TransformedLocalFiniteElementLocalBasis::Reference>) {
        assert(!!lfe_);
        return lfe_->localInterpolation();
      }
      else
        return physicalInterpolation_;
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
    PhysicalBasis basis_;
    PhysicalLocalInterpolation physicalInterpolation_;
    GeometryType type_ = {};
};

} // end namespace Dune::Functions

#endif
