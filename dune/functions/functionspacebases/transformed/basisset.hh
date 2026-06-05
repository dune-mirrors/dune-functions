// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_BASISSET_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_BASISSET_HH

#include <cassert>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/rangeutilities.hh>

#include <dune/functions/functionspacebases/transformed/derivative.hh>
#include <dune/functions/functionspacebases/transformed/derivativetraits.hh>

namespace Dune::Functions {

/**
 * \brief Derivative traits for transformations acting on the complete basis set.
 *
 * The public derivative range and the precomputation buffer are both the
 * corresponding reference local-basis types.  This is suitable for
 * non-affine-equivalent scalar elements where the element transformation mixes
 * basis functions but does not change the value/Jacobian/Hessian range type.
 */
namespace Impl {

template<class LocalBasis, class Derivative>
struct BasisSetDerivativeRange;

template<class LocalBasis>
struct BasisSetDerivativeRange<LocalBasis,Derivatives::Value>
{
  using type = typename LocalBasis::Traits::RangeType;
};

template<class LocalBasis>
struct BasisSetDerivativeRange<LocalBasis,Derivatives::Jacobian>
{
  using type = typename LocalBasis::Traits::JacobianType;
};

template<class LocalBasis>
struct BasisSetDerivativeRange<LocalBasis,Derivatives::Hessian>
{
  using type = typename LocalBasis::Traits::HessianType;
};

template<class InputRange, class Derivative, class LocalBasis, class... Stages>
struct PipelineOutputRange
{
  using type = InputRange;
};

template<class InputRange, class Derivative, class LocalBasis, class Stage, class... Stages>
struct PipelineOutputRange<InputRange,Derivative,LocalBasis,Stage,Stages...>
{
  using StageOutputRange = typename Stage::template OutputRange<Derivative,LocalBasis,InputRange>;
  using type = typename PipelineOutputRange<StageOutputRange,Derivative,LocalBasis,Stages...>::type;
};

} // namespace Impl

template<class LocalBasis, class Context>
struct BasisSetDerivativeTraits
{
  template<class Derivative>
  using DerivativeRange = typename Impl::BasisSetDerivativeRange<LocalBasis,Derivative>::type;

  template<class Derivative>
  using PrecomputeBuffer = std::vector<DerivativeRange<Derivative>>;
};

/**
 * \brief Pipeline stage that linearly mixes complete basis-function sets.
 *
 * The wrapped policy has to provide bind(context) and
 *
 * \code{.cpp}
 * template<class In, class Out>
 * void transformBasisSet(In const& in, Out& out) const;
 * \endcode
 */
template<class Transformation>
class BasisSetTransformationStage
{
  public:
    template<class Derivative, class LocalBasis, class InputRange>
    using OutputRange = InputRange;

    template<class Context>
    void bind(Context const& context)
    {
      transformation_.bind(context);
    }

    template<class Derivative, class LocalBasis, class InputRange, class OutputRange>
    void transform(Derivative,
                   LocalBasis const&,
                   typename LocalBasis::Traits::DomainType const&,
                   std::vector<InputRange> const& in,
                   std::vector<OutputRange>& out) const
    {
      transformation_.transformBasisSet(in,out);
    }

    Transformation const& transformation() const
    {
      return transformation_;
    }

  private:
    Transformation transformation_;
};

/**
 * \brief Pipeline stage applying affine scalar derivative chain rules.
 *
 * Values are passed through unchanged. Jacobians and Hessians are transformed
 * from reference to physical coordinates.
 */
template<class Geometry>
class AffineScalarDerivativePullbackStage
{
  public:
    template<class Derivative, class LocalBasis, class InputRange>
    using OutputRange = InputRange;

    template<class Context>
    void bind(Context const& context)
    {
      geometry_ = &context.geometry();
    }

    template<class LocalBasis, class InputRange, class OutputRange>
    void transform(Derivatives::Value,
                   LocalBasis const&,
                   typename LocalBasis::Traits::DomainType const&,
                   std::vector<InputRange> const& in,
                   std::vector<OutputRange>& out) const
    {
      out.assign(in.begin(),in.end());
    }

    template<class LocalBasis, class InputRange, class OutputRange>
    void transform(Derivatives::Jacobian,
                   LocalBasis const&,
                   typename LocalBasis::Traits::DomainType const& x,
                   std::vector<InputRange> const& in,
                   std::vector<OutputRange>& out) const
    {
      assert(!!geometry_);
      out.resize(in.size());
      auto&& Jinv = geometry_->jacobianInverse(x);
      for (auto i : Dune::range(in.size()))
        out[i] = in[i] * Jinv;
    }

    template<class LocalBasis, class InputRange, class OutputRange>
    void transform(Derivatives::Hessian,
                   LocalBasis const&,
                   typename LocalBasis::Traits::DomainType const& x,
                   std::vector<InputRange> const& in,
                   std::vector<OutputRange>& out) const
    {
      assert(!!geometry_);
      out.resize(in.size());
      auto&& Jinv = geometry_->jacobianInverse(x);
      auto&& JinvT = Jinv.transposed();
      for (auto i : Dune::range(in.size()))
        out[i] = JinvT * in[i] * Jinv;
    }

  private:
    Geometry const* geometry_ = nullptr;
};

/**
 * \brief Apply a sequence of typed transformation stages.
 *
 * Reference values are evaluated once in precompute(). Each stage then maps a
 * vector of its input range to a vector of its declared OutputRange. This allows
 * stages to be composed even when a later transformation changes the public
 * derivative range type.
 */
template<class... Stages>
class TransformationPipeline
{
    static_assert(sizeof...(Stages) > 0);

    template<class LocalBasis, class Context>
    struct PipelineTraits
    {
      template<class Derivative>
      using ReferenceRange = typename Impl::BasisSetDerivativeRange<LocalBasis,Derivative>::type;

      template<class Derivative>
      using PrecomputeBuffer = std::vector<ReferenceRange<Derivative>>;

      template<class Derivative>
      using DerivativeRange = typename Impl::PipelineOutputRange<
        ReferenceRange<Derivative>,
        Derivative,
        LocalBasis,
        Stages...>::type;
    };

  public:
    template<class LocalBasis, class Context>
    using Traits = PipelineTraits<LocalBasis,Context>;

    template<class Context>
    void bind(Context const& context)
    {
      std::apply([&](auto&... stage) { (stage.bind(context), ...); }, stages_);
    }

    template<class LocalBasis, class Out>
    void precompute(Derivatives::Value,
                    LocalBasis const& localBasis,
                    typename LocalBasis::Traits::DomainType const& x,
                    Out& out) const
    {
      localBasis.evaluateFunction(x,out);
    }

    template<class LocalBasis, class Out>
    void precompute(Derivatives::Jacobian,
                    LocalBasis const& localBasis,
                    typename LocalBasis::Traits::DomainType const& x,
                    Out& out) const
    {
      localBasis.evaluateJacobian(x,out);
    }

    template<class LocalBasis, class Out>
    void precompute(Derivatives::Hessian,
                    LocalBasis const& localBasis,
                    typename LocalBasis::Traits::DomainType const& x,
                    Out& out) const
    {
      localBasis.evaluateHessian(x,out);
    }

    template<class Derivative, class LocalBasis, class In, class Out>
    void finalize(Derivative derivative,
                  LocalBasis const& localBasis,
                  typename LocalBasis::Traits::DomainType const& x,
                  In const& in,
                  Out& out) const
    {
      transform<0>(derivative,localBasis,x,in,out);
    }

  private:
    template<std::size_t i, class Derivative, class LocalBasis, class InputRange, class OutputRange>
    void transform(Derivative derivative,
                   LocalBasis const& localBasis,
                   typename LocalBasis::Traits::DomainType const& x,
                   std::vector<InputRange> const& in,
                   std::vector<OutputRange>& out) const
    {
      using Stage = std::tuple_element_t<i,std::tuple<Stages...>>;
      using IntermediateRange = typename Stage::template OutputRange<Derivative,LocalBasis,InputRange>;

      if constexpr (i+1 == sizeof...(Stages))
        std::get<i>(stages_).transform(derivative,localBasis,x,in,out);
      else {
        std::vector<IntermediateRange> intermediate;
        std::get<i>(stages_).transform(derivative,localBasis,x,in,intermediate);
        transform<i+1>(derivative,localBasis,x,intermediate,out);
      }
    }

    std::tuple<Stages...> stages_;
};

} // end namespace Dune::Functions

#endif
