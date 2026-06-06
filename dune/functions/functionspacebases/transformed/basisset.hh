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

#include <dune/common/exceptions.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/functions/functionspacebases/transformed/concepts.hh>
#include <dune/functions/functionspacebases/transformed/derivative.hh>
#include <dune/functions/functionspacebases/transformed/referenceevaluation.hh>

namespace Dune::Functions {

namespace Impl {

template<class InputRange, class Derivative, class LocalBasis, class Context, class... Stages>
struct PipelineOutputRange
{
  using type = InputRange;
};

template<class InputRange, class Derivative, class LocalBasis, class Context, class Stage, class... Stages>
struct PipelineOutputRange<InputRange,Derivative,LocalBasis,Context,Stage,Stages...>
{
  using StageOutputRange = typename Stage::template OutputRange<Derivative,LocalBasis,Context,InputRange>;
  using type = typename PipelineOutputRange<StageOutputRange,Derivative,LocalBasis,Context,Stages...>::type;
};

template<class Derivative, class LocalBasis, class Geometry, class InputRange>
struct GeometryDerivativeOutputRange
{
  using type = InputRange;
};

template<class InputRange, class Derivative, class LocalBasis, class Context, class... Stages>
struct PipelineIntermediateBuffers;

template<class InputRange, class Derivative, class LocalBasis, class Context, class Stage>
struct PipelineIntermediateBuffers<InputRange,Derivative,LocalBasis,Context,Stage>
{
  using type = std::tuple<>;
};

template<class InputRange, class Derivative, class LocalBasis, class Context,
         class Stage, class NextStage, class... Stages>
struct PipelineIntermediateBuffers<
  InputRange,Derivative,LocalBasis,Context,Stage,NextStage,Stages...>
{
  using StageOutputRange = typename Stage::template OutputRange<
    Derivative,LocalBasis,Context,InputRange>;
  using Tail = typename PipelineIntermediateBuffers<
    StageOutputRange,Derivative,LocalBasis,Context,NextStage,Stages...>::type;
  using type = decltype(std::tuple_cat(
    std::declval<std::tuple<std::vector<StageOutputRange>>>(),
    std::declval<Tail>()));
};

template<class LocalBasis, class Geometry, class InputRange>
struct GeometryDerivativeOutputRange<Derivatives::Jacobian,LocalBasis,Geometry,InputRange>
{
  using type = FieldMatrix<typename LocalBasis::Traits::RangeFieldType,
                           LocalBasis::Traits::dimRange,
                           Geometry::coorddimension>;
};

template<class LocalBasis, class Geometry, class InputRange>
struct GeometryDerivativeOutputRange<Derivatives::Gradient,LocalBasis,Geometry,InputRange>
{
  static_assert(LocalBasis::Traits::dimRange == 1);
  using type = FieldVector<typename LocalBasis::Traits::RangeFieldType,Geometry::coorddimension>;
};

template<class LocalBasis, class Geometry, class InputRange>
struct GeometryDerivativeOutputRange<Derivatives::Laplacian,LocalBasis,Geometry,InputRange>
{
  using type = typename LocalBasis::Traits::RangeType;
};

} // namespace Impl

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
    template<class Derivative, class LocalBasis, class Context, class InputRange>
    using OutputRange = InputRange;

    BasisSetTransformationStage() = default;

    explicit BasisSetTransformationStage(Transformation transformation)
      : transformation_(std::move(transformation))
    {}

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
class GeometryDerivativePullbackStage
{
  public:
    template<class Derivative, class LocalBasis, class Context, class InputRange>
    using OutputRange = typename Impl::GeometryDerivativeOutputRange<
      Derivative,LocalBasis,Geometry,InputRange>::type;

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
    void transform(Derivatives::Gradient,
                   LocalBasis const&,
                   typename LocalBasis::Traits::DomainType const& x,
                   std::vector<InputRange> const& in,
                   std::vector<OutputRange>& out) const
    {
      static_assert(LocalBasis::Traits::dimRange == 1);
      assert(!!geometry_);
      out.resize(in.size());
      auto&& JinvT = geometry_->jacobianInverseTransposed(x);
      for (auto i : Dune::range(in.size()))
        JinvT.mv(in[i][0],out[i]);
    }

    template<class LocalBasis, class InputRange, class OutputRange>
      requires (Geometry::mydimension == Geometry::coorddimension)
    void transform(Derivatives::Hessian,
                   LocalBasis const&,
                   typename LocalBasis::Traits::DomainType const& x,
                   std::vector<InputRange> const& in,
                   std::vector<OutputRange>& out) const
    {
      static_assert(Geometry::mydimension == Geometry::coorddimension,
        "Physical Hessians are currently supported only for codimension-zero geometries.");
      assert(!!geometry_);
      if (!geometry_->affine())
        DUNE_THROW(NotImplemented,
          "Physical Hessians are currently supported only for affine geometries");
      out.resize(in.size());
      auto&& Jinv = geometry_->jacobianInverse(x);
      auto&& JinvT = Jinv.transposed();
      for (auto i : Dune::range(in.size()))
        out[i] = JinvT * in[i] * Jinv;
    }

    template<class LocalBasis, class InputRange, class OutputRange>
      requires (Geometry::mydimension == Geometry::coorddimension)
    void transform(Derivatives::Laplacian,
                   LocalBasis const&,
                   typename LocalBasis::Traits::DomainType const& x,
                   std::vector<InputRange> const& in,
                   std::vector<OutputRange>& out) const
    {
      static_assert(Geometry::mydimension == Geometry::coorddimension,
        "Physical Laplacians are currently supported only for codimension-zero geometries.");
      assert(!!geometry_);
      if (!geometry_->affine())
        DUNE_THROW(NotImplemented,
          "Physical Laplacians are currently supported only for affine geometries");

      out.resize(in.size());
      auto&& Jinv = geometry_->jacobianInverse(x);
      auto&& JinvT = Jinv.transposed();
      for (auto i : Dune::range(in.size())) {
        auto hessian = JinvT * in[i] * Jinv;
        out[i] = {};
        for (auto j : Dune::range(Geometry::coorddimension))
          out[i][0] += hessian[j][j];
      }
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
 *
 * \tparam Context Bind context shared by all stages.
 * \tparam Stages Ordered transformation-stage types.
 */
template<class ReferenceEvaluator, class Context, class... Stages>
class BasicTransformationPipeline
{
    static_assert(sizeof...(Stages) > 0);

    template<class LocalBasis, class BoundContext>
    struct PipelineTraits
    {
      static_assert(std::is_same_v<Context,BoundContext>);

      template<class Derivative>
      using ReferenceRange = typename ReferenceEvaluator::template Range<Derivative,LocalBasis>;

      template<class Derivative>
      struct PipelineBuffer
      {
        std::vector<ReferenceRange<Derivative>> reference;
        typename Impl::PipelineIntermediateBuffers<
          ReferenceRange<Derivative>,Derivative,LocalBasis,BoundContext,Stages...>::type intermediate;
      };

      template<class Derivative>
      using PrecomputeBuffer = PipelineBuffer<Derivative>;

      template<class Derivative>
      using DerivativeRange = typename Impl::PipelineOutputRange<
        ReferenceRange<Derivative>,
        Derivative,
        LocalBasis,
        BoundContext,
        Stages...>::type;
    };

  public:
    template<class LocalBasis, class BoundContext>
    using Traits = PipelineTraits<LocalBasis,BoundContext>;

    BasicTransformationPipeline() = default;

    explicit BasicTransformationPipeline(Stages... stages)
      : stages_(std::move(stages)...)
    {}

    BasicTransformationPipeline(ReferenceEvaluator evaluator, Stages... stages)
      : stages_(std::move(stages)...)
      , evaluator_(std::move(evaluator))
    {}

    template<class BoundContext>
    void bind(BoundContext const& context)
    {
      static_assert(std::is_same_v<Context,BoundContext>);
      std::apply([&](auto&... stage) { (stage.bind(context), ...); }, stages_);
    }

    template<class Derivative, class LocalBasis, class Out>
    void precompute(Derivative derivative,
                    LocalBasis const& localBasis,
                    typename LocalBasis::Traits::DomainType const& x,
                    Out& out) const
    {
      evaluator_.evaluate(derivative,localBasis,x,out.reference);
    }

    template<class Derivative, class LocalBasis, class In, class Out>
    void finalize(Derivative derivative,
                  LocalBasis const& localBasis,
                  typename LocalBasis::Traits::DomainType const& x,
                  In& in,
                  Out& out) const
    {
      transform<0>(derivative,localBasis,x,in.reference,in.intermediate,out);
    }

  private:
    template<std::size_t i, class Derivative, class LocalBasis,
             class InputRange, class IntermediateBuffers, class OutputRange>
    void transform(Derivative derivative,
                   LocalBasis const& localBasis,
                   typename LocalBasis::Traits::DomainType const& x,
                   std::vector<InputRange> const& in,
                   IntermediateBuffers& intermediateBuffers,
                   std::vector<OutputRange>& out) const
    {
      using Stage = std::tuple_element_t<i,std::tuple<Stages...>>;
      static_assert(Concept::TransformationStage<
        Stage,Derivative,LocalBasis,Context,InputRange>);

      if constexpr (i+1 == sizeof...(Stages))
        std::get<i>(stages_).transform(derivative,localBasis,x,in,out);
      else {
        auto& intermediate = std::get<i>(intermediateBuffers);
        std::get<i>(stages_).transform(derivative,localBasis,x,in,intermediate);
        transform<i+1>(
          derivative,localBasis,x,intermediate,intermediateBuffers,out);
      }
    }

    std::tuple<Stages...> stages_;
    ReferenceEvaluator evaluator_;
};

template<class Context, class... Stages>
using TransformationPipeline = BasicTransformationPipeline<
  ReferenceLocalBasisEvaluator,Context,Stages...>;

} // end namespace Dune::Functions

#endif
