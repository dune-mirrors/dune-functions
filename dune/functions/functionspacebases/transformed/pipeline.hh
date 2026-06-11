// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_PIPELINE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_PIPELINE_HH

#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/functions/functionspacebases/transformed/concepts.hh>
#include <dune/functions/functionspacebases/transformed/referenceevaluation.hh>

namespace Dune::Functions {

namespace Impl {

template<class InputRange, class Derivative, class LocalBasis, class Context,
         class... Stages>
struct PipelineOutputRange
{
  using type = InputRange;
};

template<class InputRange, class Derivative, class LocalBasis, class Context,
         class Stage, class... Stages>
struct PipelineOutputRange<InputRange,Derivative,LocalBasis,Context,Stage,Stages...>
{
  using StageOutputRange = typename Stage::template OutputRange<
    Derivative,LocalBasis,Context,InputRange>;
  using type = typename PipelineOutputRange<
    StageOutputRange,Derivative,LocalBasis,Context,Stages...>::type;
};

template<class InputRange, class Derivative, class LocalBasis, class Context,
         class... Stages>
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

} // namespace Impl

/**
 * \brief Adapt a pointwise range transformation to an evaluation-pipeline stage.
 *
 * The wrapped policy transforms one basis-function result at a time and may
 * change its scalar, vector, or tensor type. The adapter applies this operation
 * to the complete vector of basis-function results.
 */
template<class Transformation>
class RangeTransformationStage
{
  public:
    template<class Derivative, class LocalBasis, class Context, class InputRange>
    using OutputRange = typename Transformation::template OutputRange<
      Derivative,LocalBasis,Context,InputRange>;

    RangeTransformationStage() = default;

    explicit RangeTransformationStage(Transformation transformation)
      : transformation_(std::move(transformation))
    {}

    template<class Context>
    void bind(Context const& context)
    {
      transformation_.bind(context);
    }

    template<class Derivative, class LocalBasis, class InputRange, class OutputRange>
    void transform(Derivative derivative,
                   LocalBasis const& localBasis,
                   typename LocalBasis::Traits::DomainType const& x,
                   std::vector<InputRange> const& in,
                   std::vector<OutputRange>& out) const
      requires requires(Transformation const& transformation) {
        transformation.transform(derivative,localBasis,x,in[0],out[0]);
      }
    {
      out.resize(in.size());
      for (std::size_t i = 0; i < in.size(); ++i)
        transformation_.transform(derivative,localBasis,x,in[i],out[i]);
    }

    Transformation const& transformation() const
    {
      return transformation_;
    }

  private:
    Transformation transformation_;
};

/**
 * \brief Apply an ordered sequence of typed basis-evaluation stages.
 *
 * Reference values are evaluated once in precompute(). Finalize applies all
 * stages in order. Intermediate vectors are part of the caller-owned buffer,
 * while the cached reference vector is treated as immutable.
 */
template<class ReferenceEvaluator, class Context, class... Stages>
class BasicBasisEvaluationPipeline
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
        mutable typename Impl::PipelineIntermediateBuffers<
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

    BasicBasisEvaluationPipeline() = default;

    explicit BasicBasisEvaluationPipeline(Stages... stages)
      : stages_(std::move(stages)...)
    {}

    BasicBasisEvaluationPipeline(ReferenceEvaluator evaluator, Stages... stages)
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
                  In const& in,
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
using BasisEvaluationPipeline = BasicBasisEvaluationPipeline<
  ReferenceEvaluation<>,Context,Stages...>;

} // end namespace Dune::Functions

#endif
