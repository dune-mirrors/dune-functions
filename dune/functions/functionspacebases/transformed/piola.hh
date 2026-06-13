// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_PIOLA_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_PIOLA_HH

#include <cassert>
#include <type_traits>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/referencehelper.hh>
#include <dune/common/std/no_unique_address.hh>

#include <dune/functions/common/densevectorview.hh>
#include <dune/functions/functionspacebases/transformed/derivative.hh>
#include <dune/functions/functionspacebases/transformed/referenceevaluation.hh>

namespace Dune::Functions {

namespace Impl {

template<class LocalBasis, class Geometry, class Derivative>
struct ContravariantPiolaRange;

template<class LocalBasis, class Geometry>
struct ContravariantPiolaRange<LocalBasis,Geometry,Derivatives::Value>
{
  using type = FieldVector<typename LocalBasis::Traits::RangeFieldType,Geometry::coorddimension>;
};

template<class LocalBasis, class Geometry>
struct ContravariantPiolaRange<LocalBasis,Geometry,Derivatives::Divergence>
{
  using type = typename LocalBasis::Traits::RangeFieldType;
};

template<class LocalBasis, class Geometry>
struct ContravariantPiolaRange<LocalBasis,Geometry,Derivatives::Jacobian>
{
  using type = FieldMatrix<
    typename LocalBasis::Traits::RangeFieldType,
    Geometry::coorddimension,
    Geometry::coorddimension>;
};

template<class LocalBasis, class Geometry, class Derivative>
struct CovariantPiolaRange;

template<class LocalBasis, class Geometry, int dimension>
struct CovariantCurlRange
{};

template<class LocalBasis, class Geometry>
struct CovariantCurlRange<LocalBasis,Geometry,2>
{
  using type = typename LocalBasis::Traits::RangeFieldType;
};

template<class LocalBasis, class Geometry>
struct CovariantCurlRange<LocalBasis,Geometry,3>
{
  using type = FieldVector<typename LocalBasis::Traits::RangeFieldType,Geometry::coorddimension>;
};

template<class LocalBasis, class Geometry>
struct CovariantPiolaRange<LocalBasis,Geometry,Derivatives::Value>
{
  using type = FieldVector<
    typename LocalBasis::Traits::RangeFieldType,
    Geometry::coorddimension>;
};

template<class LocalBasis, class Geometry>
struct CovariantPiolaRange<LocalBasis,Geometry,Derivatives::Curl>
  : CovariantCurlRange<LocalBasis,Geometry,Geometry::mydimension>
{
};

template<class LocalBasis, class Geometry>
struct CovariantPiolaRange<LocalBasis,Geometry,Derivatives::Jacobian>
  : ContravariantPiolaRange<LocalBasis,Geometry,Derivatives::Jacobian>
{};

template<class LocalBasis, class Geometry, template<class,class,class> class Range>
struct PiolaDerivativeTraits
{
  template<class Derivative>
  using PrecomputeBuffer = std::vector<
    typename ReferenceLocalBasisEvaluator::template Range<Derivative,LocalBasis>>;

  template<class Derivative>
  using DerivativeRange = typename Range<LocalBasis,Geometry,Derivative>::type;
};

} // namespace Impl

/**
 * \brief Transformation policy for the contravariant Piola map.
 *
 * This policy preserves normal traces and is the canonical value
 * transformation for H(div)-conforming finite elements.  Its natural
 * derivative is Derivatives::Divergence.
 */
template<class Geometry>
class ContravariantPiolaTransformation
{
  public:
    template<class LocalBasis, class Context>
    using Traits = Impl::PiolaDerivativeTraits<
      LocalBasis,Geometry,Impl::ContravariantPiolaRange>;

    template<class Context>
    void bind(Context const& context)
    {
      geometry_ = &context.geometry();
    }

    template<class Function>
    class LocalValuedFunction
    {
      public:
        LocalValuedFunction(Function f, Geometry const& geometry)
          : f_(std::move(f))
          , geometry_(&geometry)
        {}

        template<class LocalCoordinate>
        auto operator()(LocalCoordinate const& x) const
        {
          auto globalValue = Dune::resolveRef(f_)(x);
          using Field = std::remove_cvref_t<decltype(Impl::DenseVectorView(globalValue)[0])>;
          FieldVector<Field,Geometry::mydimension> localValue;

          auto jacobianInverse = geometry_->jacobianInverse(x);
          auto globalValueDenseVector = Impl::DenseVectorView(globalValue);
          jacobianInverse.mv(globalValueDenseVector, localValue);
          localValue *= geometry_->integrationElement(x);

          return localValue;
        }

      private:
        Function f_;
        Geometry const* geometry_;
    };

    template<class Function>
    auto localFunctionPullback(Function f) const
    {
      assert(!!geometry_);
      return LocalValuedFunction<Function>(std::move(f), *geometry_);
    }

    template<class LocalBasis, class Out>
    void precompute(Derivatives::Value,
                    LocalBasis const& localBasis,
                    typename LocalBasis::Traits::DomainType const& x,
                    Out& out) const
    {
      evaluator_.evaluate(Derivatives::Value{},localBasis,x,out);
    }

    template<class LocalBasis, class Out>
    void precompute(Derivatives::Divergence,
                    LocalBasis const& localBasis,
                    typename LocalBasis::Traits::DomainType const& x,
                    Out& out) const
    {
      evaluator_.evaluate(Derivatives::Divergence{},localBasis,x,out);
    }

    template<class LocalBasis, class In, class Out>
    void finalize(Derivatives::Value,
                  LocalBasis const&,
                  typename LocalBasis::Traits::DomainType const& x,
                  In const& in,
                  Out& out) const
    {
      assert(!!geometry_);
      out.resize(in.size());

      auto jacobian = geometry_->jacobian(x);
      auto integrationElement = geometry_->integrationElement(x);

      for (auto i : Dune::range(in.size())) {
        auto tmpDenseVector = Impl::DenseVectorView(in[i]);
        auto outDenseVector = Impl::DenseVectorView(out[i]);
        jacobian.mv(tmpDenseVector, outDenseVector);
        out[i] /= integrationElement;
      }
    }

    template<class LocalBasis>
      requires (Geometry::mydimension == Geometry::coorddimension)
    void evaluateCompatibilityJacobian(
        LocalBasis const& localBasis,
        typename LocalBasis::Traits::DomainType const& x,
        std::vector<typename LocalBasis::Traits::JacobianType>& out) const
    {
      assert(!!geometry_);
      if (!geometry_->affine())
        DUNE_THROW(NotImplemented,
          "Contravariant Piola Jacobians require affine geometry");

      localBasis.evaluateJacobian(x,out);
      auto jacobianTransposed = geometry_->jacobianTransposed(x);
      auto integrationElement = geometry_->integrationElement(x);
      for (auto& jacobian : out) {
        auto referenceJacobian = jacobian;
        jacobian = 0;
        for (std::size_t k = 0; k < jacobian.M(); ++k)
          for (std::size_t l = 0; l < referenceJacobian.N(); ++l)
            for (auto&& [entry,j] : sparseRange(jacobianTransposed[l]))
              jacobian[j][k] += entry * referenceJacobian[l][k];
        jacobian /= integrationElement;
      }
    }

    template<class LocalBasis, class In, class Out>
    void finalize(Derivatives::Divergence,
                  LocalBasis const&,
                  typename LocalBasis::Traits::DomainType const& x,
                  In const& in,
                  Out& out) const
    {
      assert(!!geometry_);
      out.resize(in.size());

      auto integrationElement = geometry_->integrationElement(x);

      for (auto i : Dune::range(in.size()))
        out[i] = in[i] / integrationElement;
    }

  private:
    DUNE_NO_UNIQUE_ADDRESS ReferenceLocalBasisEvaluator evaluator_;
    Geometry const* geometry_ = nullptr;
};

/**
 * \brief Transformation policy for the covariant Piola map.
 *
 * This policy preserves tangential traces and is the canonical value
 * transformation for H(curl)-conforming finite elements.  Its natural
 * derivative is Derivatives::Curl.
 */
template<class Geometry>
class CovariantPiolaTransformation
{
  public:
    template<class LocalBasis, class Context>
    using Traits = Impl::PiolaDerivativeTraits<
      LocalBasis,Geometry,Impl::CovariantPiolaRange>;

    template<class Context>
    void bind(Context const& context)
    {
      geometry_ = &context.geometry();
    }

    template<class Function>
    class LocalValuedFunction
    {
      public:
        LocalValuedFunction(Function f, Geometry const& geometry)
          : f_(std::move(f))
          , geometry_(&geometry)
        {}

        template<class LocalCoordinate>
        auto operator()(LocalCoordinate const& x) const
        {
          auto globalValue = Dune::resolveRef(f_)(x);
          using Field = std::remove_cvref_t<decltype(Impl::DenseVectorView(globalValue)[0])>;
          FieldVector<Field,Geometry::mydimension> localValue;

          auto jacobianTransposed = geometry_->jacobianTransposed(x);
          auto globalValueDenseVector = Impl::DenseVectorView(globalValue);
          jacobianTransposed.mv(globalValueDenseVector, localValue);

          return localValue;
        }

      private:
        Function f_;
        Geometry const* geometry_;
    };

    template<class Function>
    auto localFunctionPullback(Function f) const
    {
      assert(!!geometry_);
      return LocalValuedFunction<Function>(std::move(f), *geometry_);
    }

    template<class LocalBasis, class Out>
    void precompute(Derivatives::Value,
                    LocalBasis const& localBasis,
                    typename LocalBasis::Traits::DomainType const& x,
                    Out& out) const
    {
      evaluator_.evaluate(Derivatives::Value{},localBasis,x,out);
    }

    template<class LocalBasis, class Out>
    void precompute(Derivatives::Curl,
                    LocalBasis const& localBasis,
                    typename LocalBasis::Traits::DomainType const& x,
                    Out& out) const
    {
      evaluator_.evaluate(Derivatives::Curl{},localBasis,x,out);
    }

    template<class LocalBasis, class In, class Out>
    void finalize(Derivatives::Value,
                  LocalBasis const&,
                  typename LocalBasis::Traits::DomainType const& x,
                  In const& in,
                  Out& out) const
    {
      assert(!!geometry_);
      out.resize(in.size());

      auto jacobianInverseTransposed = geometry_->jacobianInverseTransposed(x);

      for (auto i : Dune::range(in.size())) {
        auto tmpDenseVector = Impl::DenseVectorView(in[i]);
        auto outDenseVector = Impl::DenseVectorView(out[i]);
        jacobianInverseTransposed.mv(tmpDenseVector, outDenseVector);
      }
    }

    template<class LocalBasis, class In, class Out>
      requires (Geometry::mydimension == 2 || Geometry::mydimension == 3)
    void finalize(Derivatives::Curl,
                  LocalBasis const&,
                  typename LocalBasis::Traits::DomainType const& x,
                  In const& in,
                  Out& out) const
    {
      assert(!!geometry_);
      out.resize(in.size());

      auto integrationElement = geometry_->integrationElement(x);

      if constexpr (Geometry::mydimension == 2) {
        for (auto i : Dune::range(in.size()))
          out[i] = in[i] / integrationElement;
      }
      else {
        static_assert(Geometry::mydimension == 3,
          "CovariantPiolaTransformation::finalize(Curl) supports only intrinsic dimension 2 or 3.");

        auto jacobian = geometry_->jacobian(x);
        for (auto i : Dune::range(in.size())) {
          auto referenceCurlDenseVector = Impl::DenseVectorView(in[i]);
          auto outDenseVector = Impl::DenseVectorView(out[i]);
          jacobian.mv(referenceCurlDenseVector, outDenseVector);
          out[i] /= integrationElement;
        }
      }
    }

    template<class LocalBasis>
      requires (Geometry::mydimension == Geometry::coorddimension)
    void evaluateCompatibilityJacobian(
        LocalBasis const& localBasis,
        typename LocalBasis::Traits::DomainType const& x,
        std::vector<typename LocalBasis::Traits::JacobianType>& out) const
    {
      assert(!!geometry_);
      if (!geometry_->affine())
        DUNE_THROW(NotImplemented,
          "Covariant Piola Jacobians require affine geometry");

      localBasis.evaluateJacobian(x,out);
      auto jacobianInverseTransposed = geometry_->jacobianInverseTransposed(x);
      for (auto& jacobian : out) {
        auto referenceJacobian = jacobian;
        jacobian = 0;
        for (std::size_t j = 0; j < jacobian.N(); ++j)
          for (std::size_t k = 0; k < jacobian.M(); ++k)
            for (auto&& [entry,l] : sparseRange(jacobianInverseTransposed[j]))
              jacobian[j][k] += entry * referenceJacobian[l][k];
      }
    }

  private:
    DUNE_NO_UNIQUE_ADDRESS ReferenceLocalBasisEvaluator evaluator_;
    Geometry const* geometry_ = nullptr;
};

} // end namespace Dune::Functions

#endif
