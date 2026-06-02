// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_PIOLA_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_PIOLA_HH

#include <cassert>
#include <vector>

#include <dune/common/rangeutilities.hh>

#include <dune/functions/common/densevectorview.hh>
#include <dune/functions/functionspacebases/transformed/derivative.hh>
#include <dune/functions/functionspacebases/transformed/derivativetraits.hh>

namespace Dune::Functions {

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
    using Traits = StandardDerivativeTraits<LocalBasis,Geometry>;

    template<class Context>
    void bind(Context const& context)
    {
      geometry_ = &context.geometry();
    }

    template<class Function>
    class LocalValuedFunction
    {
      public:
        LocalValuedFunction(Function const& f, Geometry const& geometry)
          : f_(&f)
          , geometry_(&geometry)
        {}

        template<class LocalCoordinate>
        auto operator()(LocalCoordinate const& x) const
        {
          auto globalValue = (*f_)(x);
          auto localValue = globalValue;

          auto jacobianInverse = geometry_->jacobianInverse(x);
          auto globalValueDenseVector = Impl::DenseVectorView(globalValue);
          auto localValueDenseVector = Impl::DenseVectorView(localValue);
          jacobianInverse.mv(globalValueDenseVector, localValueDenseVector);
          localValue *= geometry_->integrationElement(x);

          return localValue;
        }

      private:
        Function const* f_;
        Geometry const* geometry_;
    };

    template<class Function>
    auto localFunctionPullback(Function const& f) const
    {
      assert(!!geometry_);
      return LocalValuedFunction<Function>(f, *geometry_);
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
    void precompute(Derivatives::Divergence,
                    LocalBasis const& localBasis,
                    typename LocalBasis::Traits::DomainType const& x,
                    Out& out) const
    {
      localBasis.evaluateJacobian(x,out);
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
        auto tmp = in[i];
        auto tmpDenseVector = Impl::DenseVectorView(tmp);
        auto outDenseVector = Impl::DenseVectorView(out[i]);
        jacobian.mv(tmpDenseVector, outDenseVector);
        out[i] /= integrationElement;
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

      for (auto i : Dune::range(in.size())) {
        out[i] = {};
        for (auto j : Dune::range(Geometry::coorddimension))
          out[i] += in[i][j][j];
        out[i] /= integrationElement;
      }
    }

  private:
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
    using Traits = StandardDerivativeTraits<LocalBasis,Geometry>;

    template<class Context>
    void bind(Context const& context)
    {
      geometry_ = &context.geometry();
    }

    template<class Function>
    class LocalValuedFunction
    {
      public:
        LocalValuedFunction(Function const& f, Geometry const& geometry)
          : f_(&f)
          , geometry_(&geometry)
        {}

        template<class LocalCoordinate>
        auto operator()(LocalCoordinate const& x) const
        {
          auto globalValue = (*f_)(x);
          auto localValue = globalValue;

          auto jacobianTransposed = geometry_->jacobianTransposed(x);
          auto globalValueDenseVector = Impl::DenseVectorView(globalValue);
          auto localValueDenseVector = Impl::DenseVectorView(localValue);
          jacobianTransposed.mv(globalValueDenseVector, localValueDenseVector);

          return localValue;
        }

      private:
        Function const* f_;
        Geometry const* geometry_;
    };

    template<class Function>
    auto localFunctionPullback(Function const& f) const
    {
      assert(!!geometry_);
      return LocalValuedFunction<Function>(f, *geometry_);
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
    void precompute(Derivatives::Curl,
                    LocalBasis const& localBasis,
                    typename LocalBasis::Traits::DomainType const& x,
                    Out& out) const
    {
      localBasis.evaluateJacobian(x,out);
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
        auto tmp = in[i];
        auto tmpDenseVector = Impl::DenseVectorView(tmp);
        auto outDenseVector = Impl::DenseVectorView(out[i]);
        jacobianInverseTransposed.mv(tmpDenseVector, outDenseVector);
      }
    }

    template<class LocalBasis, class In, class Out>
    void finalize(Derivatives::Curl,
                  LocalBasis const&,
                  typename LocalBasis::Traits::DomainType const& x,
                  In const& in,
                  Out& out) const
    {
      assert(!!geometry_);
      out.resize(in.size());

      auto integrationElement = geometry_->integrationElement(x);

      if constexpr (Geometry::coorddimension == 2) {
        for (auto i : Dune::range(in.size())) {
          out[i] = in[i][1][0] - in[i][0][1];
          out[i] /= integrationElement;
        }
      }
      else {
        static_assert(Geometry::coorddimension == 3,
          "CovariantPiolaTransformation::finalize(Curl) supports only dimension 2 or 3.");

        auto jacobian = geometry_->jacobian(x);
        using CurlRange = typename Traits<LocalBasis,void>::template DerivativeRange<Derivatives::Curl>;

        for (auto i : Dune::range(in.size())) {
          CurlRange referenceCurl;
          referenceCurl[0] = in[i][2][1] - in[i][1][2];
          referenceCurl[1] = in[i][0][2] - in[i][2][0];
          referenceCurl[2] = in[i][1][0] - in[i][0][1];

          auto referenceCurlDenseVector = Impl::DenseVectorView(referenceCurl);
          auto outDenseVector = Impl::DenseVectorView(out[i]);
          jacobian.mv(referenceCurlDenseVector, outDenseVector);
          out[i] /= integrationElement;
        }
      }
    }

  private:
    Geometry const* geometry_ = nullptr;
};

} // end namespace Dune::Functions

#endif
