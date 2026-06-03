// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_PULLBACK_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_PULLBACK_HH

#include <array>
#include <cassert>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/functions/functionspacebases/transformed/derivative.hh>
#include <dune/functions/functionspacebases/transformed/derivativetraits.hh>

namespace Dune::Functions {

/**
 * \brief Transformation policy for affine geometry pullbacks of derivatives.
 *
 * This policy implements the chain-rule transformations used by
 * PullbackTransformedLocalBasis.  Values are converted to the configured
 * public value range, Jacobians and partial derivatives are transformed with
 * the inverse element Jacobian, and Hessians are transformed as
 * \f$J^{-T} H J^{-1}\f$.
 *
 * The \p DerivativeTraits template determines the public range types.  Use
 * StandardDerivativePullback for vector-valued bases and
 * ScalarDerivativePullback when a scalar component range is desired.
 */
struct StandardDerivativeTraitsFactory
{
  template<class LocalBasis, class Geometry>
  using Traits = StandardDerivativeTraits<LocalBasis,Geometry>;
};

struct ScalarDerivativeTraitsFactory
{
  template<class LocalBasis, class Geometry>
  using Traits = ScalarDerivativeTraits<LocalBasis,Geometry>;
};

template<class ConcreteDerivativeTraits>
struct FixedDerivativeTraitsFactory
{
  template<class LocalBasis, class Geometry>
  using Traits = ConcreteDerivativeTraits;
};

template<class Geometry,
         class DerivativeTraitsFactory = StandardDerivativeTraitsFactory>
class DerivativePullback
{
  public:
    template<class LocalBasis, class Context>
    using Traits = typename DerivativeTraitsFactory::template Traits<LocalBasis,Geometry>;

    template<class Context>
    void bind(Context const& context)
    {
      geometry_ = &context.geometry();
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
      if constexpr (requires { localBasis.evaluateHessian(x,out); }) {
        localBasis.evaluateHessian(x,out);
      }
      else {
        out.resize(localBasis.size());

        std::vector<typename LocalBasis::Traits::RangeType> partialValues;
        std::array<unsigned int,LocalBasis::Traits::dimDomain> order;

        for (auto i : Dune::range(LocalBasis::Traits::dimDomain)) {
          for (auto j : Dune::range(i,LocalBasis::Traits::dimDomain)) {
            order.fill(0);
            ++order[i];
            ++order[j];
            localBasis.partial(order,x,partialValues);

            for (auto k : Dune::range(out.size())) {
              for (auto r : Dune::range(LocalBasis::Traits::dimRange)) {
                out[k][r][i][j] = partialValues[k][r];
                out[k][r][j][i] = partialValues[k][r];
              }
            }
          }
        }
      }
    }

    template<class LocalBasis, class Out>
    void precompute(Derivatives::Laplacian,
                    LocalBasis const& localBasis,
                    typename LocalBasis::Traits::DomainType const& x,
                    Out& out) const
    {
      precompute(Derivatives::Hessian{}, localBasis, x, out);
    }

    template<class LocalBasis, class In, class Out>
    void finalize(Derivatives::Value,
                  LocalBasis const&,
                  typename LocalBasis::Traits::DomainType const&,
                  In const& in,
                  Out& out) const
    {
      out.resize(in.size());
      using T = typename Traits<LocalBasis,void>::template DerivativeRange<Derivatives::Value>;
      for (auto i : Dune::range(in.size()))
        out[i] = T(in[i]);
    }

    template<class LocalBasis, class In, class Out>
    void finalize(Derivatives::Jacobian,
                  LocalBasis const&,
                  typename LocalBasis::Traits::DomainType const& x,
                  In const& in,
                  Out& out) const
    {
      assert(!!geometry_);
      out.resize(in.size());

      if constexpr (requires{out[0][0][0];}) {
        auto&& Jinv = geometry_->jacobianInverse(x);
        for (auto i : Dune::range(in.size()))
          out[i] = in[i] * Jinv;
      }
      else {
        auto&& JinvT = geometry_->jacobianInverseTransposed(x);
        for (auto i : Dune::range(in.size()))
          JinvT.mv(in[i][0], out[i]);
      }
    }

    template<class LocalBasis, class In, class Out>
    void finalize(Derivatives::Hessian,
                  LocalBasis const&,
                  typename LocalBasis::Traits::DomainType const& x,
                  In const& in,
                  Out& out) const
    {
      assert(!!geometry_);
      out.resize(in.size());

      auto&& Jinv = geometry_->jacobianInverse(x);
      auto&& JinvT = Jinv.transposed();

      if constexpr (requires{out[0][0][0][0];}) {
        if constexpr (requires{in[0][0][0][0];}) {
          for (auto i : Dune::range(in.size()))
            for (auto r : Dune::range(LocalBasis::Traits::dimRange))
              out[i][r] = JinvT * in[i][r] * Jinv;
        }
        else {
          static_assert(LocalBasis::Traits::dimRange == 1,
            "A scalar HessianType can only be used with scalar local basis ranges.");
          for (auto i : Dune::range(in.size()))
            out[i][0] = JinvT * in[i] * Jinv;
        }
      }
      else {
        if constexpr (requires{in[0][0][0][0];}) {
          for (auto i : Dune::range(in.size()))
            out[i] = JinvT * in[i][0] * Jinv;
        }
        else {
          for (auto i : Dune::range(in.size()))
            out[i] = JinvT * in[i] * Jinv;
        }
      }
    }

    template<class LocalBasis, class In, class Out>
    void finalize(Derivatives::Laplacian,
                  LocalBasis const&,
                  typename LocalBasis::Traits::DomainType const& x,
                  In const& in,
                  Out& out) const
    {
      assert(!!geometry_);
      out.resize(in.size());

      auto&& Jinv = geometry_->jacobianInverse(x);
      auto&& JinvT = Jinv.transposed();

      if constexpr (requires{out[0][0];}) {
        for (auto i : Dune::range(in.size())) {
          out[i] = typename Traits<LocalBasis,void>::template DerivativeRange<Derivatives::Laplacian>{};
          if constexpr (requires{in[0][0][0][0];}) {
            for (auto r : Dune::range(LocalBasis::Traits::dimRange)) {
              auto hessian = JinvT * in[i][r] * Jinv;
              for (auto j : Dune::range(Geometry::coorddimension))
                out[i][r] += hessian[j][j];
            }
          }
          else {
            static_assert(LocalBasis::Traits::dimRange == 1,
              "A scalar HessianType can only be used with scalar local basis ranges.");
            auto hessian = JinvT * in[i] * Jinv;
            for (auto j : Dune::range(Geometry::coorddimension))
              out[i][0] += hessian[j][j];
          }
        }
      }
      else {
        for (auto i : Dune::range(in.size())) {
          auto hessian = [&] {
            if constexpr (requires{in[0][0][0][0];})
              return JinvT * in[i][0] * Jinv;
            else
              return JinvT * in[i] * Jinv;
          }();
          out[i] = 0;
          for (auto j : Dune::range(Geometry::coorddimension))
            out[i] += hessian[j][j];
        }
      }
    }

  private:
    Geometry const* geometry_ = nullptr;
};

template<class Geometry>
using StandardDerivativePullback = DerivativePullback<Geometry,StandardDerivativeTraitsFactory>;

template<class Geometry>
using ScalarDerivativePullback = DerivativePullback<Geometry,ScalarDerivativeTraitsFactory>;

} // end namespace Dune::Functions

#endif
