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
 * transformation for H(div)-conforming finite elements.  The Jacobian
 * transformation is the affine-geometry variant used by the existing
 * GlobalValuedLocalFiniteElement implementation.
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
    void precompute(Derivatives::Partial,
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

      auto jacobianTransposed = geometry_->jacobianTransposed(x);
      auto integrationElement = geometry_->integrationElement(x);

      for (auto i : Dune::range(in.size())) {
        auto tmp = in[i];
        auto tmpDenseVector = Impl::DenseVectorView(tmp);
        auto outDenseVector = Impl::DenseVectorView(out[i]);
        jacobianTransposed.mtv(tmpDenseVector, outDenseVector);
        out[i] /= integrationElement;
      }
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

      auto jacobianTransposed = geometry_->jacobianTransposed(x);
      auto integrationElement = geometry_->integrationElement(x);

      for (auto i : Dune::range(in.size())) {
        out[i] = 0;
        for (std::size_t k=0; k<out[i].M(); ++k)
          for (std::size_t l=0; l<in[i].N(); ++l)
            for (auto&& [jacobianTransposed_l_j, j] : sparseRange(jacobianTransposed[l]))
              out[i][j][k] += jacobianTransposed_l_j * in[i][l][k];
        out[i] /= integrationElement;
      }
    }

    template<class LocalBasis, class In, class Out>
    void finalize(Derivatives::Partial d,
                  LocalBasis const& localBasis,
                  typename LocalBasis::Traits::DomainType const& x,
                  In const& in,
                  Out& out) const
    {
      assert(0 <= d.i && d.i < Geometry::coorddimension);

      std::vector<typename Traits<LocalBasis,void>::template DerivativeRange<Derivatives::Jacobian>> jacobian;
      finalize(Derivatives::Jacobian{}, localBasis, x, in, jacobian);

      out.resize(jacobian.size());
      for (auto i : Dune::range(out.size()))
        for (auto j : Dune::range(out[i].size()))
          out[i][j] = jacobian[i][j][d.i];
    }

  private:
    Geometry const* geometry_ = nullptr;
};

/**
 * \brief Transformation policy for the covariant Piola map.
 *
 * This policy preserves tangential traces and is the canonical value
 * transformation for H(curl)-conforming finite elements.  The Jacobian
 * transformation is the affine-geometry variant used by the existing
 * GlobalValuedLocalFiniteElement implementation.
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
    void precompute(Derivatives::Partial,
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
    void finalize(Derivatives::Jacobian,
                  LocalBasis const&,
                  typename LocalBasis::Traits::DomainType const& x,
                  In const& in,
                  Out& out) const
    {
      assert(!!geometry_);
      out.resize(in.size());

      auto jacobianInverseTransposed = geometry_->jacobianInverseTransposed(x);

      for (auto i : Dune::range(in.size())) {
        out[i] = 0;
        for (std::size_t j=0; j<out[i].N(); ++j)
          for (std::size_t k=0; k<out[i].M(); ++k)
            for (auto&& [jacobianInverseTransposed_j_l, l] : sparseRange(jacobianInverseTransposed[j]))
              out[i][j][k] += jacobianInverseTransposed_j_l * in[i][l][k];
      }
    }

    template<class LocalBasis, class In, class Out>
    void finalize(Derivatives::Partial d,
                  LocalBasis const& localBasis,
                  typename LocalBasis::Traits::DomainType const& x,
                  In const& in,
                  Out& out) const
    {
      assert(0 <= d.i && d.i < Geometry::coorddimension);

      std::vector<typename Traits<LocalBasis,void>::template DerivativeRange<Derivatives::Jacobian>> jacobian;
      finalize(Derivatives::Jacobian{}, localBasis, x, in, jacobian);

      out.resize(jacobian.size());
      for (auto i : Dune::range(out.size()))
        for (auto j : Dune::range(out[i].size()))
          out[i][j] = jacobian[i][j][d.i];
    }

  private:
    Geometry const* geometry_ = nullptr;
};

} // end namespace Dune::Functions

#endif
