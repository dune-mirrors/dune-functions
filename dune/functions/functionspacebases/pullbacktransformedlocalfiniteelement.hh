// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PULLBACKTRANSFORMEDLOCALFINITEELEMENT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PULLBACKTRANSFORMEDLOCALFINITEELEMENT_HH

#include <array>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

namespace Dune::Functions {

namespace Derivatives
{
  template <class D>
  D local (D const& d) { return d; }

  struct Value {};
  struct Jacobian {};
  struct Gradient {};
  struct Divergence {};
  struct Hessian {};
  struct Laplacian {};
  struct Partial { int i; };

} // end namespace Derivatives

template <class LocalBasis, class Geometry>
struct StandardDerivativeTraits
{
  template <class Derivative>
  struct Traits;

  template <>
  struct Traits<Derivatives::Value> {
    using type = typename LocalBasis::Traits::RangeType;
  };

  template <>
  struct Traits<Derivatives::Jacobian> {
    using K = typename LocalBasis::Traits::RangeFieldType;
    static constexpr int dimRange = typename LocalBasis::Traits::dimRange;
    using type = FieldMatrix<K,dimRange,Geometry::coordinatedimension>;
  };

  // template <>
  // struct Traits<Derivatives::Hessian> {
  //   using type = typename LocalBasis::Traits::HessianType;
  // };
};

template <class LocalBasis, class Geometry>
struct ScalarDerivativeTraits
{
  template <class Derivative>
  struct Traits;

  template <>
  struct Traits<Derivatives::Value> {
    using type = typename LocalBasis::Traits::RangeType::value_type;
  };

  template <>
  struct Traits<Derivatives::Jacobian> {
    using K = typename LocalBasis::Traits::RangeFieldType;
    using type = FieldVector<K,Geometry::coordinatedimension>;
  };

  // template <>
  // struct Traits<Derivatives::Hessian> {
  //   using type = typename LocalBasis::Traits::HessianType;
  // };
};

/**
 * \brief Implementation of a dune-localfunctions LocalBasis that applies a
 * linear basis transformation
 *
 * \tparam FEImplementation The finite element implementation
 * \tparam LocalBasisTraits Traits of the local basis
 * \tparam LocalToGlobalTraits Transformation of the LocalBasisTraits
 */
template <class LocalBasis, class Geometry,
          class DerivativeTraits = StandardDerivativeTraits<LocalBasis,Geometry>>
class PullbackTransformedLocalBasis
{
  private:

    template <class D>
    struct PrecomputeBufferTraits;

    template <>
    struct PrecomputeBufferTraits<Derivatives::Value> {
      using type = std::vector<typename LocalBasis::Traits::RangeType>;
    };

    template <>
    struct PrecomputeBufferTraits<Derivatives::Jacobian> {
      using type = std::vector<typename LocalBasis::Traits::JacobianType>;
    };

  public:

    template <class D>
    using PrecomputeBuffer = typename PrecomputeBufferTraits<D>::type;

    template <class D>
    using DerivativeRange = typename DerivativeTraits::template Traits<D>::type;

    using Domain = typename LocalBasis::Traits::DomainType;
    using Range = DerivativeRange<Derivatives::Value>;

  public:
    /**
     * \brief Number of shape functions
     * This need not to be equal to the size of the reference local basis
     */
    auto size() const
    {
      return localBasis_->size();
    }

    void bind (LocalBasis const& localBasis)
    {
      localBasis_ = &localBasis;
    }

    void bind (Geometry const& geometry)
    {
      geometry_ = &geometry;
    }

    void precompute(Derivatives::Value,
                    Domain const& x,
                    PrecomputeBuffer<Derivatives::Value>& out) const
    {
      assert(!!localBasis_);
      localBasis_->evaluateFunction(x,out);
    }

    void precompute(Derivatives::Jacobian,
                    Domain const& x,
                    PrecomputeBuffer<Derivatives::Jacobian>& out) const
    {
      assert(!!localBasis_);
      localBasis_->evaluateJacobian(x,out);
    }

    void finalize(Derivatives::Value d,
                  Domain const& x,
                  PrecomputeBuffer<Derivatives::Value> const& in,
                  std::vector<DerivativeRange<Derivatives::Value>>& out) const
    {
      out.resize(in.size());
      using T = DerivativeRange<Derivatives::Value>;
      for (auto i : Dune::range(in.size()))
        out[i] = T(in[i]);
    }

    // auto const& finalize(Derivatives::Value d,
    //                      Domain const& x,
    //                      PrecomputeBuffer<Derivatives::Value> const& in) const
    // {
    //   return in;
    // }

    void finalize(Derivatives::Jacobian d,
                  Domain const& x,
                  PrecomputeBuffer<Derivatives::Jacobian> const& in,
                  std::vector<DerivativeRange<Derivatives::Jacobian>>& out) const
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

    // auto finalize(Derivatives::Jacobian d,
    //               Domain const& x,
    //               PrecomputeBuffer<Derivatives::Jacobian> const& in) const
    // {
    //   assert(!!geometry_);
    //   return transformedRangeView(in,
    //     [Jinv=geometry_->jacobianInverse(x)](auto const& jac) { return jac * Jinv; });
    // }

    void evaluate(Derivatives::Value d,
                  Domain const& x,
                  std::vector<Range>& out) const
    {
      precompute(d,x,out);
    }

    void evaluate(Derivatives::Jacobian d,
                  Domain const& x,
                  std::vector<DerivativeRange<Derivatives::Jacobian>>& out) const
    {
      precompute(d,x,jacobianBuffer_);
      finalize(d,x,jacobianBuffer_,out);
    }

  private:
    LocalBasis const* localBasis_ = nullptr;
    Geometry const* geometry_ = nullptr;

    mutable PrecomputeBuffer<Derivatives::Jacobian> jacobianBuffer_ = {};
};


template <class LocalFiniteElement, class Geometry,
          class DerivativeTraits = StandardDerivativeTraits<typename LocalFiniteElement::Traits::LocalBasisType,Geometry>>
class PullbackTransformedLocalFiniteElement
{
  using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;

  public:
    void bind (LocalFiniteElement const& lfe)
    {
      globalizedBasis_.bind(lfe.localBasis());
    }

    void bind (Geometry const& geometry)
    {
      globalizedBasis_.bind(geometry);
    }

    auto const& globalizedBasis() const
    {
      return globalizedBasis_;
    }

  protected:
    PullbackTransformedLocalBasis<LocalBasis, Geometry, DerivativeTraits> globalizedBasis_;
};

} // end namespace Dune::Functions

#endif
