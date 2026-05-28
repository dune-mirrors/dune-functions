// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PULLBACKTRANSFORMEDLOCALFINITEELEMENT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PULLBACKTRANSFORMEDLOCALFINITEELEMENT_HH

#include <array>
#include <cassert>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

namespace Dune::Functions {

namespace Derivatives
{
  struct Value {};
  struct Jacobian {};
  struct Gradient {};
  struct Divergence {};
  struct Hessian {};
  struct Laplacian {};
  struct Partial { int i; };

} // end namespace Derivatives

namespace Impl {

template <class LocalBasis, class Geometry, class Derivative>
struct StandardDerivativeRange;

template <class LocalBasis, class Geometry>
struct StandardDerivativeRange<LocalBasis,Geometry,Derivatives::Value>
{
  using type = typename LocalBasis::Traits::RangeType;
};

template <class LocalBasis, class Geometry>
struct StandardDerivativeRange<LocalBasis,Geometry,Derivatives::Jacobian>
{
  using K = typename LocalBasis::Traits::RangeFieldType;
  using type = FieldMatrix<K,LocalBasis::Traits::dimRange,Geometry::coorddimension>;
};

template <class LocalBasis, class Geometry, class Derivative>
struct ScalarDerivativeRange;

template <class LocalBasis, class Geometry>
struct ScalarDerivativeRange<LocalBasis,Geometry,Derivatives::Value>
{
  using type = typename LocalBasis::Traits::RangeType::value_type;
};

template <class LocalBasis, class Geometry>
struct ScalarDerivativeRange<LocalBasis,Geometry,Derivatives::Jacobian>
{
  using K = typename LocalBasis::Traits::RangeFieldType;
  using type = FieldVector<K,Geometry::coorddimension>;
};

template <class LocalBasis, class Derivative>
struct PullbackPrecomputeBuffer;

template <class LocalBasis>
struct PullbackPrecomputeBuffer<LocalBasis,Derivatives::Value>
{
  using type = std::vector<typename LocalBasis::Traits::RangeType>;
};

template <class LocalBasis>
struct PullbackPrecomputeBuffer<LocalBasis,Derivatives::Jacobian>
{
  using type = std::vector<typename LocalBasis::Traits::JacobianType>;
};

} // end namespace Impl

template <class LocalBasis, class Geometry>
struct StandardDerivativeTraits
{
  template <class Derivative>
  using Range = Impl::StandardDerivativeRange<LocalBasis,Geometry,Derivative>;
};

template <class LocalBasis, class Geometry>
struct ScalarDerivativeTraits
{
  template <class Derivative>
  using Range = Impl::ScalarDerivativeRange<LocalBasis,Geometry,Derivative>;
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
  public:

    template <class D>
    using PrecomputeBuffer = typename Impl::PullbackPrecomputeBuffer<LocalBasis,D>::type;

    template <class D>
    using DerivativeRange = typename DerivativeTraits::template Range<D>::type;

    using Domain = typename LocalBasis::Traits::DomainType;
    using Range = DerivativeRange<Derivatives::Value>;

  public:
    /**
     * \brief Number of shape functions
     * This need not to be equal to the size of the reference local basis
     */
    std::size_t size () const
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

    void precompute (Derivatives::Value,
                     Domain const& x,
                     PrecomputeBuffer<Derivatives::Value>& out) const
    {
      assert(!!localBasis_);
      localBasis_->evaluateFunction(x,out);
    }

    void precompute (Derivatives::Jacobian,
                     Domain const& x,
                     PrecomputeBuffer<Derivatives::Jacobian>& out) const
    {
      assert(!!localBasis_);
      localBasis_->evaluateJacobian(x,out);
    }

    void finalize (Derivatives::Value d,
                   Domain const& x,
                   PrecomputeBuffer<Derivatives::Value> const& in,
                   std::vector<DerivativeRange<Derivatives::Value>>& out) const
    {
      out.resize(in.size());
      using T = DerivativeRange<Derivatives::Value>;
      for (auto i : Dune::range(in.size()))
        out[i] = T(in[i]);
    }

    void finalize (Derivatives::Jacobian d,
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

    void evaluate (Derivatives::Value d,
                   Domain const& x,
                   std::vector<Range>& out) const
    {
      if constexpr(std::is_same_v<Range,typename LocalBasis::Traits::RangeType>)
        precompute(d,x,out);
      else {
        precompute(d,x,valueBuffer_);
        finalize(d,x,valueBuffer_,out);
      }
    }

    void evaluate (Derivatives::Jacobian d,
                   Domain const& x,
                   std::vector<DerivativeRange<Derivatives::Jacobian>>& out) const
    {
      precompute(d,x,jacobianBuffer_);
      finalize(d,x,jacobianBuffer_,out);
    }

    int order () const
    {
      return localBasis_->order();
    }

  private:
    LocalBasis const* localBasis_ = nullptr;
    Geometry const* geometry_ = nullptr;

    mutable PrecomputeBuffer<Derivatives::Value> valueBuffer_ = {};
    mutable PrecomputeBuffer<Derivatives::Jacobian> jacobianBuffer_ = {};
};


template <class LocalFiniteElement, class Geometry,
          class DerivativeTraits = StandardDerivativeTraits<typename LocalFiniteElement::Traits::LocalBasisType,Geometry>>
class PullbackTransformedLocalFiniteElement
{
    using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;

  public:
    using Basis = PullbackTransformedLocalBasis<LocalBasis, Geometry, DerivativeTraits>;

    void bind (LocalFiniteElement const& lfe)
    {
      basis_.bind(lfe.localBasis());
    }

    void bind (Geometry const& geometry)
    {
      basis_.bind(geometry);
      type_ = geometry.type();
    }

    Basis const& basis () const
    {
      return basis_;
    }

    std::size_t size () const
    {
      return basis_.size();
    }

    GeometryType type () const
    {
      return type_;
    }

  protected:
    Basis basis_;
    GeometryType type_;
};

} // end namespace Dune::Functions

#endif
