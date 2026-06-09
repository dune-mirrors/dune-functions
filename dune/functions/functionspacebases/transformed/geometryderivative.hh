// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_GEOMETRYDERIVATIVE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_GEOMETRYDERIVATIVE_HH

#include <cassert>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/functions/functionspacebases/transformed/derivative.hh>
#include <dune/functions/functionspacebases/transformed/pipeline.hh>

namespace Dune::Functions {

namespace Impl {

template<class Derivative, class LocalBasis, class Geometry, class InputRange>
struct GeometryDerivativeOutputRange
{
  using type = InputRange;
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
 * \brief Pointwise transformation from reference to physical scalar derivatives.
 */
template<class Geometry>
class GeometryDerivativeTransformation
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
                   InputRange const& in,
                   OutputRange& out) const
    {
      out = in;
    }

    template<class LocalBasis, class InputRange, class OutputRange>
    void transform(Derivatives::Jacobian,
                   LocalBasis const&,
                   typename LocalBasis::Traits::DomainType const& x,
                   InputRange const& in,
                   OutputRange& out) const
    {
      assert(!!geometry_);
      out = in * geometry_->jacobianInverse(x);
    }

    template<class LocalBasis, class InputRange, class OutputRange>
    void transform(Derivatives::Gradient,
                   LocalBasis const&,
                   typename LocalBasis::Traits::DomainType const& x,
                   InputRange const& in,
                   OutputRange& out) const
    {
      static_assert(LocalBasis::Traits::dimRange == 1);
      assert(!!geometry_);
      geometry_->jacobianInverseTransposed(x).mv(in[0],out);
    }

    template<class LocalBasis, class InputRange, class OutputRange>
      requires (Geometry::mydimension == Geometry::coorddimension)
    void transform(Derivatives::Hessian,
                   LocalBasis const&,
                   typename LocalBasis::Traits::DomainType const& x,
                   InputRange const& in,
                   OutputRange& out) const
    {
      assert(!!geometry_);
      if (!geometry_->affine())
        DUNE_THROW(NotImplemented,
          "Physical Hessians are currently supported only for affine geometries");
      auto&& Jinv = geometry_->jacobianInverse(x);
      out = Jinv.transposed() * in * Jinv;
    }

    template<class LocalBasis, class InputRange, class OutputRange>
      requires (Geometry::mydimension == Geometry::coorddimension)
    void transform(Derivatives::Laplacian,
                   LocalBasis const&,
                   typename LocalBasis::Traits::DomainType const& x,
                   InputRange const& in,
                   OutputRange& out) const
    {
      assert(!!geometry_);
      if (!geometry_->affine())
        DUNE_THROW(NotImplemented,
          "Physical Laplacians are currently supported only for affine geometries");

      auto&& Jinv = geometry_->jacobianInverse(x);
      auto hessian = Jinv.transposed() * in * Jinv;
      out = {};
      for (auto j : Dune::range(Geometry::coorddimension))
        out[0] += hessian[j][j];
    }

  private:
    Geometry const* geometry_ = nullptr;
};

template<class Geometry>
using GeometryDerivativeStage =
  RangeTransformationStage<GeometryDerivativeTransformation<Geometry>>;

} // end namespace Dune::Functions

#endif
