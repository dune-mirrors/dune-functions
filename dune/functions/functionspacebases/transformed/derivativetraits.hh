// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_DERIVATIVETRAITS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_DERIVATIVETRAITS_HH

#include <array>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/functions/functionspacebases/transformed/derivative.hh>

namespace Dune::Functions::Impl {

template <class LocalBasis, int domainDimension>
using HessianMatrix = FieldMatrix<typename LocalBasis::Traits::RangeFieldType,domainDimension,domainDimension>;

template <class LocalBasis, int domainDimension>
using HessianTensor = std::array<HessianMatrix<LocalBasis,domainDimension>,LocalBasis::Traits::dimRange>;

template <class LocalBasis, class Geometry, int dimension>
struct CurlDerivativeRange;

template <class LocalBasis, class Geometry>
struct CurlDerivativeRange<LocalBasis,Geometry,2>
{
  using type = typename LocalBasis::Traits::RangeFieldType;
};

template <class LocalBasis, class Geometry>
struct CurlDerivativeRange<LocalBasis,Geometry,3>
{
  using type = FieldVector<typename LocalBasis::Traits::RangeFieldType,Geometry::coorddimension>;
};

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

template <class LocalBasis, class Geometry>
struct StandardDerivativeRange<LocalBasis,Geometry,Derivatives::Partial>
{
  using type = typename LocalBasis::Traits::RangeType;
};

template <class LocalBasis, class Geometry>
struct StandardDerivativeRange<LocalBasis,Geometry,Derivatives::Divergence>
{
  using type = typename LocalBasis::Traits::RangeFieldType;
};

template <class LocalBasis, class Geometry>
struct StandardDerivativeRange<LocalBasis,Geometry,Derivatives::Curl>
{
  using type = typename CurlDerivativeRange<LocalBasis,Geometry,Geometry::coorddimension>::type;
};

template <class LocalBasis, class Geometry>
struct StandardDerivativeRange<LocalBasis,Geometry,Derivatives::Hessian>
{
  using type = HessianTensor<LocalBasis,Geometry::coorddimension>;
};

template <class LocalBasis, class Geometry>
struct StandardDerivativeRange<LocalBasis,Geometry,Derivatives::Laplacian>
{
  using type = typename LocalBasis::Traits::RangeType;
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

template <class LocalBasis, class Geometry>
struct ScalarDerivativeRange<LocalBasis,Geometry,Derivatives::Partial>
{
  using type = typename LocalBasis::Traits::RangeType::value_type;
};

template <class LocalBasis, class Geometry>
struct ScalarDerivativeRange<LocalBasis,Geometry,Derivatives::Hessian>
{
  using type = HessianMatrix<LocalBasis,Geometry::coorddimension>;
};

template <class LocalBasis, class Geometry>
struct ScalarDerivativeRange<LocalBasis,Geometry,Derivatives::Laplacian>
{
  using type = typename LocalBasis::Traits::RangeType::value_type;
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

template <class LocalBasis>
struct PullbackPrecomputeBuffer<LocalBasis,Derivatives::Partial>
{
  using type = std::vector<typename LocalBasis::Traits::JacobianType>;
};

template <class LocalBasis>
struct PullbackPrecomputeBuffer<LocalBasis,Derivatives::Divergence>
{
  using type = std::vector<typename LocalBasis::Traits::JacobianType>;
};

template <class LocalBasis>
struct PullbackPrecomputeBuffer<LocalBasis,Derivatives::Curl>
{
  using type = std::vector<typename LocalBasis::Traits::JacobianType>;
};

template <class LocalBasis>
  requires requires { typename LocalBasis::Traits::HessianType; }
struct PullbackPrecomputeBuffer<LocalBasis,Derivatives::Hessian>
{
  using type = std::vector<typename LocalBasis::Traits::HessianType>;
};

template <class LocalBasis>
  requires (!requires { typename LocalBasis::Traits::HessianType; })
struct PullbackPrecomputeBuffer<LocalBasis,Derivatives::Hessian>
{
  using type = std::vector<HessianTensor<LocalBasis,LocalBasis::Traits::dimDomain>>;
};

template <class LocalBasis>
struct PullbackPrecomputeBuffer<LocalBasis,Derivatives::Laplacian>
{
  using type = typename PullbackPrecomputeBuffer<LocalBasis,Derivatives::Hessian>::type;
};

} // end namespace Dune::Functions::Impl

namespace Dune::Functions {

/**
 * \brief Derivative range traits for vector-valued pullback transformed bases.
 *
 * This traits class maps a derivative tag to the range type produced by
 * transformed local bases after the geometry pullback has been applied.
 * The nested alias template Range<Derivative>::type is the customization point
 * used by PullbackTransformedLocalBasis and draft transformation adapters.
 *
 * \tparam LocalBasis Reference local basis type.
 * \tparam Geometry Bound element geometry type.
 */
template <class LocalBasis, class Geometry>
struct StandardDerivativeTraits
{
  //! Range traits wrapper for compatibility with the initial pullback interface.
  template <class Derivative>
  using Range = Impl::StandardDerivativeRange<LocalBasis,Geometry,Derivative>;

  //! Buffer type used by reference-element precomputation.
  template <class Derivative>
  using PrecomputeBuffer = typename Impl::PullbackPrecomputeBuffer<LocalBasis,Derivative>::type;

  //! Public output range type for the quantity selected by Derivative.
  template <class Derivative>
  using DerivativeRange = typename Range<Derivative>::type;
};

/**
 * \brief Derivative range traits for scalar components of transformed bases.
 *
 * This traits class is useful when a vector-valued reference basis is used as a
 * source for scalar-valued component shape functions.  Values and partial
 * derivatives are represented by the scalar value_type of the reference range,
 * Jacobians are represented by physical-coordinate gradients, Hessians are
 * represented by matrices in physical coordinates, and Laplacians are
 * represented by scalar values.
 *
 * \tparam LocalBasis Reference local basis type.
 * \tparam Geometry Bound element geometry type.
 */
template <class LocalBasis, class Geometry>
struct ScalarDerivativeTraits
{
  //! Range traits wrapper for compatibility with the initial pullback interface.
  template <class Derivative>
  using Range = Impl::ScalarDerivativeRange<LocalBasis,Geometry,Derivative>;

  //! Buffer type used by reference-element precomputation.
  template <class Derivative>
  using PrecomputeBuffer = typename Impl::PullbackPrecomputeBuffer<LocalBasis,Derivative>::type;

  //! Public output range type for the quantity selected by Derivative.
  template <class Derivative>
  using DerivativeRange = typename Range<Derivative>::type;
};

} // end namespace Dune::Functions

#endif
