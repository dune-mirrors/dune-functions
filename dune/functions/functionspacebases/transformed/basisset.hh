// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_BASISSET_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_BASISSET_HH

#include <vector>

#include <dune/functions/functionspacebases/transformed/derivative.hh>
#include <dune/functions/functionspacebases/transformed/derivativetraits.hh>

namespace Dune::Functions {

/**
 * \brief Derivative traits for transformations acting on the complete basis set.
 *
 * The public derivative range and the precomputation buffer are both the
 * corresponding reference local-basis types.  This is suitable for
 * non-affine-equivalent scalar elements where the element transformation mixes
 * basis functions but does not change the value/Jacobian/Hessian range type.
 */
namespace Impl {

template<class LocalBasis, class Derivative>
struct BasisSetDerivativeRange;

template<class LocalBasis>
struct BasisSetDerivativeRange<LocalBasis,Derivatives::Value>
{
  using type = typename LocalBasis::Traits::RangeType;
};

template<class LocalBasis>
struct BasisSetDerivativeRange<LocalBasis,Derivatives::Jacobian>
{
  using type = typename LocalBasis::Traits::JacobianType;
};

template<class LocalBasis>
struct BasisSetDerivativeRange<LocalBasis,Derivatives::Hessian>
{
  using type = typename LocalBasis::Traits::HessianType;
};

} // namespace Impl

template<class LocalBasis, class Context>
struct BasisSetDerivativeTraits
{
  template<class Derivative>
  using DerivativeRange = typename Impl::BasisSetDerivativeRange<LocalBasis,Derivative>::type;

  template<class Derivative>
  using PrecomputeBuffer = std::vector<DerivativeRange<Derivative>>;
};

/**
 * \brief Base class for transformations that linearly mix basis-function indices.
 *
 * The derived class has to implement
 *
 * \code{.cpp}
 * template<class In, class Out>
 * void transformBasisSet(In const& in, Out& out) const;
 * \endcode
 *
 * This helper implements the common staged precomputation for values,
 * Jacobians, Hessians, and Laplacians by evaluating the corresponding reference
 * basis quantity first and then applying the same basis-set transformation to
 * the whole result vector.
 */
template<class Derived>
class BasisSetTransformationMixin
{
  public:
    template<class LocalBasis, class Context>
    using Traits = BasisSetDerivativeTraits<LocalBasis,Context>;

    template<class Context>
    void bind(Context const&)
    {}

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
      localBasis.evaluateHessian(x,out);
    }

    template<class Derivative, class LocalBasis, class In, class Out>
    void finalize(Derivative,
                  LocalBasis const&,
                  typename LocalBasis::Traits::DomainType const&,
                  In const& in,
                  Out& out) const
    {
      static_cast<Derived const&>(*this).transformBasisSet(in,out);
    }
};

} // end namespace Dune::Functions

#endif
