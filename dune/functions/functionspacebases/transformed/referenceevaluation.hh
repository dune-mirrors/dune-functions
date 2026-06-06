// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_REFERENCEEVALUATION_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_REFERENCEEVALUATION_HH

#include <array>
#include <vector>

#include <dune/common/rangeutilities.hh>

#include <dune/functions/functionspacebases/transformed/derivative.hh>
#include <dune/functions/functionspacebases/transformed/derivativetraits.hh>

namespace Dune::Functions {

namespace Impl {

template<class LocalBasis, class Derivative>
struct ReferenceDerivativeRange;

template<class LocalBasis>
struct ReferenceDerivativeRange<LocalBasis,Derivatives::Value>
{
  using type = typename LocalBasis::Traits::RangeType;
};

template<class LocalBasis>
struct ReferenceDerivativeRange<LocalBasis,Derivatives::Jacobian>
{
  using type = typename LocalBasis::Traits::JacobianType;
};

template<class LocalBasis>
struct ReferenceDerivativeRange<LocalBasis,Derivatives::Gradient>
  : ReferenceDerivativeRange<LocalBasis,Derivatives::Jacobian>
{};

template<class LocalBasis>
struct ReferenceDerivativeRange<LocalBasis,Derivatives::Divergence>
  : ReferenceDerivativeRange<LocalBasis,Derivatives::Jacobian>
{};

template<class LocalBasis>
struct ReferenceDerivativeRange<LocalBasis,Derivatives::Curl>
  : ReferenceDerivativeRange<LocalBasis,Derivatives::Jacobian>
{};

template<class LocalBasis>
struct ReferenceDerivativeRange<LocalBasis,Derivatives::Hessian>
{
  using type = typename PullbackPrecomputeBuffer<LocalBasis,Derivatives::Hessian>::type::value_type;
};

template<class LocalBasis>
struct ReferenceDerivativeRange<LocalBasis,Derivatives::Laplacian>
  : ReferenceDerivativeRange<LocalBasis,Derivatives::Hessian>
{};

} // namespace Impl

/**
 * \brief Evaluate the reference quantity needed for a transformed derivative.
 */
struct ReferenceLocalBasisEvaluator
{
  template<class Derivative, class LocalBasis>
  using Range = typename Impl::ReferenceDerivativeRange<LocalBasis,Derivative>::type;

  template<class LocalBasis, class Out>
  void evaluate(Derivatives::Value,
                LocalBasis const& localBasis,
                typename LocalBasis::Traits::DomainType const& x,
                Out& out) const
  {
    localBasis.evaluateFunction(x,out);
  }

  template<class Derivative, class LocalBasis, class Out>
    requires (std::is_same_v<Derivative,Derivatives::Jacobian>
           || std::is_same_v<Derivative,Derivatives::Gradient>
           || std::is_same_v<Derivative,Derivatives::Divergence>
           || std::is_same_v<Derivative,Derivatives::Curl>)
  void evaluate(Derivative,
                LocalBasis const& localBasis,
                typename LocalBasis::Traits::DomainType const& x,
                Out& out) const
  {
    localBasis.evaluateJacobian(x,out);
  }

  template<class Derivative, class LocalBasis, class Out>
    requires (std::is_same_v<Derivative,Derivatives::Hessian>
           || std::is_same_v<Derivative,Derivatives::Laplacian>)
  void evaluate(Derivative,
                LocalBasis const& localBasis,
                typename LocalBasis::Traits::DomainType const& x,
                Out& out) const
  {
    if constexpr (requires { localBasis.evaluateHessian(x,out); })
      localBasis.evaluateHessian(x,out);
    else {
      out.resize(localBasis.size());
      std::vector<typename LocalBasis::Traits::RangeType> partialValues;
      std::array<unsigned int,LocalBasis::Traits::dimDomain> order;

      for (auto i : Dune::range(LocalBasis::Traits::dimDomain))
        for (auto j : Dune::range(i,LocalBasis::Traits::dimDomain)) {
          order.fill(0);
          ++order[i];
          ++order[j];
          localBasis.partial(order,x,partialValues);

          for (auto k : Dune::range(out.size()))
            for (auto r : Dune::range(LocalBasis::Traits::dimRange)) {
              out[k][r][i][j] = partialValues[k][r];
              out[k][r][j][i] = partialValues[k][r];
            }
        }
    }
  }
};

} // namespace Dune::Functions

#endif
