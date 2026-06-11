// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_REFERENCEEVALUATION_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_REFERENCEEVALUATION_HH

#include <array>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/std/no_unique_address.hh>

#include <dune/functions/functionspacebases/transformed/derivative.hh>

namespace Dune::Functions {

namespace Impl {

template <class C>
static constexpr int rank = std::rank_v<C>;

template <class K, std::size_t N>
static constexpr int rank<std::array<K,N>> = 1+rank<K>;

template <class K, int N>
static constexpr int rank<Dune::FieldVector<K,N>> = 1+rank<K>;

template <class K, int N, int M>
static constexpr int rank<Dune::FieldMatrix<K,N,M>> = 2+rank<K>;



template<class LocalBasis>
using ReferenceHessianMatrix = FieldMatrix<
  typename LocalBasis::Traits::RangeFieldType,
  LocalBasis::Traits::dimDomain,
  LocalBasis::Traits::dimDomain>;

template<class LocalBasis>
using ReferenceHessianTensor = std::array<
  ReferenceHessianMatrix<LocalBasis>,
  LocalBasis::Traits::dimRange>;

template<class LocalBasis, bool hasHessianType>
struct ReferenceHessianRange;

template<class LocalBasis>
struct ReferenceHessianRange<LocalBasis,true>
{
  using type = typename LocalBasis::Traits::HessianType;
};

template<class LocalBasis>
struct ReferenceHessianRange<LocalBasis,false>
{
  using RangeType = typename LocalBasis::Traits::RangeType;
  using type = std::conditional_t<
    rank<RangeType> == 0, ReferenceHessianMatrix<LocalBasis>, std::conditional_t<
    rank<RangeType> == 1, ReferenceHessianTensor<LocalBasis>, void>>;
};

template<class LocalBasis, class Derivative>
struct ReferenceOperatorRange;

template<class LocalBasis>
struct ReferenceOperatorRange<LocalBasis,Derivatives::Value>
{
  using type = typename LocalBasis::Traits::RangeType;
};

template<class LocalBasis>
struct ReferenceOperatorRange<LocalBasis,Derivatives::Jacobian>
{
  using type = typename LocalBasis::Traits::JacobianType;
};

template<class LocalBasis>
struct ReferenceOperatorRange<LocalBasis,Derivatives::Gradient>
{
  static_assert(LocalBasis::Traits::dimRange == 1);
  using type = FieldVector<
    typename LocalBasis::Traits::RangeFieldType,
    LocalBasis::Traits::dimDomain>;
};

template<class LocalBasis>
struct ReferenceOperatorRange<LocalBasis,Derivatives::Divergence>
{
  using type = typename LocalBasis::Traits::RangeFieldType;
};

template<class LocalBasis>
struct ReferenceOperatorRange<LocalBasis,Derivatives::DivDiv>
{
  using type = typename LocalBasis::Traits::RangeFieldType;
};

template<class LocalBasis>
struct ReferenceOperatorRange<LocalBasis,Derivatives::Curl>
{
  using Field = typename LocalBasis::Traits::RangeFieldType;
  using type = std::conditional_t<
    LocalBasis::Traits::dimDomain == 2,
    Field,
    FieldVector<Field,3>>;
};

template<class LocalBasis>
struct ReferenceOperatorRange<LocalBasis,Derivatives::Hessian>
{
  using type = typename ReferenceHessianRange<
    LocalBasis,
    requires { typename LocalBasis::Traits::HessianType; }>::type;
};

template<class LocalBasis>
struct ReferenceOperatorRange<LocalBasis,Derivatives::Laplacian>
{
  using type = typename LocalBasis::Traits::RangeType;
};

} // namespace Impl

/**
 * \brief Evaluate named differential operators on a reference local basis.
 *
 * Each derivative tag denotes the corresponding mathematical operator on the
 * reference element. Selecting a different reference operator as input for a
 * physical transformation is the responsibility of ReferenceEvaluation.
 */
struct ReferenceLocalBasisEvaluator
{
  template<class Derivative, class LocalBasis>
  using Range = typename Impl::ReferenceOperatorRange<LocalBasis,Derivative>::type;

  template<class LocalBasis, class Out>
  void evaluate(Derivatives::Value,
                LocalBasis const& localBasis,
                typename LocalBasis::Traits::DomainType const& x,
                Out& out) const
  {
    localBasis.evaluateFunction(x,out);
  }

  template<class LocalBasis, class Out>
  void evaluate(Derivatives::Jacobian,
                LocalBasis const& localBasis,
                typename LocalBasis::Traits::DomainType const& x,
                Out& out) const
  {
    localBasis.evaluateJacobian(x,out);
  }

  template<class LocalBasis, class Out>
  void evaluate(Derivatives::Gradient,
                LocalBasis const& localBasis,
                typename LocalBasis::Traits::DomainType const& x,
                Out& out) const
  {
    if constexpr (requires { localBasis.evaluateGradient(x,out); })
      localBasis.evaluateGradient(x,out);
    else {
      std::vector<Range<Derivatives::Jacobian,LocalBasis>> jacobians;
      localBasis.evaluateJacobian(x,jacobians);
      out.resize(jacobians.size());
      for (auto i : Dune::range(out.size()))
        out[i] = jacobians[i][0];
    }
  }

  template<class LocalBasis, class Out>
  void evaluate(Derivatives::Divergence,
                LocalBasis const& localBasis,
                typename LocalBasis::Traits::DomainType const& x,
                Out& out) const
  {
    if constexpr (requires { localBasis.evaluateDivergence(x,out); })
      localBasis.evaluateDivergence(x,out);
    else {
      std::vector<Range<Derivatives::Jacobian,LocalBasis>> jacobians;
      localBasis.evaluateJacobian(x,jacobians);
      out.resize(jacobians.size());
      for (auto i : Dune::range(out.size())) {
        out[i] = {};
        for (auto j : Dune::range(LocalBasis::Traits::dimDomain))
          out[i] += jacobians[i][j][j];
      }
    }
  }

  template<class LocalBasis, class Out>
  void evaluate(Derivatives::DivDiv,
                LocalBasis const& localBasis,
                typename LocalBasis::Traits::DomainType const& x,
                Out& out) const
    requires requires { localBasis.evaluateDivDiv(x,out); }
  {
    localBasis.evaluateDivDiv(x,out);
  }

  template<class LocalBasis, class Out>
    requires (LocalBasis::Traits::dimDomain == 2
           || LocalBasis::Traits::dimDomain == 3)
  void evaluate(Derivatives::Curl,
                LocalBasis const& localBasis,
                typename LocalBasis::Traits::DomainType const& x,
                Out& out) const
  {
    if constexpr (requires { localBasis.evaluateCurl(x,out); })
      localBasis.evaluateCurl(x,out);
    else {
      std::vector<Range<Derivatives::Jacobian,LocalBasis>> jacobians;
      localBasis.evaluateJacobian(x,jacobians);
      out.resize(jacobians.size());
      if constexpr (LocalBasis::Traits::dimDomain == 2)
        for (auto i : Dune::range(out.size()))
          out[i] = jacobians[i][1][0] - jacobians[i][0][1];
      else
        for (auto i : Dune::range(out.size())) {
          out[i][0] = jacobians[i][2][1] - jacobians[i][1][2];
          out[i][1] = jacobians[i][0][2] - jacobians[i][2][0];
          out[i][2] = jacobians[i][1][0] - jacobians[i][0][1];
        }
    }
  }

  template<class LocalBasis, class Out>
  void evaluate(Derivatives::Hessian,
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

          for (auto k : Dune::range(out.size())) {
            if constexpr (Impl::rank<typename Out::value_type> == 2) {
              out[k][i][j] = partialValues[k][0];
              out[k][j][i] = partialValues[k][0];
            }
            else
              for (auto r : Dune::range(LocalBasis::Traits::dimRange)) {
                out[k][r][i][j] = partialValues[k][r];
                out[k][r][j][i] = partialValues[k][r];
              }
          }
        }
    }
  }

  template<class LocalBasis, class Out>
  void evaluate(Derivatives::Laplacian,
                LocalBasis const& localBasis,
                typename LocalBasis::Traits::DomainType const& x,
                Out& out) const
  {
    if constexpr (requires { localBasis.evaluateLaplacian(x,out); })
      localBasis.evaluateLaplacian(x,out);
    else {
      std::vector<Range<Derivatives::Hessian,LocalBasis>> hessians;
      evaluate(Derivatives::Hessian{},localBasis,x,hessians);
      out.resize(hessians.size());
      for (auto i : Dune::range(out.size())) {
        out[i] = {};
        if constexpr (Impl::rank<Range<Derivatives::Hessian,LocalBasis>> == 2)
          for (auto j : Dune::range(LocalBasis::Traits::dimDomain))
            out[i][0] += hessians[i][j][j];
        else
          for (auto r : Dune::range(LocalBasis::Traits::dimRange))
            for (auto j : Dune::range(LocalBasis::Traits::dimDomain))
              out[i][r] += hessians[i][r][j][j];
      }
    }
  }
};

//! Keep the requested derivative as reference evaluation operator.
struct IdentityReferenceOperator
{
  template<class Derivative>
  using Operator = Derivative;
};

/**
 * \brief Adapt requested derivatives to the reference operators needed by a transformation.
 */
template<class OperatorPolicy = IdentityReferenceOperator,
         class Evaluator = ReferenceLocalBasisEvaluator>
class ReferenceEvaluation
{
  public:
    template<class Derivative>
    using Operator = typename OperatorPolicy::template Operator<Derivative>;

    template<class Derivative, class LocalBasis>
    using Range = typename Evaluator::template Range<Operator<Derivative>,LocalBasis>;

    ReferenceEvaluation() = default;

    explicit ReferenceEvaluation(Evaluator evaluator)
      : evaluator_(std::move(evaluator))
    {}

    template<class Derivative, class LocalBasis, class Out>
    void evaluate(Derivative,
                  LocalBasis const& localBasis,
                  typename LocalBasis::Traits::DomainType const& x,
                  Out& out) const
    {
      evaluator_.evaluate(Operator<Derivative>{},localBasis,x,out);
    }

  private:
    DUNE_NO_UNIQUE_ADDRESS Evaluator evaluator_;
};

} // namespace Dune::Functions

#endif
