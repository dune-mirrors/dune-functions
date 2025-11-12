// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMEDFINITEELEMENTMIXIN_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMEDFINITEELEMENTMIXIN_HH

#include <array>
#include <type_traits>
#include <utility>
#include <vector>

#include <array>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/std/type_traits.hh>

namespace Dune::Functions::Impl {
  namespace TypeDeduction{

    template<class Traits, class Dummy = int >
    struct JacobianType{
      using type = bool;
    };

    template<class Traits>
    struct JacobianType<Traits, decltype(std::declval<typename Traits::JacobianType>(), 0)>
    {
      using type = typename Traits::JacobianType;
    };


    template<class Traits, class Dummy = int >
    struct HessianType{
      using type = bool;
    };

    template<class Traits>
    struct HessianType<Traits, decltype(std::declval<typename Traits::HessianType>(), 0)>
    {
      using type = typename Traits::HessianType;
    };

    template<class Traits, class Dummy = int >
    struct DivergenceType{
      using type = bool;
    };

    template<class Traits>
    struct DivergenceType<Traits, decltype(std::declval<typename Traits::DivergenceType>(), 0)>
    {
      using type = typename Traits::DivergenceType;
    };
  }
/**
 * \brief Implementation of a dune-localfunctions LocalBasis that applies a
 * linear basis transformation
 *
 * \tparam FEImplementation The finite element implementation
 * \tparam ReferenceLocalBasisTraits LocalBasisTraits of the reference local basis
 */
template<class FEImplementation, class ReferenceLocalBasisTraits>
class TransformedLocalBasis
{
  using JT = TypeDeduction::JacobianType<ReferenceLocalBasisTraits>::type;

  using HT = TypeDeduction::HessianType<ReferenceLocalBasisTraits>::type;

  using DT = TypeDeduction::DivergenceType<ReferenceLocalBasisTraits>::type;


  public:
    using Traits = ReferenceLocalBasisTraits;

    TransformedLocalBasis(FEImplementation const& feImpl)
      : feImpl_(&feImpl)
    {}

  public:
    /**
     * \brief Number of shape functions
     * This need not to be equal to the size of the reference local basis
     */
    auto size() const
    {
      return feImpl_->size();
    }

    //! \brief Evaluate all shape functions
    void evaluateFunction(const typename Traits::DomainType &x,
                          std::vector<typename Traits::RangeType> &out) const
    {
      feImpl_->referenceLocalBasis().evaluateFunction(x, rangeBuffer_);
      out.resize(size());
      feImpl_->transform(rangeBuffer_, out);
    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     *
     * \param x Point in the reference element where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
    void evaluateJacobian(const typename Traits::DomainType &x,
                          std::vector<typename Traits::JacobianType> &out) const
    {
      feImpl_->referenceLocalBasis().evaluateJacobian(x, jacobianBuffer_);
      out.resize(size());
      feImpl_->transform(jacobianBuffer_, out);
    }

    /**
     * \brief Evaluate Hessian of all shape functions
     *
     * \note SFINA protected: this method is only available, if the wrapped LocalBasis
     * 1. exports a <HessianType> and
     * 2. provides a evaluateHessian() method with a corresponding signature
     *
     * \param x Point in the reference element where to evaluation the Hessians
     * \param[out] out The Hessians of all shape functions at the point x
     */
    template<class TT,
             std::enable_if_t<std::is_same_v<TT,HT>, int> = 0>
    void evaluateHessian(const typename Traits::DomainType &x,
                         std::vector<TT> &out) const
    {
      feImpl_->referenceLocalBasis().evaluateHessian(x, hessianBuffer_);
      out.resize(size());
      feImpl_->transform(hessianBuffer_, out);
    }

    template<class TT,
             std::enable_if_t<std::is_same_v<TT,HT>, int> = 0>
    void evaluateDivergence(const typename Traits::DomainType &x,
                         std::vector<TT> &out) const
    {
      feImpl_->referenceLocalBasis().evaluateDivergence(x, divergenceBuffer_);
      out.resize(size());
      feImpl_->transform(divergenceBuffer_, out);
    }
    /**
     * \brief Evaluate partial derivatives of any order of all shape functions
     *
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out The desired partial derivatives
     */
    void partial(std::array<unsigned int, Traits::dimDomain> const &order,
                 const typename Traits::DomainType &x,
                 std::vector<typename Traits::RangeType> &out) const
    {
      feImpl_->referenceLocalBasis().partial(order, x, rangeBuffer_);
      out.resize(size());
      feImpl_->transform(rangeBuffer_, out);
    }

    //! \brief Polynomial order of the shape functions
    auto order() const { return feImpl_->referenceLocalBasis().order(); }

  private:
    FEImplementation const* feImpl_;
    mutable std::vector<typename Traits::RangeType> rangeBuffer_;
    mutable std::vector<JT> jacobianBuffer_;
    mutable std::vector<HT> hessianBuffer_;
    mutable std::vector<DT> divergenceBuffer_;
};



/**
 * \brief A mixin for implementing a LocalFiniteElement using a
 * linear basis transformation
 *
 * \tparam FEImplementation The finite element implementation
 * \tparam ReferenceLocalBasisTraits LocalBasisTraits of the reference local basis
 *
 * The derived class should implement localCoefficients() and localInterpolation()
 * manually and provide transform() and referenceLocalBasis() as below:
 *
 * \code{.cpp}
 * template<class ReferenceLocalBasisTraits>
 * class TransformedFEExample : TransformedFiniteElementMixin<TransformedFEExample<ReferenceLocalBasisTraits>, ReferenceLocalBasisTraits>
 * {
 *     friend class TransformedLocalBasis<TransformedFEExample, ReferenceLocalBasisTraits>;
 *
 *   protected:
 *
 *     auto const& referenceLocalBasis();
 *
 *     template<class InputValues, class OutputValues>
 *     void transform(InputValues const &inValues, OutputValues &outValues) const;
 *
 *   public:
 *     auto localCoefficients();
 *     auto localInterpolation();
 *     auto size();
 *
 * };
 * \endcode
 */
template<class FEImplementation, class ReferenceLocalBasisTraits>
class TransformedFiniteElementMixin
{
  public:
    TransformedFiniteElementMixin()
      : tlb_(this->asImpl())
    {}

    TransformedFiniteElementMixin(TransformedFiniteElementMixin const& other)
    : TransformedFiniteElementMixin()
    {}

    const FEImplementation& asImpl() const
    {
      return *(static_cast<FEImplementation const*>(this));
    }

    auto const& localBasis() const
    {
      return tlb_;
    }

  protected:
    TransformedLocalBasis<FEImplementation, ReferenceLocalBasisTraits> tlb_;
};

} // end namespace Dune::Functions::Impl

#endif
