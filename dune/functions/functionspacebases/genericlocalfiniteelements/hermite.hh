// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_FUNCTIONS_GENERICLOCALFINITEELEMENTS_HERMITE_HH
#define DUNE_FUNCTIONS_GENERICLOCALFINITEELEMENTS_HERMITE_HH
#include <algorithm>
#include <array>
#include <numeric>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune::Functions
{
template<class DF, int n, class D, class RF, int m, class R, class J, class H>
struct H2LocalBasisTraits : public LocalBasisTraits<DF, n, D, RF, m, R, J> {
    /** \brief Type to represent the Hessian

        When \f$ \hat\phi : \mbox{IR}^n \to \mbox{IR}\f$ then HessianType
        is an 2D-array of m x m components where entry H[i][j] contains
        the derivative  \f$\partial_i \partial_j \hat\phi \f$.
      */
    using HessianType = H;
};

namespace Impl
{

// *****************************************************************************
// * Some helper functions for building polynomial bases from monomials
// *****************************************************************************


// Evaluation of 1d monomial values
template<class K>
static constexpr auto evaluateMonomialValues(const Dune::FieldVector<K,1>& x)
{
  using Range = Dune::FieldVector<K,1>;
  constexpr std::size_t maxOrder=3;
  constexpr std::size_t size = (maxOrder+1);
  auto xPowers = std::array<double,maxOrder+1>{};
  xPowers[0] = 1.0;
  for(auto k: Dune::range(maxOrder))
    xPowers[k+1] = xPowers[k]*x[0];
  auto y = Dune::FieldVector<Range,size>{};
  for(auto order : Dune::range(maxOrder+1))
    y[order] = xPowers[order];
  return y;
}

// Evaluation of 1d monomial jacobians
template<class K>
static constexpr auto evaluateMonomialJacobians(const Dune::FieldVector<K,1>& x)
{
  using Jacobian = Dune::FieldMatrix<K,1,1>;
  constexpr std::size_t maxOrder=3;
  constexpr std::size_t size = (maxOrder+1);
  auto xPowers = std::array<double,maxOrder+1>{};
  xPowers[0] = 1.0;
  for(auto k: Dune::range(maxOrder))
    xPowers[k+1] = xPowers[k]*x[0];
  auto y = Dune::FieldVector<Jacobian,size>{};
  for(auto order : Dune::range(std::size_t(1), maxOrder+1))
    y[order][0][0] = order*xPowers[order-1];
  return y;
}

// Evaluation of 1d monomial hessians
template<class K>
static constexpr auto evaluateMonomialHessians(const Dune::FieldVector<K,1>& x)
{
  using Hessian = Dune::FieldMatrix<K,1,1>;
  constexpr std::size_t maxOrder=3;
  constexpr std::size_t size = (maxOrder+1);
  auto xPowers = std::array<double,maxOrder+1>{};
  xPowers[0] = 1.0;
  for(auto k: Dune::range(maxOrder))
    xPowers[k+1] = xPowers[k]*x[0];
  auto y = Dune::FieldVector<Hessian,size>{};
  for(auto order : Dune::range(std::size_t(2), maxOrder+1))
    y[order][0][0] = order * (order-1)*xPowers[order-2];
  return y;
}

// Evaluation of 2d monomial values
template<class K>
static constexpr auto evaluateMonomialValues(const Dune::FieldVector<K,2>& x)
{
  using Range = Dune::FieldVector<K,1>;
  constexpr std::size_t maxOrder=3;
  constexpr std::size_t dim=2;
  constexpr std::size_t size = (maxOrder+1)*(maxOrder+2)/2;
  auto xPowers = std::array<std::array<double,maxOrder+1>,dim>{};
  for(auto j: Dune::range(dim))
  {
    xPowers[j][0] = 1.0;
    for(auto k: Dune::range(maxOrder))
      xPowers[j][k+1] = xPowers[j][k]*x[j];
  }
  auto y = Dune::FieldVector<Range,size>{};
  std::size_t index=0;
  for(auto order : Dune::range(maxOrder+1))
  {
    for(auto k : Dune::range(order+1))
    {
      y[index] = xPowers[0][order-k]*xPowers[1][k];
      ++index;
    }
  }
  return y;
}

// Evaluation of 2d monomial jacobians
template<class K>
static constexpr auto evaluateMonomialJacobians(const Dune::FieldVector<K,2>& x)
{
  using Jacobian = Dune::FieldMatrix<K,1,2>;
  constexpr std::size_t maxOrder=3;
  constexpr std::size_t dim=2;
  constexpr std::size_t size = (maxOrder+1)*(maxOrder+2)/2;
  auto xPowers = std::array<std::array<double,maxOrder+1>,dim>{};
  for(auto j: Dune::range(dim))
  {
    xPowers[j][0] = 1.0;
    for(auto k: Dune::range(maxOrder))
      xPowers[j][k+1] = xPowers[j][k]*x[j];
  }
  auto y = Dune::FieldVector<Jacobian,size>{};
  std::size_t index=0;
  for(auto order : Dune::range(maxOrder+1))
  {
    for(auto k : Dune::range(order+1))
    {
      if (order-k>0)
        y[index][0][0] = (order-k)*xPowers[0][order-k-1]*xPowers[1][k];
      if (k>0)
        y[index][0][1] = k*xPowers[0][order-k]*xPowers[1][k-1];
      ++index;
    }
  }
  return y;
}

// Evaluation of 2d monomial jacobians
template<class K>
static constexpr auto evaluateMonomialHessians(const Dune::FieldVector<K,2>& x)
{
  using Hessian = Dune::FieldMatrix<K,2,2>;
  constexpr std::size_t maxOrder=3;
  constexpr std::size_t dim=2;
  constexpr std::size_t size = (maxOrder+1)*(maxOrder+2)/2;
  auto xPowers = std::array<std::array<double,maxOrder+1>,dim>{};
  for(auto j: Dune::range(dim))
  {
    xPowers[j][0] = 1.0;
    for(auto k: Dune::range(maxOrder))
      xPowers[j][k+1] = xPowers[j][k]*x[j];
  }
  auto y = Dune::FieldVector<Hessian,size>{};
  std::size_t index=0;
  for(auto order : Dune::range(maxOrder+1))
  {
    for(auto k : Dune::range(order+1))
    {
      if (k < order - 2)
        y[index][0][0] = (order-k-1)*(order-k)*xPowers[0][order-k-2]*xPowers[1][k];
      if (k > 0 and k < order - 1){
        auto mixed = k*(order-k)*xPowers[0][order-k-2]*xPowers[1][k-1];
        y[index][0][1]= mixed;
        y[index][1][0]= mixed;
      }
      if (k > 1)
        y[index][0][0] = k*(k-1)*xPowers[0][order-k]*xPowers[1][k-2];

      ++index;
    }
  }
  return y;
}


// binomial coefficient for evaluation of polynomial basis size.
template <int n, int k>
constexpr int cnk(){
  if constexpr (n == 0 or k == 0 or n == k)
    return 1;
  else
    return cnk<n-1, k>() + cnk<n-1, k-1>();
  }


// Evaluation of 3-dimensial monomial values up to order p
template<int p = 3, class K>
static constexpr auto evaluateMonomialValues(const Dune::FieldVector<K,3>& x)
{
  using Range = Dune::FieldVector<K,1>;
  constexpr int dim = 3;
  constexpr std::size_t size = cnk<p + dim, dim>();
  auto xPowers = std::array<std::array<double, p + 1>, dim>{};
  for(auto j: Dune::range(dim))
  {
    xPowers[j][0] = 1.0;
    for(auto k: Dune::range(p))
      xPowers[j][k+1] = xPowers[j][k]*x[j];
  }

  auto y = Dune::FieldVector<Range,size>{};
  std::size_t index=0;
  for(auto order : Dune::range(p + 1))
  {
    for(auto k : Dune::range(order + 1))
    {
      for (auto l : Dune::range(order - k + 1))
      {
        y[index] = xPowers[0][order - k - l] * xPowers[1][l] * xPowers[2][k];
        ++index;
      }
    }
  }
  return y;
}

// Evaluation of 2d monomial jacobians
template<int p = 3, class K>
static constexpr auto evaluateMonomialJacobians(const Dune::FieldVector<K,3>& x)
{
  using Jacobian = Dune::FieldMatrix<K,1,3>;
  constexpr std::size_t maxOrder = p;
  constexpr std::size_t dim = 3;
  constexpr std::size_t size = cnk<p + dim, dim>();
  auto xPowers = std::array<std::array<double, maxOrder + 1>, dim>{};
  for(auto j: Dune::range(dim))
  {
    xPowers[j][0] = 1.0;
    for(auto k: Dune::range(maxOrder))
      xPowers[j][k+1] = xPowers[j][k]*x[j];
  }
  auto y = Dune::FieldVector<Jacobian,size>{};
  std::size_t index=0;
  for(auto order : Dune::range(p + 1))
  {
    for(auto k : Dune::range(order + 1))
    {
      for (auto l : Dune::range(order - k + 1))
      {
      if (order-k-l>0)
        y[index][0][0] = (order-k-l)*xPowers[0][order-k-l-1]*xPowers[1][l] * xPowers[2][k];
      if (l>0)
        y[index][0][1] = l*xPowers[0][order-k-l]*xPowers[1][l-1]*xPowers[2][k];
      if (k>0)
        y[index][0][2] = k * xPowers[0][order-k-l]*xPowers[1][l]*xPowers[2][k-1];
      ++index;
      }
    }
  }
  return y;
}


/**
 * \brief Implementation of hermite Polynomials
 * \tparam D Type to represent the field in the domain
   \tparam R Type to represent the field in the range
   \tparam dim Dimension of the domain simplex
*/
template<class D, class R, unsigned int dim, bool reduced>
class HermiteLocalBasis
{
  public:
    using Traits = H2LocalBasisTraits<D, dim, FieldVector<D, dim>, R, 1, FieldVector<R, 1>,
                                      FieldMatrix<R, 1, dim>, FieldMatrix<R, dim, dim>>;

  private:
    /**
    * @brief Get the Hermite Coefficients Matrix
    *
    * @tparam F Field type
    * @tparam dim dimesion of domain of Reference triangle
    * @return HermiteVecMatrix<F,dim> where size of the underlying matrix depends on dim
    */
    static constexpr auto getHermiteCoefficients()
    {
      static_assert(dim > 0 and dim < 4 and not(reduced and dim != 2));

      if constexpr (dim == 1) {
        return Dune::FieldMatrix<D, 4, 4>({{1, 0, -3, 2}, {0, 1, -2, 1}, {0, 0, 3, -2}, {0, 0, -1, 1}});
      } else if constexpr (dim == 2) {
        if constexpr (reduced) {
          auto w = std::array<D, 9>{1. / 3,  1. / 18, 1. / 18, 1. / 3, -1. / 9,
                                    1. / 18, 1. / 3,  1. / 18, -1. / 9};
          return Dune::FieldMatrix<D, 9, 10>({
              {1, 0, 0, -3, -13 + w[0] * 27, -3, 2, 13 - w[0] * 27, 13 - w[0] * 27, 2},
              {0, 1, 0, -2, -3 + w[1] * 27, 0, 1, 3 - w[1] * 27, 2 - w[1] * 27, 0},
              {0, 0, 1, 0, -3 + w[2] * 27, -2, 0, 2 - w[2] * 27, 3 - w[2] * 27, 1},
              {0, 0, 0, 3, -7 + w[3] * 27, 0, -2, 7 - w[3] * 27, 7 - w[3] * 27, 0},
              {0, 0, 0, -1, 2 + w[4] * 27, 0, 1, -2 - w[4] * 27, -2 - w[4] * 27, 0},
              {0, 0, 0, 0, -1 + w[5] * 27, 0, 0, 2 - w[5] * 27, 1 - w[5] * 27, 0},
              {0, 0, 0, 0, -7 + w[6] * 27, 3, 0, 7 - w[6] * 27, 7 - w[6] * 27, -2},
              {0, 0, 0, 0, -1 + w[7] * 27, 0, 0, 1 - w[7] * 27, 2 - w[7] * 27, 0},
              {0, 0, 0, 0, 2 + w[8] * 27, -1, 0, -2 - w[8] * 27, -2 - w[8] * 27, 1},
          });
        } else
          return Dune::FieldMatrix<D, 10,10>({{1, 0, 0, -3, -13, -3, 2, 13, 13, 2},
                                      {0, 1, 0, -2, -3, 0, 1, 3, 2, 0},
                                      {0, 0, 1, 0, -3, -2, 0, 2, 3, 1}, // l_2
                                      {0, 0, 0, 3, -7, 0, -2, 7, 7, 0},
                                      {0, 0, 0, -1, 2, 0, 1, -2, -2, 0},
                                      {0, 0, 0, 0, -1, 0, 0, 2, 1, 0},
                                      {0, 0, 0, 0, -7, 3, 0, 7, 7, -2}, // l_6
                                      {0, 0, 0, 0, -1, 0, 0, 1, 2, 0},
                                      {0, 0, 0, 0, 2, -1, 0, -2, -2, 1},
                                      {0, 0, 0, 0, 27, 0, 0, -27, -27, 0}}); // l_9, inner dof
      } else if constexpr (dim == 3) {
        return Dune::FieldMatrix<D, 20,20>({{1, 0,  0,  0, -3, -13, -3, -13, -13, -3, // deg 0 to 2
                                    2, 13, 13, 2, 13, 33,  13, 13,  13,  2}, // deg 3
                                    {0, 1, 0, 0,/*xx*/ -2, /*xy*/-3,/*yy*/ 0,/*xz*/ -3,/*yz*/ 0,/*zz*/ 0, 1, 3, 2, 0, 3, 4, 0, 2, 0, 0},
                                    {0, 0, 1, 0, 0, -3, -2, 0, -3, 0, 0, 2, 3, 1, 0, 4, 3, 0, 2, 0},
                                    {0, 0, 0, 1, 0, 0, 0, -3, -3, -2, 0, 0, 0, 0, 2, 4, 2, 3, 3, 1},
                                    {0,  0, 0, 0, 3, -7, 0, -7, 0, 0, // l_4
                                    -2, 7, 7, 0, 7, 7,  0, 7,  0, 0},
                                    {0, 0, 0, 0, -1, 2, 0, 2, 0, 0, 1, -2, -2, 0, -2, -2, 0, -2, 0, 0},
                                    {0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0},
                                    {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0},
                                    {0, 0, 0, 0,  0, -7, 3, 0, -7, 0, // l_8
                                    0, 7, 7, -2, 0, 7,  7, 0, 7,  0},
                                    {0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0},
                                    {0, 0, 0, 0, 0, 2, -1, 0, 2, 0, 0, -2, -2, 1, 0, -2, -2, 0, -2, 0},
                                    {0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0},
                                    {0, 0, 0, 0, 0, 0, 0, -7, -7, 3, // l_12
                                    0, 0, 0, 0, 7, 7, 7, 7,  7,  -2},
                                    {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0},
                                    {0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0},
                                    {0, 0, 0, 0, 0, 0, 0, 2, 2, -1, 0, 0, 0, 0, -2, -2, -2, -2, -2, 1},
                                    // l_16, from here on inner dofs
                                    {0, 0,   0,   0, 0, 27,  0, 0, 0, 0, // bottom
                                    0, -27, -27, 0, 0, -27, 0, 0, 0, 0},
                                    {0, 0, 0, 0, 0,   0,   0, 27,  0, 0, // front
                                    0, 0, 0, 0, -27, -27, 0, -27, 0, 0},
                                    {0, 0, 0, 0, 0, 0,   0,   0, 27,  0, // left
                                    0, 0, 0, 0, 0, -27, -27, 0, -27, 0},
                                    {0, 0, 0, 0, 0, 0,  0, 0, 0, 0, // right
                                    0, 0, 0, 0, 0, 27, 0, 0, 0, 0}});
      }
    }

    static constexpr auto referenceBasisCoefficients = getHermiteCoefficients();

  public:
    static_assert(not reduced || dim == 2, "Reduced Hermite element only implemented for 2d");
    static constexpr unsigned int coeffSize = (dim == 1)   ? 4
                                              : (dim == 2) ? ((reduced) ? 9 : 10)
                                                           : 20;
    HermiteLocalBasis()
    {
      if (not (dim <= 3))
        DUNE_THROW(Dune::NotImplemented, "only implemented for dim <= 3");
    }

    static constexpr unsigned int size() { return coeffSize; }

    unsigned int order() const { return 3; }

    void evaluateFunction(const typename Traits::DomainType &in,
                          std::vector<typename Traits::RangeType> &out) const
    {
      out.resize(size());
      auto monomialValues = evaluateMonomialValues(in);
      referenceBasisCoefficients.mv(monomialValues, out);
    }

    /** \brief Evaluate Jacobians of all shape functions at a given point
     *
     * \param[in]  in  The evaluation point
     * \param[out] out Jacobians of all shape functions at that point
     */
    void evaluateJacobian(const typename Traits::DomainType &in,
                          std::vector<typename Traits::JacobianType> &out) const
    {
      out.resize(size());
      auto monomialValues = evaluateMonomialJacobians(in);
      referenceBasisCoefficients.mv(monomialValues, out);
    }

    /** \brief Evaluate Hessians of all shape functions at a given point
     *
     * \param[in]  in  The evaluation point
     * \param[out] out Hessians of all shape functions at that point
     */
    void evaluateHessian(const typename Traits::DomainType &in,
                          std::vector<typename Traits::HessianType> &out) const
    {
      out.resize(size());
      auto monomialValues = evaluateMonomialHessians(in);
      referenceBasisCoefficients.mv(monomialValues, out);
    }

    /** \brief Evaluate partial derivatives of all shape functions at a given point
     *
     * \param[in] order The partial derivative to be computed, as a multi-index
     * \param[in] in  The evaluation point
     * \param[out] out Jacobians of all shape functions at that point
     */
    void partial(std::array<unsigned int, dim> const &order, const typename Traits::DomainType &in,
                 std::vector<typename Traits::RangeType> &out) const
    {
      out.resize(size());
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0)
        evaluateFunction(in, out);
      else if (totalOrder == 1){
        thread_local std::vector<typename Traits::JacobianType> jacobians;
        evaluateJacobian(in,jacobians);
        std::size_t which = std::max_element(order.begin(), order.end()) - order.begin();
        for (auto i : Dune::range(size()))
          out[i] = jacobians[i][0][which];
      }
      else
        DUNE_THROW(RangeError, "partial() not implemented for given order");
    }
};

/** \brief Associations of the Hermite degrees of freedom to subentities of the
 * reference simplex
 *
 * \tparam dim Dimension of the reference simplex
 */
template<unsigned int dim, bool reduced>
class HermiteLocalCoefficients
{
  public:
    using size_type = std::size_t;

  private:
    static constexpr size_type innerDofCodim = (dim == 2) ? 0 : 1;

  public:
    HermiteLocalCoefficients() : localKeys_(size())
    {
      static_assert(dim <= 3, "HermiteLocalCoefficients only implemented for dim<=3!");
      size_type numberOfVertices = dim + 1;
      size_type numberOfInnerDofs = (dim - 1) * (dim - 1); // probably incorrect for dim > 3
      for (size_type i = 0; i < numberOfVertices; ++i)     // subentities: vertices
      {
        for (size_type k = 0; k < numberOfVertices; ++k) // dim + 1 dofs per subentity
          localKeys_[numberOfVertices * i + k] = LocalKey(i, dim, k);
      }
      if constexpr (not reduced)
        for (size_type i = 0; i < numberOfInnerDofs; ++i) // subentities: element
          localKeys_[numberOfVertices * numberOfVertices + i] =
              LocalKey(i, innerDofCodim, 0); // inner dofs
    }

    //! number of coefficients
    static constexpr size_type size()
    {
      return dim == 1 ? 4 : dim == 2 ? ((reduced) ? 9 : 10) : 20;
    }

    //! get i'th index
    LocalKey const &localKey(size_type i) const { return localKeys_[i]; }

  private:
    std::vector<LocalKey> localKeys_;
};

/** \brief
 * Note: The wrapper classes in dune/functions/functionspacebases do not use this class.
 *
 * Evaluate the degrees of freedom of a Hermite basis. This class provides a template
 * hook pattern, which allows switching between the local(reference) Interpolation and the
 * global(physical) Interpolation if the interpolation function f provides a free derivative()
 * function. If so, we have a local interpolation if derivative(f) returns the local derivative
 * and global interpolation if derivative(f) returns the global derivative.
 *
 *  \tparam LocalBasis The corresponding set of shape functions
 */
template<int dim, bool reduced>
class HermiteLocalInterpolation
{
    using size_type = std::size_t;
    static constexpr size_type size = dim == 1 ? 4 : dim == 2 ? ((reduced) ? 9 : 10) : 20;
    static constexpr size_type numberOfVertices = dim + 1;
    static constexpr size_type innerDofCodim = (dim == 2) ? 0 : 1; // probably wrong for dim > 3
    static constexpr size_type numberOfInnerDofs =
        (dim - 1) * (dim - 1); // probably wrong for dim > 3

  public:
    /** \brief Evaluate a given function at the Lagrange nodes
     *
     * \tparam F Type of function to evaluate
     * \tparam C Type used for the values of the function
     * \param[in] ff Function to evaluate
     * \param[out] out Array of function values
     */
    template<typename F, typename C>
    void interpolate(const F &f, std::vector<C> &out) const
    {
      static_assert(dim <= 3, "HermiteLocalInterpolation only implemented for dim<=3!");

      auto &&df = derivative(f);
      out.resize(size);

      auto refElement = Dune::ReferenceElements<double, dim>::simplex();
      for (size_type i = 0; i < numberOfVertices; ++i) {
        auto const &x = refElement.position(i, dim);
        out[numberOfVertices * i] = f(x);
        auto const &grad = df(x);

        for (size_type d = 0; d < dim; ++d)
          out[numberOfVertices * i + d + 1] = getPartialDerivative(grad, d);
      }
      if constexpr (not reduced)
        for (size_type i = 0; i < numberOfInnerDofs; ++i) {
          out[numberOfVertices * numberOfVertices + i] = f(refElement.position(i, innerDofCodim));
        }
    }

  protected:
    template<class DerivativeType, class FieldType>
    FieldType getPartialDerivative(DerivativeType const &df, std::size_t i) const
    {
      DUNE_THROW(Dune::NotImplemented, "Derivative Type is neither FieldMatrix<double,1,d> nor "
                                       "FieldVector<double,d>");
    }

    template<class FieldType, int d>
    FieldType getPartialDerivative(Dune::FieldVector<FieldType, d> const &df, std::size_t i) const
    {
      return df[i];
    }

    template<class FieldType, int d>
    FieldType getPartialDerivative(Dune::FieldMatrix<FieldType, 1, d> const &df,
                                   std::size_t i) const
    {
      if (df.N() == 1)
        return df[0][i];
      else if (df.M() == 1)
        return df[i][0];
      else
        DUNE_THROW(Dune::NotImplemented, "Derivative of scalar function is a matrix!");
    }
};

} // namespace Impl

/** \brief Hermite finite element for simplices, as defined on the reference Element.
 * Note, that this is a non affine-equivalent finite element, that requires an additional transformation to the relate reference basis with the pullbacks of global basis.
 * For more Details, see <dune/functions/functionspacebases/hermitebasis.hh>.
 *
 * \tparam D Type used for domain coordinates
 * \tparam R Type used for function values
 * \tparam dim dimension of the reference element
 */
template<class D, class R, unsigned int dim, bool reduced = false>
class HermiteLocalFiniteElement
{

  public:
    /** \brief Export number types, dimensions, etc.
     */

    using Traits = LocalFiniteElementTraits<
        Impl::HermiteLocalBasis<D, R, dim, reduced>, Impl::HermiteLocalCoefficients<dim, reduced>,
        Impl::HermiteLocalInterpolation<dim, reduced>>;

    /** \brief Returns the local basis, i.e., the set of shape functions
     */
    const typename Traits::LocalBasisType &localBasis() const { return basis_; }

    /** \brief Returns the assignment of the degrees of freedom to the element
     * subentities
     */
    const typename Traits::LocalCoefficientsType &localCoefficients() const
    {
      return coefficients_;
    }

    /** \brief Returns object that evaluates degrees of freedom
     */
    const typename Traits::LocalInterpolationType &localInterpolation() const
    {
      return interpolation_;
    }

    /** \brief The number of shape functions */
    static constexpr std::size_t size() { return dim == 1 ? 4 : dim == 2 ? 10 : 20; }

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type() { return GeometryTypes::simplex(dim); }

  private:
    typename Traits::LocalBasisType basis_;
    typename Traits::LocalCoefficientsType coefficients_;
    typename Traits::LocalInterpolationType interpolation_;
};

} // namespace Dune

#endif // DUNE_C1ELEMENTS_HERMITE_HH
