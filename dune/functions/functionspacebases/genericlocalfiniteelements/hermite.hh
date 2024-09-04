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

#include <dune/functions/analyticfunctions/monomialset.hh>
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

  template<class KCoeff, int sizePolynom, class KMonom, int sizeMonom>
  void multiplyWithCoefficentMatrix(Dune::FieldMatrix<KCoeff, sizePolynom, sizeMonom> const& coefficients,
                                   Dune::FieldVector<KMonom, sizeMonom> const& monomialValues,
                                   std::vector<Dune::FieldVector<typename Dune::PromotionTraits<KCoeff, KMonom>::PromotedType, 1>>& polynomialValues)
  {
    polynomialValues.resize(sizePolynom);
    std::fill(std::begin(polynomialValues), std::end(polynomialValues), 0.);
    for (auto&& i : Dune::range(sizePolynom))
      for (auto&& j : Dune::range(sizeMonom))
        polynomialValues[i][0] += coefficients[i][j]*monomialValues[j];
  }

  template<class KCoeff, int sizePolynom, class KMonom, int sizeMonom, int dim>
  void multiplyWithCoefficentMatrix(Dune::FieldMatrix<KCoeff, sizePolynom, sizeMonom> const& coefficients,
                                   Dune::FieldMatrix<KMonom, sizeMonom, dim> const& monomialJacobians,
                                   std::vector<Dune::FieldMatrix<typename Dune::PromotionTraits<KCoeff, KMonom>::PromotedType, 1,dim>>& polynomialJacobians)
  {
    polynomialJacobians.resize(sizePolynom);
    std::fill(std::begin(polynomialJacobians), std::end(polynomialJacobians), 0.);
    for (auto&& i : Dune::range(sizePolynom))
      for (auto&& j : Dune::range(sizeMonom))
        polynomialJacobians[i][0] += coefficients[i][j]*monomialJacobians[j];
  }

  template<class KCoeff, int sizePolynom, class KMonom, int sizeMonom, int dim>
  void multiplyWithCoefficentMatrix(Dune::FieldMatrix<KCoeff, sizePolynom, sizeMonom> const& coefficients,
                                   std::array<Dune::FieldMatrix<KMonom, dim, dim>, sizeMonom> const& monomialHessians,
                                   std::vector<Dune::FieldMatrix<typename Dune::PromotionTraits<KCoeff, KMonom>::PromotedType, dim, dim>>& polynomialHessians)
  {
    polynomialHessians.resize(sizePolynom);
    std::fill(std::begin(polynomialHessians), std::end(polynomialHessians), 0.);
    for (auto&& i : Dune::range(sizePolynom))
      for (auto&& j : Dune::range(sizeMonom))
        polynomialHessians[i] += coefficients[i][j]*monomialHessians[j];
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
      * @return FieldMatrix<F, (possibly reduced) size, size>
      *  where size is the dimension of the cubic polynomial space
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
      MonomialSet<typename Traits::RangeFieldType, dim, 3> monomials;

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
        auto monomialValues = monomials(in);
        // thread_local auto monomialValues = evaluateMonomialValues(in);
        multiplyWithCoefficentMatrix(referenceBasisCoefficients, monomialValues, out);
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
        auto monomialValues = derivative(monomials)(in);
        multiplyWithCoefficentMatrix(referenceBasisCoefficients, monomialValues, out);
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
        auto monomialValues = derivative(derivative(monomials))(in);
        multiplyWithCoefficentMatrix(referenceBasisCoefficients, monomialValues, out);
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
        void_t>;

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
