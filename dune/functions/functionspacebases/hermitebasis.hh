// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HERMITEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HERMITEBASIS_HH

#include <algorithm>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <type_traits>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/lineartransformedlocalfiniteelement.hh>
#include <dune/functions/analyticfunctions/monomialset.hh>


namespace Dune
{
namespace Functions
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
  /**
   * \brief Multiply the evaluations of the monomialSet (see dune/functions/analyticfunctions/monomialset.hh) with a coefficient matrix, here FieldMatrix.
   * \tparam KCoeff The Field type of the coefficient matrix.
   * \tparam KMonom The range type of the monomials.
   * \tparam sizePolynom The number of polynomials to evaluate
   * \tparam sizeMonom The number of monomials handed as input.
   * \note This function also turns the return types from the dune-functions::DifferentiableFunction interface into those of the dune-localfunctions::LocalBasis interface.
   */
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

  struct DummyInterpolation{

    template<typename F, typename C>
    void interpolate(const F &f, std::vector<C> &out) const;
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
        Impl::DummyInterpolation>;

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

  namespace Impl
  {

    // Helper function returning an unordered range
    // of global indices associated to the element.
    // This could be implemented cheaper internally in
    // the MCMGMapper by storing a precomputed
    // container of all subsentities addressed by the layout.
    template<class GridView>
    auto subIndexSet(Dune::MultipleCodimMultipleGeomTypeMapper<GridView> const &mapper,
                    const typename GridView::template Codim<0>::Entity &element)
    {
      using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
      using Index = typename Mapper::Index;
      constexpr auto dimension = GridView::dimension;
      auto subIndices = std::vector<Index>();
      auto referenceElement = Dune::referenceElement<double, dimension>(element.type());
      for (auto codim : Dune::range(dimension + 1)) {
        for (auto subEntity : Dune::range(referenceElement.size(codim))) {
          std::size_t c = mapper.layout()(referenceElement.type(subEntity, codim), dimension);
          if (c > 0) {
            std::size_t firstIndex = mapper.subIndex(element, subEntity, codim);
            for (auto j : Dune::range(firstIndex, firstIndex + c)) {
              subIndices.push_back(j);
            }
          }
        }
      }
      return subIndices;
    }

    // Helper function computing an average mesh size per subentity
    // by averaging over the adjacent elements. This only considers
    // the subentities handled by the given mapper and returns a
    // vector of mesh sizes indixed according to the mapper.
    template<class FieldType = double, class Mapper>
    auto computeAverageSubEntityMeshSize(const Mapper& mapper)
    {
      constexpr auto dimension = Mapper::GridView::dimension;

      std::vector<unsigned int> adjacentElements(mapper.size(), 0);
      std::vector<FieldType> subEntityMeshSize(mapper.size(), 0.0);
      for(const auto& element : Dune::elements(mapper.gridView()))
      {
        auto A = element.geometry().volume();
        for(auto i : Impl::subIndexSet(mapper, element))
        {
          subEntityMeshSize[i] += A;
          ++(adjacentElements[i]);
        }
      }
      for(auto i : Dune::range(mapper.size()))
        subEntityMeshSize[i] = std::pow(subEntityMeshSize[i]/adjacentElements[i], 1./dimension);
      return subEntityMeshSize;
    }

    /** \brief Helper struct describing the Traits of the Global State of an Hermite Element*/
    template<class M, class R>
    struct HermiteGlobalStateTraits
    {
      using Mapper = M;
      using GlobalState = std::tuple<Mapper, std::vector<R>>;
      using LocalState = std::vector<R>;
    };


    /** \brief This class implements the Transformation for Hermite finite elements, that 'corrects' the non affine-equivalence.
    *   It is bindable to an Element, stateful and offers a transform(...) method.
    *   Its state can change upon binding and can be accessed via localState().
    *   Additionally, it offers a GlobalValuedInterpolation class to be used for interpolation.
    *   \tparam Element   The Grid Element
    *   \tparam R         The Fieldtype of the Finite Element
    *   \tparam reduced   If True, use the reduced Hermite element aka Discrete Kirchhoff Triangle
    */
    template<class Element, class R, class GlobalStateTraits, bool reduced = false>
    class HermiteTransformator
    {
      static constexpr int dim = Element::mydimension;
      static_assert(dim > 0 && dim < 4);
      static_assert(!(reduced && (dim != 2))); // TODO is there a reduced 3d version ?
      static constexpr std::size_t numberOfVertices = dim + 1;
      static constexpr std::size_t numberOfInnerDofs = reduced ? 0 : (dim - 1) * (dim - 1);
      static constexpr std::size_t numberOfVertexDofs = numberOfVertices * numberOfVertices;
    public:
      using GlobalState = typename GlobalStateTraits::GlobalState;
      using LocalState = typename GlobalStateTraits::LocalState;
      using size_type = std::size_t;

      HermiteTransformator(GlobalState const &globalState) : globalState_(&globalState) {}

      /** Binds the Transformator to an element.
      */
      void bind(Element const &e)
      {
        localState_ = std::apply(
            [&e](auto &&mapper, auto &&data) {
              LocalState localState;
              for (auto const &index : range(e.subEntities(dim)))
                localState.push_back(data[mapper.subIndex(e, index, dim)]);
              return localState;
            },
            *globalState_);

        fillMatrix(e.geometry(), localState_);
      }

      LocalState const &localState() const { return localState_; }

      //! The size of the transformed finite element.
      static constexpr size_type size()
      {
        if constexpr (dim == 1)
          return 4;
        else if constexpr (dim == 2) {
          if constexpr (reduced)
            return 9;
          else
            return 10;
        } else // dim == 3
          return 20;
      }

      /** Applies the transformation. Note that we do not distinguish for
      * Scalar/Vector/Matrix Type,
      * but only assume the Values to be Elements of a Vectorspace.
      * We assume random access containers. */
      template<typename InputValues, typename OutputValues, class LocalCoordinate>
      void transform(InputValues const &inValues, OutputValues &outValues, LocalCoordinate const& x) const
      {
        assert(inValues.size() == numberOfVertexDofs + numberOfInnerDofs);
        assert(reduced || (outValues.size() == inValues.size()));
        auto inIt = inValues.begin();
        auto outIt = outValues.begin();

        for (auto vertex : Dune::range(numberOfVertices)) {
          *outIt = *inIt; // value dof is not transformed
          outIt++, inIt++;
          // transform the gradient dofs together
          for (auto &&[row_i, i] : sparseRange(subMatrices_[vertex])) {
            outIt[i] = 0.;
            for (auto &&[val_i_j, j] : sparseRange(row_i))
              outIt[i] += val_i_j * inIt[j];
          }
          // increase pointer by size of gradient = dim
          outIt += dim, inIt += dim;
        }

        if constexpr (not reduced)
          // copy all remaining inner dofs
          std::copy(inIt, inValues.end(), outIt);
      }

    private:
      /**
      * \brief Fill the transformationmatrix m
      *
      * \tparam Geometry  the Geometry class
      * \param geometry   the geometry of the element we are bound to
      *
      *
      *      |1          0|} repeat for each vertex
      * m =  |  J/h       |}
      *      |0          1|} repeat for each inner dof i.e. (dim-1)^2 times
      *  where h is the mesh size average over the local vertex patch
      */
      template<class Geometry>
      void fillMatrix(Geometry const &geometry, LocalState const &averageSubEntityMeshSize)
      {
        auto const &refElement = Dune::ReferenceElements<typename Geometry::ctype, dim>::simplex();
        for (std::size_t i = 0; i < numberOfVertices; ++i) // dim + 1 vertices
        {
          subMatrices_[i] = geometry.jacobian(refElement.position(i, dim));
          subMatrices_[i] /= averageSubEntityMeshSize[i];
        }
      }

      // one transformation per vertex
      std::array<Dune::FieldMatrix<R, dim, dim>, numberOfVertices> subMatrices_;
      GlobalState const *globalState_;
      LocalState localState_;

    public:
      /**
      * \brief Class that evaluates the push forwards of the global nodes of a
      * LocalFunction. It stretches the LocalInterpolation interface, because we
      * evaluate the derivatives of f.
      *
      */
      class GlobalValuedInterpolation
      {
        using size_type = std::size_t;
        using ctype = typename Element::Geometry::ctype;
        static constexpr size_type numberOfVertices = dim + 1;
        static constexpr size_type innerDofCodim = (dim == 2) ? 0 : 1; // probably wrong for dim > 3
        static constexpr size_type numberOfInnerDofs =
            (dim - 1) * (dim - 1); // probably wrong for dim > 3

      public:
        GlobalValuedInterpolation(HermiteTransformator const &t) : transformator_(&t) {}

        /** \brief bind the Interpolation to an element and a localInterpolation.*/
        template<class LocalValuedLocalInterpolation>
        void bind([[maybe_unused]] Element const &element,
                  [[maybe_unused]] LocalValuedLocalInterpolation const &localInterpolation)
        {
          localState_ =&( transformator_->localState());
        }
      public:
        /** \brief Evaluate a given function and its derivatives at the nodes
        *
        * \tparam F Type of function to evaluate
        * \tparam C Type used for the values of the function
        * \param[in] f Function to evaluate
        * \param[out] out Array of function values
        */
        template<typename F, typename C>
        void interpolate(const F &f, std::vector<C> &out) const
        {
          auto df = derivative(f);
          out.resize(transformator_->size());

          auto const &refElement = Dune::ReferenceElements<ctype, dim>::simplex();
          // Iterate over vertices, dim dofs per vertex
          for (int i = 0; i < dim + 1; ++i) {
            auto x = refElement.position(i, dim);

            auto derivativeValue = df(x);
            out[i * numberOfVertices] = f(x);
            for (int d = 0; d < dim; ++d)
              out[i * numberOfVertices + d + 1] = getPartialDerivative(derivativeValue,d) * (*localState_)[i];
          }

          if constexpr (not reduced)
            for (size_type i = 0; i < numberOfInnerDofs; ++i) {
              out[numberOfVertices * numberOfVertices + i] =
                  f(refElement.position(i, innerDofCodim));
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


        HermiteTransformator const *transformator_;
        LocalState const* localState_;
      };

    };

  } // namespace Impl

  /**
  * \brief A pre-basis for a Hermitebasis
  *
  * \ingroup FunctionSpaceBasesImplementations
  *
  * \tparam GV  The grid view that the FE basis is defined on
  * \tparam R   Range type used for shape function values
  * \note This only works for simplex grids
  */
  template<typename GV, typename R, bool reduced = false>
  class HermitePreBasis : public LeafPreBasisMapperMixin<GV>
  {
    using Base = LeafPreBasisMapperMixin<GV>;
    using SubEntityMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;

    // helper methods to assign each subentity the number of dofs. Used by the LeafPreBasisMapperMixin.
    static constexpr auto cubicHermiteMapperLayout(Dune::GeometryType type, int gridDim)
    {
      if (type.isVertex())
        return 1 + gridDim; // one evaluation dof and gridDim derivative dofs per vertex
      if (gridDim == 1)    // in 1d there are no other dofs
        return 0;
      // in 2d we have one inner dof (i.e. on the triangle) or non for the reduced case
      // and in 3d we have one dof on each face (i.e. on each triangle)
      if ((type.isTriangle()) and (not reduced))
        return 1;
      else
        return 0; // this case is only entered for the interior of the 3d element. There are no dofs.
    }

  public:
    //! The grid view that the FE basis is defined on
    using GridView = GV;

    //! Type used for indices and size information
    using size_type = std::size_t;

  private:
    static const size_type dim = GV::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

    // the following typedefs configure the whole transformation framework
    //! Type used for the generic LocalFiniteElement
    using LocalFE = HermiteLocalFiniteElement<typename GridView::ctype, R, dim, reduced>;
    //! The Traits for the global state of the Hermite elements
    using GlobalStateTraits =
        Impl::HermiteGlobalStateTraits<Dune::MultipleCodimMultipleGeomTypeMapper<GridView>, R>;

    using GlobalState = typename GlobalStateTraits::GlobalState;
    //! Type used for the transformator which turns the generic LocalFiniteElement
    //! into a LocalFiniteElement that is affine equivalent to the global
    //! FiniteElement
    using HermiteTrafo = Impl::HermiteTransformator<Element, R, GlobalStateTraits, reduced>;

  public:
    //! Template mapping root tree path to type of created tree node
    using Node = Impl::TransformedNode<GridView, HermiteTrafo, LocalFE>;

    static constexpr size_type maxMultiIndexSize = 1;
    static constexpr size_type minMultiIndexSize = 1;
    static constexpr size_type multiIndexBufferSize = 1;

  public:
    //! Constructor for a given grid view object
    HermitePreBasis(const GV &gv)
        : Base(gv, cubicHermiteMapperLayout), globalState_({gv, mcmgVertexLayout()}, {})
    {
      updateState(gv);
      if (dim > 3)
        DUNE_THROW(Dune::NotImplemented, "HermitePreBasis only implemented for dim <= 3");
    }

    //! Update the stored grid view, to be called if the grid has changed
    void update(GridView const &gv)
    {
      Base::update(gv);
      updateState(gv);
    }

    /**
    * \brief Create tree node
    */
    Node makeNode() const { return Node(globalState_); }

  protected:
    void updateState(GridView const &gridView)
    {
      auto &[mapper, data] = globalState_;

      mapper.update(gridView);
      data = Impl::computeAverageSubEntityMeshSize<R>(mapper);
    }

    GlobalState globalState_;

  }; // class HermitePreBasis

  namespace BasisFactory
  {

    /**
    * \brief construct a PreBasisFactory for the full cubic Hermite Finite Element
    *
    * \tparam R RangeFieldType
    * \return the PreBasisFactory
    */
    template<typename R = double>
    auto hermite()
    {
      return [=](auto const &gridView) {
        return HermitePreBasis<std::decay_t<decltype(gridView)>, R>(gridView);
      };
    }

    /**
    * \brief construct a PreBasisFactory for the reduced cubic Hermite Finite Element
    *
    * \tparam R RangeFieldType
    * \return the PreBasisFactory
    */
    template<typename R = double>
    auto reducedHermite()
    {
      return [=](auto const &gridView) {
        return HermitePreBasis<std::decay_t<decltype(gridView)>, R, true>(gridView);
      };
    }

  } // namespace BasisFactory

} // namespace Functions
} // namespace Dune

#endif
