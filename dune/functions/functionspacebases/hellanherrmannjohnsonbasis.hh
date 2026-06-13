// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HELLANHERRMANNJOHNSONBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HELLANHERRMANNJOHNSONBASIS_HH

#include <algorithm>
#include <array>
#include <cassert>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/referencehelper.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/functions/common/innerproduct.hh>
#include <dune/functions/common/mapperutilities.hh>
#include <dune/functions/common/pullback.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/functionaldescriptor.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/transformed/bindcontext.hh>
#include <dune/functions/functionspacebases/transformed/derivative.hh>
#include <dune/functions/functionspacebases/transformed/localfiniteelement.hh>
#include <dune/functions/functionspacebases/transformed/piola.hh>

/**
 * \file hellanHerrmannjohnsonbasis.hh
 * \brief This file provides an implementation of the Hellan-Herrmann-Johnson finite element on triangles and tetrahedra.
 *
 * For reference, see [...].
 *
 * It contains in the following order:
 *     - A GlobalBasis typedef HellanHerrmannJohnsonBasis
 *     - A template HellanHerrmannJohnsonLocalFiniteElement providing an implementation
 *       of the LocalFiniteElement interface, along with its subparts (Impl namespace)
 *     - A template HellanHerrmannJohnsonNode
 *     - A template HellanHerrmannJohnsonPreBasis
 *     - Two factories hhj() and hellanHerrmannJohnson() in the BasisFactory namespace
 */
namespace Dune::Functions
{

  template<class GV, unsigned int k, class R>
  struct HellanHerrmannJohnsonPreBasis;

  /** \brief Nodal basis of a scalar cubic Hermite finite element space
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \note This only works for simplex grids.
   *
   * All arguments passed to the constructor will be forwarded to the constructor
   * of HellanHerrmannJohnsonPreBasis.
   *
   * \tparam GV The GridView that the space is defined on
   * \tparam k  The polynomial order of the element
   * \tparam R The range type of the local basis
   */
  template <class GV, unsigned int k, class R = double>
  using HellanHerrmannJohnsonBasis = DefaultGlobalBasis<HellanHerrmannJohnsonPreBasis<GV, k, R> >;

  namespace Impl
  {
    /** \brief Associations of the Hermite degrees of freedom to subentities of the
     * reference simplex
     *
     * \tparam dim Dimension of the reference simplex
     */
    template<int dim, unsigned int k>
    class HellanHerrmannJohnsonLocalCoefficients
    {
    public:
      using size_type = std::size_t;

      HellanHerrmannJohnsonLocalCoefficients()
        : localKeys_(size())
      {
        std::size_t idx = 0;
        if constexpr (dim == 2) {
          // edge DOFs
          int edgeSize = k+1;
          for (unsigned int s = 0; s < 3; ++s) // three edges
            for (int i = 0; i < edgeSize; ++i)
              localKeys_[idx++] = LocalKey(s,1,i);

          // cell DOFs
          int cellSize = ((k+1)*k)/2 * 3;
          for (unsigned int s = 0; s < 1; ++s) // one cell
            for (int i = 0; i < cellSize; ++i)
              localKeys_[idx++] = LocalKey(s,0,i);
        }
        else if constexpr (dim == 3) {
          // face DOFs
          int faceSize = (k+1)*(k+2)/2;
          for (unsigned int s = 0; s < 4; ++s) // four faces
            for (int i = 0; i < faceSize; ++i)
              localKeys_[idx++] = LocalKey(s,1,i);

          // cell DOFs
          int cellSize = (k+1)*(k+1)*(k+2);
          for (unsigned int s = 0; s < 1; ++s) // one cell
            for (int i = 0; i < cellSize; ++i)
              localKeys_[idx++] = LocalKey(s,0,i);
        }
      }

      /** \brief number of coefficients
       */
      static constexpr size_type size()
      {
        if constexpr (dim == 2)
          return 3*(k+1)*(k+2)/2;
        else if constexpr (dim == 3)
          return (k+1)*(k+2)*(k+3);
        else
          return 0;
      }

      /** \brief get i'th index
       */
      LocalKey const& localKey(size_type i) const
      {
        return localKeys_[i];
      }

    private:
      std::vector<LocalKey> localKeys_;
    };


    /** \brief Implementation of Hellan-Herrmann-Johnson basis function on the reference element
     * \tparam D Type to represent the field in the domain
     * \tparam R Type to represent the field in the range
     * \tparam dim Dimension of the domain simplex (limited to dim=2)
     * \tparam k The polynomial order of the basis
     */
    template<class D, class R, int dim, unsigned int k>
    class HellanHerrmannJohnsonReferenceLocalBasis
    {
    public:
      struct Traits {
        static constexpr int dimDomain = dim;
        static constexpr int dimRange = dim;
        using DomainFieldType = D;
        using DomainType = FieldVector<D,dim>;
        using RangeFieldType = R;
        using RangeType = FieldMatrix<R,dim,dim>;
        using DivDivType = R;
      };

    public:
      HellanHerrmannJohnsonReferenceLocalBasis()
      {
        static_assert(dim == 2 || dim == 3, "HellanHerrmannJohnsonReferenceLocalBasis only implemented for dim=2,3");
      }

      /** The number of basis functions in the basis
       */
      static constexpr unsigned int size()
      {
        return HellanHerrmannJohnsonLocalCoefficients<dim,k>::size();
      }

      /** The polynomial order of the basis
       */
      unsigned int order() const
      {
        return k;
      }

      /** \brief Evaluate function values of all shape functions at a given point
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Values of all shape functions at that point
       */
      void evaluateFunction(const typename Traits::DomainType& x,
                            std::vector<typename Traits::RangeType>& out) const;


      /** \brief Evaluate `div(div(phi))` of all shape functions `phi` at a given point
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Second derivative of all shape functions at that point
       */
      void evaluateDivDiv(const typename Traits::DomainType& x,
                          std::vector<typename Traits::DivDivType>& out) const;

    private:
      template <class Range>
      static Range sym(R a00, R a01, R a11)
      {
        return Range({{a00,a01},{a01,a11}});
      }

      template <class Range>
      static Range sym(R a00, R a01, R a02, R a11, R a12, R a22)
      {
        return Range({{a00,a01,a02},{a01,a11,a12},{a02,a12,a22}});
      }
    };

    /** \brief Reference interpolation for the HHJ degrees of freedom.
     *
     * The callable passed to interpolate() is reference-valued. Pulling a
     * physical tensor field back to that representation is handled separately
     * by HellanHerrmannJohnsonTransformation.
     */
    template<class D, int dim, unsigned int k>
    class HellanHerrmannJohnsonLocalInterpolation
    {
      using size_type = std::size_t;

      static constexpr unsigned int size()
      {
        return HellanHerrmannJohnsonLocalCoefficients<dim,k>::size();
      }

    public:

      /** \brief Evaluate a given function and its derivatives at the nodes
       *
       * \tparam F Type of function to evaluate
       * \tparam C Type used for the values of the function
       * \param[in] f Function to evaluate
       * \param[out] out Array of function values
       */
      template<class F, class C>
      void interpolate(const F& f, std::vector<C>& out) const
      {
        out.clear();
        out.resize(size(), C(0));

        if constexpr (dim == 2)
          interpolate2d(f, out);
        else if constexpr (dim == 3)
          interpolate3d(f, out);
        else
          DUNE_THROW(Dune::NotImplemented, "HHJ-interpolation only implemented for dim in {2,3}");
      }

    private:
      template<class F, class C>
      void interpolate2d(const F& local_f, std::vector<C>& out) const
      {
        auto const& refElem = ReferenceElements<D,dim>::simplex();
        auto const& edgeQuadRule = Dune::QuadratureRules<D,dim-1>::rule(refElem.type(0,1), 2*k+1);
        auto const& cellQuadRule = Dune::QuadratureRules<D,dim>::rule(refElem.type(), 2*k);

        using T = FieldMatrix<D,dim,dim>;
        std::array<T,3> directions {T({{0, 1}, {1, 0}}),
                                    T({{-2, 1}, {1, 0}}),
                                    T({{0, -1}, {-1, 2}})};

        if constexpr (k == 0) {
          for (int i = 0; i < refElem.size(1); ++i) {
            auto n = refElem.integrationOuterNormal(i); n/= n.two_norm();
            auto geoInCell = refElem.template geometry<1>(i);
            auto vol = geoInCell.volume();

            for (auto const& [x,w] : edgeQuadRule) {
              auto dx = geoInCell.integrationElement(x) * w;
              out[i] += Impl::pullback(local_f(geoInCell.global(x)),n) * vol * dx;
            }
          }
        }
        else {
          Dune::Impl::LagrangeSimplexLocalBasis<D,C, dim-1, k> edgeLagrangebasis;
          thread_local std::vector< typename Dune::Impl::LagrangeSimplexLocalBasis<D,C, dim-1, k>::Traits::RangeType> edgeValues;
          static constexpr std::size_t edgeSize = edgeLagrangebasis.size();

          Dune::Impl::LagrangeSimplexLocalBasis<D,C,dim, k - 1> cellLagrangebasis;
          thread_local std::vector< typename Dune::Impl::LagrangeSimplexLocalBasis<D,C,dim, k>::Traits::RangeType > cellValues;
          static constexpr std::size_t cellSize = cellLagrangebasis.size();

          for (int i = 0; i < refElem.size(1); ++i) {
            auto n = refElem.integrationOuterNormal(i); n/= n.two_norm();
            auto geoInCell = refElem.template geometry<1>(i);
            auto vol = geoInCell.volume();

            for (auto const& [x,w] : edgeQuadRule) {
              edgeLagrangebasis.evaluateFunction(x, edgeValues);
              auto dx = geoInCell.integrationElement(x) * w;
              auto nVn = Impl::pullback(local_f(geoInCell.global(x)),n) * vol * dx;
              for (std::size_t j = 0; j < edgeSize; j++)
                out[edgeSize*i + j] +=  edgeValues[j] * nVn;
            }
          }

          int startCell = refElem.size(1) * edgeSize;
          auto geoInCell = refElem.template geometry<0>(0);

          for (auto const& [x,w] : cellQuadRule) {
            cellLagrangebasis.evaluateFunction(x, cellValues);
            auto dx = geoInCell.integrationElement(x) * w;
            auto V = local_f(geoInCell.global(x));
            for (std::size_t j = 0; j < cellSize; ++j ) {
              for (std::size_t i = 0; i < directions.size(); ++i) {
                out[startCell + j*3 + i] += cellValues[j] * Impl::innerProduct(V,directions[i]) * dx;
              }
            }
          }
        }
      }

      template<class F, class C>
      void interpolate3d(const F& local_f, std::vector<C>& out) const
      {
        auto const& refElem = ReferenceElements<D,dim>::simplex();
        auto const& faceQuadRule = Dune::QuadratureRules<D,dim-1>::rule(refElem.type(0,1), 2*k+1);
        auto const& cellQuadRule = Dune::QuadratureRules<D,dim>::rule(refElem.type(), 2*k);

        using T = FieldMatrix<D,dim,dim>;
        std::array<T,4> directions1{T({{0, 1, 1}, {1, 0, 1}, {1, 1, 0}}),
                                    T({{-6, 1, 1}, {1, 0, 1}, {1, 1, 0}}),
                                    T({{0, 1, 1}, {1, -6, 1}, {1, 1, 0}}),
                                    T({{0, 1, 1}, {1,  0, 1}, {1, 1, -6}})};
        std::array<T,2> directions2{T({{0, 0, -1}, {0, 0, 1}, {-1, 1, 0}}),
                                    T({{0, -1, 0}, {-1, 0, 1}, {0, 1, 0}})};

        if constexpr (k == 0) {
          for (int i = 0; i < refElem.size(1); ++i) {
            auto n = refElem.integrationOuterNormal(i); n/= n.two_norm();
            auto geoInCell = refElem.template geometry<1>(i);
            auto vol = geoInCell.volume();

            for (auto const& [x,w] : faceQuadRule) {
              auto dx = geoInCell.integrationElement(x) * w;
              out[i] += Impl::pullback(local_f(geoInCell.global(x)),n) * vol * dx;
            }
          }

          int startCell = refElem.size(1);
          auto geoInCell = refElem.template geometry<0>(0);

          for (auto const& [x,w] : cellQuadRule) {
            auto dx = geoInCell.integrationElement(x) * w;
            auto V = local_f(geoInCell.global(x));
            for (std::size_t i = 0; i < directions2.size(); ++i) {
              out[startCell + i] += Impl::innerProduct(V,directions2[i]) * dx;
            }
          }
        }
        else {
          Dune::Impl::LagrangeSimplexLocalBasis<D,C, dim-1, k> faceLagrangebasis;
          thread_local std::vector<typename Dune::Impl::LagrangeSimplexLocalBasis<D,C, dim-1, k>::Traits::RangeType> faceValues;
          static constexpr std::size_t faceSize = faceLagrangebasis.size();

          for (int i = 0; i < refElem.size(1); ++i) {
            auto n = refElem.integrationOuterNormal(i); n/= n.two_norm();
            auto geoInCell = refElem.template geometry<1>(i);
            auto vol = geoInCell.volume();

            for (auto const& [x,w] : faceQuadRule) {
              faceLagrangebasis.evaluateFunction(x, faceValues);
              auto dx = geoInCell.integrationElement(x) * w;
              auto nVn = Impl::pullback(local_f(geoInCell.global(x)),n) * vol * dx;
              for (std::size_t j = 0; j < faceSize; j++)
                out[faceSize*i + j] += faceValues[j] * nVn;
            }
          }

          Dune::Impl::LagrangeSimplexLocalBasis<D,C,dim, k - 1> cellLagrangebasis1;
          thread_local std::vector<typename Dune::Impl::LagrangeSimplexLocalBasis<D,C,dim, k-1>::Traits::RangeType > cellValues1;
          static constexpr std::size_t cellSize1 = cellLagrangebasis1.size();

          Dune::Impl::LagrangeSimplexLocalBasis<D,C,dim, k> cellLagrangebasis2;
          thread_local std::vector<typename Dune::Impl::LagrangeSimplexLocalBasis<D,C,dim, k>::Traits::RangeType > cellValues2;
          static constexpr std::size_t cellSize2 = cellLagrangebasis2.size();

          int startCell1 = refElem.size(1) * faceSize;
          int n1 = directions1.size();
          int startCell2 = startCell1 + cellSize1 * n1;
          int n2 = directions2.size();
          auto geoInCell = refElem.template geometry<0>(0);

          for (auto const& [x,w] : cellQuadRule) {
            cellLagrangebasis1.evaluateFunction(x, cellValues1);
            cellLagrangebasis2.evaluateFunction(x, cellValues2);
            auto dx = geoInCell.integrationElement(x) * w;
            auto V = local_f(geoInCell.global(x));
            for (std::size_t j = 0; j < cellSize1; ++j) {
              for (std::size_t i = 0; i < directions1.size(); ++i) {
                out[startCell1 + n1*j + i] += cellValues1[j] * Impl::innerProduct(V,directions1[i]) * dx;
              }
            }
            for (std::size_t j = 0; j < cellSize2; ++j) {
              for (std::size_t i = 0; i < directions2.size(); ++i) {
                out[startCell2 + n2*j + i] += cellValues2[j] * Impl::innerProduct(V,directions1[i]) * dx;
              }
            }
          }
        }
      }

    };


    /** \brief Hellan-Herrmann-Johnson finite element on the reference simplex.
     * For more Details, see <dune/functions/functionspacebases/hellanherrmannjohnsonbasis.hh>.
     *
     * \tparam D Type used for domain coordinates
     * \tparam R Type used for function values
     * \tparam dim dimension of the reference element
     */
    template<class D, class R, int dim, unsigned int k>
    class HellanHerrmannJohnsonReferenceLocalFiniteElement
    {
    public:
      HellanHerrmannJohnsonReferenceLocalFiniteElement()
      {
        static_assert(dim==2 || dim==3,
          "HellanHerrmannJohnsonReferenceLocalFiniteElement only implemented for dim=2,3");
      }

      /** \brief Export number types, dimensions, etc.
       */
      using size_type = std::size_t;
      using Traits = LocalFiniteElementTraits<
          HellanHerrmannJohnsonReferenceLocalBasis<D, R, dim, k>,
          HellanHerrmannJohnsonLocalCoefficients<dim, k>,
          HellanHerrmannJohnsonLocalInterpolation<D, dim, k>>;


      /** \brief Returns object that evaluates degrees of freedom
       */
      const typename Traits::LocalBasisType& localBasis() const
      {
        return basis_;
      }

      /** \brief Returns the assignment of the degrees of freedom to the element
       * subentities
       */
      const typename Traits::LocalCoefficientsType& localCoefficients() const
      {
        return coefficients_;
      }

      /** \brief Returns object that evaluates degrees of freedom
       */
      const typename Traits::LocalInterpolationType& localInterpolation() const
      {
        return interpolation_;
      }

      /** \brief The reference element that the local finite element is defined on
       */
      static constexpr GeometryType type()
      {
        return GeometryTypes::simplex(dim);
      }

      /** The size of the transformed finite element.
       */
      static constexpr size_type size()
      {
        return HellanHerrmannJohnsonLocalCoefficients<dim,k>::size();
      }

    private:
      typename Traits::LocalBasisType basis_;
      typename Traits::LocalCoefficientsType coefficients_;
      typename Traits::LocalInterpolationType interpolation_;
    };

  } // end namespace Impl



  // *****************************************************************************
  // This is the reusable part of the basis. It contains
  //
  //   HellanHerrmannJohnsonPreBasis
  //   HellanHerrmannJohnsonNode
  //
  // The pre-basis allows to create the others and is the owner of possible shared
  // state. These components do _not_ depend on the global basis and local view
  // and can be used without a global basis.
  // *****************************************************************************

  template<class GV, unsigned int k, class R>
  class HellanHerrmannJohnsonNode
    : public LeafBasisNode
  {
  public:
    using size_type = std::size_t;
    using Element = typename GV::template Codim<0>::Entity;
    using Context = ElementBindContext<Element>;
    using ReferenceFiniteElement = Impl::HellanHerrmannJohnsonReferenceLocalFiniteElement<
      typename GV::ctype,R,GV::dimension,k>;
    using Transformation = DoubleContravariantTransformation<typename Element::Geometry>;
    using FiniteElement = TransformedLocalFiniteElement<ReferenceFiniteElement, Context,
      Transformation, LocalBasisMode::physical, LocalInterpolationMode::transformed>;

    HellanHerrmannJohnsonNode() = default;

    //! Return current element, throw if unbound
    Element const& element() const
    {
      return context_.element();
    }

    /** \brief Return the LocalFiniteElement for the element we are bound to
     *
     * The LocalFiniteElement implements the corresponding interfaces of the
     * dune-localfunctions module
     */
    FiniteElement const& finiteElement() const
    {
      return finiteElement_;
    }

    //! Bind to element.
    void bind(Element const& e)
    {
      context_.bind(e);
      finiteElement_.bind(referenceFiniteElement_,context_);
      this->setSize(finiteElement_.size());
    }

    //! The order of the local basis.
    unsigned int order() const { return k; }

  protected:
    Context context_;
    ReferenceFiniteElement referenceFiniteElement_;
    FiniteElement finiteElement_;
  };


  /**
   * \brief A pre-basis for a Hellan-Herrmann-Johnson
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \tparam GV  The grid view that the FE basis is defined on
   * \tparam R   Range type used for shape function values
   * \note This only works for simplex grids
   */
  template<class GV, unsigned int k, class R>
  class HellanHerrmannJohnsonPreBasis
    : public LeafPreBasisMixin< HellanHerrmannJohnsonPreBasis<GV,k,R> >
  {
    using Base = LeafPreBasisMixin< HellanHerrmannJohnsonPreBasis<GV,k,R> >;
    using Element = typename GV::template Codim<0>::Entity;
    using D = typename GV::ctype;
    static const int dim = GV::dimension;

  private:
    static constexpr int edgeSize = (k+1);
    static constexpr int faceSize = (k+1)*(k+2)/2;
    static constexpr int cellSize = (dim == 2 ? k*(k+1)/2 * 3 : (k+2)*(k+1)*(k+1));

    // helper methods to assign each subentity the number of dofs. Used by the LeafPreBasisMapperMixin.
    static constexpr std::size_t hellanHerrmannJohnsonLayout(Dune::GeometryType type, int gridDim)
    {
      // currently only implemented for dim == 2
      if (type.isLine() && gridDim == 2)
        return edgeSize;
      else if (type.isTriangle() && gridDim == 3)
        return faceSize;
      else if (int(type.dim()) == gridDim)
        return cellSize;
      else
        return 0;
    }

    using FaceDOFPermutation = Experimental::LagrangeFaceDOFPermutation<typename GV::Grid::GlobalIdSet>;

  public:

    //! The grid view that the FE basis is defined on
    using GridView = GV;

    //! Type used for indices and size information
    using size_type = std::size_t;

    //! Template mapping root tree path to type of created tree node
    using Node = HellanHerrmannJohnsonNode<GridView, k, R>;

  public:

    //! Constructor for a given grid view object
    HellanHerrmannJohnsonPreBasis(const GridView& gv)
      : gridView_(gv)
      , mapper_(gridView_, hellanHerrmannJohnsonLayout)
      , faceDOFPermutation_(gridView_.grid().globalIdSet(), k+2)
    {
      static_assert(dim==2 || dim==3, "HellanHerrmannJohnsonPreBasis only implemented for dim=2,3");
    }

    //! Initialize the global indices
    void initializeIndices()
    {}

    //! Obtain the grid view that the basis is defined on
    const GridView& gridView() const
    {
      return gridView_;
    }

    //! Update the stored grid view, to be called if the grid has changed
    void update (const GridView& gv)
    {
      gridView_ = gv;
      mapper_.update(gridView_);
      faceDOFPermutation_ = FaceDOFPermutation(gridView_.grid().globalIdSet(), k+2);
    }

    //! Create tree node
    Node makeNode() const
    {
      return Node{};
    }

    //! Get the total dimension of the space spanned by this basis
    size_type dimension() const
    {
      return mapper_.size();
    }

    //! Get the maximal number of DOFs associated to node for any element
    size_type maxNodeSize() const
    {
      constexpr unsigned int dim = GV::dimension;
      return dim*(dim+1)/2 * binomial(k+dim,dim);
    }

    //! Fill the range referenced by the iterator `it` with the global indices
    template<class Node, class It>
    It indices(const Node& node, It it) const
    {
      const auto& element = node.element();
      const auto& localCoefficients = node.finiteElement().localCoefficients();

      // Precompute orientations of all faces
      auto faceOrientations = faceDOFPermutation_.computeFaceOrientations(element);
      for(auto localIndex : Dune::range(localCoefficients.size()))
      {
        auto localKey = localCoefficients.localKey(localIndex);
        auto globalIndex = mapper_.subIndex(element, localKey.subEntity(), localKey.codim());
        globalIndex += faceDOFPermutation_.permuteFaceDOF(localKey, faceOrientations);
        *it = {{ (size_type)(globalIndex) }};
        ++it;
      }
      return it;
    }

  protected:
    GridView gridView_;

    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper_;
    FaceDOFPermutation faceDOFPermutation_;
  };


  namespace BasisFactory
  {

    /**
     * \brief construct a PreBasisFactory for the Hellan-Herrmann-Johnson Finite Element
     *
     * \tparam k  The polynomial order of the basis
     * \tparam R  Type of the basis range
     * \return a factory function to create the HellanHerrmannJohnsonPreBasis
     * \relates HellanHerrmannJohnsonBasis
     */
    template<unsigned int k, class R = double>
    auto hhj()
    {
      return []<class GV>(GV const &gridView) {
        return HellanHerrmannJohnsonPreBasis<GV, k, R>(gridView);
      };
    }

    /**
     * \brief construct a PreBasisFactory for the Hellan-Herrmann-Johnson Finite Element
     *
     * \tparam k  The polynomial order of the basis
     * \tparam R  Type of the basis range
     * \return a factory function to create the HellanHerrmannJohnsonPreBasis
     * \relates HellanHerrmannJohnsonBasis
     */
    template<unsigned int k, class R = double>
    auto hellanHerrmannJohnson()
    {
      return []<class GV>(GV const &gridView) {
        return HellanHerrmannJohnsonPreBasis<GV, k, R>(gridView);
      };
    }

  } // end namespace BasisFactory
} // end namespace Dune::Functions
#include <dune/functions/functionspacebases/hellanherrmannjohnsonbasis_inc.hh>
#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HELLANHerrmannJOHNSONBASIS_HH
