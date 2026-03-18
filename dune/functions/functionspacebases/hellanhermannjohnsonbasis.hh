// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HELLANHERRMANNJOHNSONBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HELLANHERRMANNJOHNSONBASIS_HH

#include <algorithm>
#include <array>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/functions/common/innerproduct.hh>
#include <dune/functions/common/mapperutilities.hh>
#include <dune/functions/common/multidot.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/functionaldescriptor.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/nodes.hh>

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
        static_assert(dim == 2, "HellanHerrmannJohnsonLocalCoefficients only implemented for dim=2");

        std::size_t idx = 0;
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

      /** \brief number of coefficients
       */
      static constexpr size_type size()
      {
        if constexpr (dim==2)
          return 3*(k+1)*(k+2)/2;
        else
          return 0;
      }

      /** \brief get i'th index
       */
      LocalKey const &localKey(size_type i) const
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
    private:
      using Range = typename Traits::RangeType;

    public:
      HellanHerrmannJohnsonReferenceLocalBasis()
      {
        static_assert(dim == 2, "HellanHerrmannJohnsonReferenceLocalBasis only implemented for dim=2");
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
    };

    /**
    * \brief Implementation of a dune-localfunctions LocalBasis that applies a
    * linear basis transformation
    *
    * \tparam FEImplementation The finite element implementation
    * \tparam ReferenceLocalBasisTraits LocalBasisTraits of the reference local basis
    */
    template<class E, class R, unsigned int k>
    class HellanHerrmannJohnsonLocalBasis
    {
      using Element = E;
      using Geometry = typename E::Geometry;
      static constexpr int dim = Geometry::mydimension;
      using D = typename Geometry::ctype;

      using ReferenceLocalBasis = HellanHerrmannJohnsonReferenceLocalBasis<D, R, dim, k>;

      public:
        struct Traits {
          static constexpr int dimDomain = dim;
          static constexpr int dimRange = Geometry::coorddimension;
          using DomainFieldType = D;
          using DomainType = FieldVector<D,dimDomain>;
          using RangeFieldType = R;
          using RangeType = FieldMatrix<R,dimRange,dimRange>;
          using DivDivType = R;
        };

      public:
        /**
        * \brief Number of shape functions
        * This need not to be equal to the size of the reference local basis
        */
        auto size() const
        {
          return refLocalBasis_.size();
        }

        void bind(Element const& element)
        {
          geometry_.emplace(element.geometry());
        }

        //! \brief Evaluate all shape functions
        void evaluateFunction(const typename Traits::DomainType& x,
                              std::vector<typename Traits::RangeType>& out) const
        {
          refLocalBasis_.evaluateFunction(x, rangeBuffer_);
          out.resize(size());
          transformRange(x, rangeBuffer_, out);
        }

        /**
        * \brief Evaluate div(div(phi)) of all shape functions phi
        *
        * \param x Point in the reference element where to evaluation the derivatives
        * \param[out] out The div(div(phi)) of all shape functions phi at the point x
        */
        void evaluateDivDiv(const typename Traits::DomainType& x,
                            std::vector<typename Traits::DivDivType>& out) const
        {
          refLocalBasis_.evaluateDivDiv(x, divDivBuffer_);
          out.resize(size());
          transformDivDiv(x, divDivBuffer_, out);
        }

        //! \brief Polynomial order of the shape functions
        auto order() const { return refLocalBasis_.order(); }

      protected:
        // contravariant Piola transform
        void transformRange(typename Traits::DomainType const& x,
              std::vector<typename ReferenceLocalBasis::Traits::RangeType> const& inValues,
              std::vector<typename Traits::RangeType>& outValues) const
        {
          assert(inValues.size() == size());
          assert(outValues.size() == inValues.size());
          //...

          auto Jt = geometry_->jacobianTransposed(x);
          auto dx = geometry_->integrationElement(x);

          for (std::size_t i = 0; i < inValues.size(); ++i)
          {
            outValues[i] = multiDot(inValues[i],Jt,Jt);
            outValues[i] /= dx*dx;
          }
        }

        void transformDivDiv(typename Traits::DomainType const& x,
              std::vector<typename ReferenceLocalBasis::Traits::DivDivType> const& inValues,
              std::vector<typename Traits::DivDivType>& outValues) const
        {
          assert(inValues.size() == size());
          assert(outValues.size() == inValues.size());

          auto dx = geometry_->integrationElement(x);

          for (std::size_t i = 0; i < inValues.size(); ++i)
          {
            outValues[i] = dx * inValues[i];
          }
        }


      private:
        ReferenceLocalBasis refLocalBasis_ = {};
        std::optional<Geometry> geometry_ = std::nullopt;
        mutable std::vector<typename ReferenceLocalBasis::Traits::RangeType> rangeBuffer_ = {};
        mutable std::vector<typename ReferenceLocalBasis::Traits::DivDivType> divDivBuffer_ = {};
    };


    /** \brief Class that evaluates the push forwards of the global nodes of a
     * LocalFunction. It stretches the LocalInterpolation interface, because we
     * evaluate the derivatives of f.
     *
     */
    template<class E, unsigned int k>
    class HellanHerrmannJohnsonLocalInterpolation
    {
      using Element = E;
      using Geometry = typename E::Geometry;
      static constexpr int dim = Geometry::mydimension;
      static constexpr int dimRange = Geometry::coorddimension;

      using D = typename Geometry::ctype;

      using size_type = std::size_t;

      static constexpr unsigned int size()
      {
        return HellanHerrmannJohnsonLocalCoefficients<dim,k>::size();
      }

    public:

      HellanHerrmannJohnsonLocalInterpolation()
      {}

      /** \brief bind the Interpolation to an element and a local state.
       */
      void bind(Element const& element)
      {
        geometry_.emplace(element.geometry());
      }

      template <class F>
      struct LocalValuedFunction
      {
        LocalValuedFunction(const F& f, const Geometry& geometry)
          : f_(f), geometry_(geometry)
        {}

        // Apply the inverse Piola transform
        auto operator() (const typename Geometry::LocalCoordinate& xi) const
        {
          auto Jit = geometry_.jacobianInverseTransposed(xi);
          auto dx = geometry_.integrationElement(xi);
          return multiDot(f_(xi), Jit, Jit) * (dx*dx);
        }

        F const& f_;
        Geometry const& geometry_;
      };

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

        auto refElem = referenceElement(*geometry_);

        auto local_f = LocalValuedFunction{f, *geometry_};
        auto const& edgeQuadRule = Dune::QuadratureRules<D,dim-1>::rule(refElem.type(0,1), 2*k+1);
        auto const& cellQuadRule = Dune::QuadratureRules<D,dim>::rule(refElem.type(), 2*k);

        static_assert(dim == 2, "HHJ-interpolation only implemented for dim == 2");
        using T = FieldMatrix<D,dim,dim>;
        std::array<T,3> directions {T({{0, 1}, {1, 0}}),
                                    T({{-2, 1}, {1, 0}}),
                                    T({{0, -1}, {-1, 2}})};

        // 1. integral moments over the edges
        if constexpr (k == 0) {
          for (int i = 0; i < refElem.size(1); ++i) {
            auto n = refElem.integrationOuterNormal(i); n/= n.two_norm();
            auto geoInCell = refElem.template geometry<1>(i);
            auto vol = geoInCell.volume();

            for (auto const& [x,w] : edgeQuadRule) {
              auto dx = geoInCell.integrationElement(x) * w;
              out[i] += multiDot(local_f(geoInCell.global(x)),n,n) * vol * dx;
            }
          }
        }
        else{
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
              auto nVn = multiDot(local_f(geoInCell.global(x)),n,n) * vol * dx;
              for (std::size_t j = 0; j < edgeSize; j++)
                out[edgeSize*i + j] +=  edgeValues[j] * nVn;
            }
          }

          static int startCell = refElem.size(1) * edgeSize;
          auto geoInCell = refElem.template geometry<0>(0);

          for (auto const& [x,w] : cellQuadRule) {
            cellLagrangebasis.evaluateFunction(x, cellValues);
            auto dx = geoInCell.integrationElement(x) * w;
            for (std::size_t j = 0; j < cellSize; ++j ){
              for (std::size_t i = 0; i < directions.size(); ++i) {
                out[startCell + j*3 + i] += cellValues[j] * innerProduct(local_f(geoInCell.global(x)),directions[i]) * dx;
              }
            }
          }
        }
      }

    private:
      std::optional<Geometry> geometry_ = std::nullopt;
    };


    /** \brief Hellan-Herrmann-Johnson finite element for simplices, as defined on the reference Element.
     * For more Details, see <dune/functions/functionspacebases/hellanHerrmannjohnsonbasis.hh>.
     *
     * \tparam D Type used for domain coordinates
     * \tparam R Type used for function values
     * \tparam dim dimension of the reference element
     */
    template<class E, class R, unsigned int k>
    class HellanHerrmannJohnsonLocalFiniteElement
    {
      using Element = E;
      using Geometry = typename E::Geometry;
      using D = typename Geometry::ctype;
      static constexpr int dim = Geometry::mydimension;

    public:
      HellanHerrmannJohnsonLocalFiniteElement()
      {
        static_assert(dim==2, "HellanHerrmannJohnsonLocalFiniteElement only implemented for dim=2");
      }

      /** \brief Export number types, dimensions, etc.
       */
      using size_type = std::size_t;
      using Traits = LocalFiniteElementTraits<
          HellanHerrmannJohnsonLocalBasis<E, R, k>,
          HellanHerrmannJohnsonLocalCoefficients<dim, k>,
          HellanHerrmannJohnsonLocalInterpolation<E, k>>;


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

      /** Binds the Finite Element to an element.
       */
      void bind(Element const& e)
      {
        basis_.bind(e);
        interpolation_.bind(e);
      }

    protected: // implementation of the expected interface of TransformedLocalBasis

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
    using FiniteElement = Impl::HellanHerrmannJohnsonLocalFiniteElement<Element, R, k>;

    HellanHerrmannJohnsonNode()
      : element_(nullptr)
    {}

    //! Return current element, throw if unbound
    Element const& element() const
    {
      return *element_;
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
      element_ = &e;
      finiteElement_.bind(*element_);
      this->setSize(finiteElement_.size());
    }

    //! The order of the local basis.
    unsigned int order() const { return k; }

  protected:
    FiniteElement finiteElement_;
    Element const* element_;
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
    : public LeafPreBasisMapperMixin<GV, Impl::EdgeTwist<typename GV::IndexSet>>
  {
    using Twist = Impl::EdgeTwist<typename GV::IndexSet>;
    using Base = LeafPreBasisMapperMixin<GV, Twist>;
    using Element = typename GV::template Codim<0>::Entity;
    using D = typename GV::ctype;
    static const std::size_t dim = GV::dimension;

    // helper methods to assign each subentity the number of dofs. Used by the LeafPreBasisMapperMixin.
    static constexpr std::size_t hellanHerrmannJohnsonLayout(Dune::GeometryType type, int gridDim)
    {
      if constexpr(k == 0)
        return type.isLine() ? 1 : 0;
      else if constexpr(k == 1)
        return type.isLine() ? 2 : int(type.dim()) == gridDim ? 3 : 0;
      else if constexpr(k == 2)
        return type.isLine() ? 3 : int(type.dim()) == gridDim ? 9 : 0;
      else{

        return type.isLine() ? edgeSize : int(type.dim()) == gridDim ? cellSize : 0;
      }
    }
  private:
    static constexpr int edgeSize = (k+1);
    static constexpr int cellSize = ((k+1)*k)/2*3;
  public:
    //! The grid view that the FE basis is defined on
    using GridView = GV;

    //! Type used for indices and size information
    using size_type = std::size_t;

    //! Template mapping root tree path to type of created tree node
    using Node = HellanHerrmannJohnsonNode<GridView, k, R>;

  public:

    //! Constructor for a given grid view object
    HellanHerrmannJohnsonPreBasis(const GV &gv)
      : Base(gv, hellanHerrmannJohnsonLayout, Twist{gv.indexSet(), hellanHerrmannJohnsonLayout(GeometryTypes::line,dim)})
    {
      static_assert(dim==2, "HellanHerrmannJohnsonPreBasis only implemented for dim=2");
    }

    //! Update the stored grid view, to be called if the grid has changed
    void update(GridView const &gv)
    {
      Base::update(gv);
    }

    /**
     * \brief Create tree node
     */
    Node makeNode() const
    {
      return Node{};
    }
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
#include <dune/functions/functionspacebases/hellanhermannjohnsonbasis_inc.hh>
#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HELLANHerrmannJOHNSONBASIS_HH
