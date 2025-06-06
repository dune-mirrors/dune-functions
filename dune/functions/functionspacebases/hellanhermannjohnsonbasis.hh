// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HELLANHERMANNJOHNSONBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HELLANHERMANNJOHNSONBASIS_HH

#include <algorithm>
#include <array>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/densetensor.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/tensordot.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/functions/common/mapperutilities.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/functionaldescriptor.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/nodes.hh>

/**
 * \file hellanhermannjohnsonbasis.hh
 * \brief This file provides an implementation of the Hellan-Hermann-Johnson finite element on triangles and tetrahedra.
 *
 * For reference, see [...].
 *
 * It contains in the following order:
 *     - A GlobalBasis typedef HellanHermannJohnsonBasis
 *     - A template HellanHermannJohnsonLocalFiniteElement providing an implementation
 *       of the LocalFiniteElement interface, along with its subparts (Impl namespace)
 *     - A template HellanHermannJohnsonNode
 *     - A template HellanHermannJohnsonPreBasis
 *     - Two factories hhj() and hellanHermannJohnson() in the BasisFactory namespace
 */
namespace Dune::Functions
{

  template<class GV, unsigned int k, class R>
  struct HellanHermannJohnsonPreBasis;

  /** \brief Nodal basis of a scalar cubic Hermite finite element space
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \note This only works for simplex grids. The Hermite basis is only implemented for degree 3.
   * \note The Hermite Finite element has the following properties:
   *   - Its global space is in \f$ C^1 \f$ if 1d otherwise in \f$ H^1 \f$.
   *   - The reduced version (only 2d) is part of the Discrete Kirchhoff Triangle.
   *   - Its interpolation evaluates derivatives, i.e. you cannot interpolate into a lambda function.
   *   - Strongly enforcing boundary conditions is not as simple as with langrange bases
   *   - It global space is not nested, i.e. the space on a refined grid is not a subspace of the
   *     space on the coarser grid.
   * All arguments passed to the constructor will be forwarded to the constructor
   * of HellanHermannJohnsonPreBasis.
   *
   * \tparam GV The GridView that the space is defined on
   * \tparam k  The polynomial order of the element
   * \tparam R The range type of the local basis
   */
  template <class GV, unsigned int k, class R = double>
  using HellanHermannJohnsonBasis = DefaultGlobalBasis<HellanHermannJohnsonPreBasis<GV, k, R> >;

  namespace Impl
  {
    /** \brief Associations of the Hermite degrees of freedom to subentities of the
     * reference simplex
     *
     * \tparam dim Dimension of the reference simplex
     */
    template<int dim, unsigned int k>
    class HellanHermannJohnsonLocalCoefficients
    {
    public:
      using size_type = std::size_t;

      HellanHermannJohnsonLocalCoefficients()
        : localKeys_(size())
      {
        static_assert(dim == 2, "HellanHermannJohnsonLocalCoefficients only implemented for dim=2");

        if constexpr(k == 0) {
          localKeys_[0] = LocalKey(0,1,0);
          localKeys_[1] = LocalKey(1,1,0);
          localKeys_[2] = LocalKey(2,1,0);
        } else if constexpr(k == 1) {
          localKeys_[0] = LocalKey(0,1,0);
          localKeys_[1] = LocalKey(0,1,1);
          localKeys_[2] = LocalKey(1,1,0);
          localKeys_[3] = LocalKey(1,1,1);
          localKeys_[4] = LocalKey(2,1,0);
          localKeys_[5] = LocalKey(2,1,1);
          localKeys_[6] = LocalKey(0,0,0);
          localKeys_[7] = LocalKey(0,0,1);
          localKeys_[8] = LocalKey(0,0,2);
        } else if constexpr(k == 2) {
          std::size_t idx = 0;
          // edge DOFs
          for (unsigned int s = 0; s < 3; ++s)
            for (unsigned int i = 0; i < 3; ++i)
              localKeys_[idx++] = LocalKey(s,1,i);

          // cell DOFs
          for (unsigned int s = 0; s < 1; ++s)
            for (unsigned int i = 0; i < 9; ++i)
              localKeys_[idx++] = LocalKey(s,0,i);
        }
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


    /** \brief Implementation of Hellan-Hermann-Johnson basis function on the reference element
     * \tparam D Type to represent the field in the domain
     * \tparam R Type to represent the field in the range
     * \tparam dim Dimension of the domain simplex (limited to dim=2)
     * \tparam k The polynomial order of the basis
     */
    template<class D, class R, int dim, unsigned int k>
    class HellanHermannJohnsonReferenceLocalBasis
    {
    public:
      struct Traits {
        static constexpr int dimDomain = dim;
        static constexpr int dimRange = dim;
        using DomainFieldType = D;
        using DomainType = FieldVector<D,dim>;
        using RangeFieldType = R;
        using RangeType = DenseTensor<R,dim,dim>;
        using DivDivType = R;
      };

    public:
      HellanHermannJohnsonReferenceLocalBasis()
      {
        static_assert(dim == 2, "HellanHermannJohnsonReferenceLocalBasis only implemented for dim=2");
      }

      /** The number of basis functions in the basis
       */
      static constexpr unsigned int size()
      {
        return HellanHermannJohnsonLocalCoefficients<dim,k>::size();
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
                            std::vector<typename Traits::RangeType>& out) const
      {
        out.resize(size());
        // ...
        using Range = typename Traits::RangeType;
        if constexpr(k == 0) {
          R half = -R(1)/2;
          out[0] = Range({{0,-half},{-half,1}});
          out[1] = Range({{1,-half},{-half,0}});
          out[2] = Range({{0, half},{ half,0}});
        }
        else if constexpr(k == 1) {
          out[0] = sym<Range>(0, 3*x[0]+3*x[1]-2, -6*x[0]-6*x[1]+4);
          out[1] = sym<Range>(0, 1-3*x[0], 6*x[0]-2);
          out[2] = sym<Range>(-6*x[0]-6*x[1]+4, 3*x[0]+3*x[1]-2, 0);
          out[3] = sym<Range>(6*x[1]-2, 1-3*x[1], 0);
          out[4] = sym<Range>(0, 3*x[0]-1, 0);
          out[5] = sym<Range>(0, 3*x[1]-1, 0);
          out[6] = sym<Range>(3*x[0], -15*x[0]/2-15*x[1]/2+6, 3*x[1]);
          out[7] = sym<Range>(-3*x[0], 3*x[0]+3*x[1]/2-R(3)/2, 0);
          out[8] = sym<Range>(0, -3*x[0]/2-3*x[1]+R(3)/2, 3*x[1]);
        }
        else if constexpr(k == 2) {
          R xx = x[0]*x[0];
          R xy = x[0]*x[1];
          R yy = x[1]*x[1];

          // edge 0
          out[0] = sym<Range>(0, -15*xx - 30*xy + 18*x[0] - 15*yy + 18*x[1] - R(9)/2,
                              30*xx + 60*xy - 36*x[0] + 30*yy - 36*x[1] + 9);
          out[1] = sym<Range>(0, -15*xx + 12*x[0] - R(3)/2, 30*xx - 24*x[0] + 3);
          out[2] = sym<Range>(0, 15*xx/2 + 15*xy/2 - 15*x[0]/2 - 15*yy/4 + 3*x[1]/2 + R(3)/4,
                              -15*xx - 15*xy + 15*x[0] + 15*yy/2 - 3*x[1] - R(3)/2);

          // edge 1
          out[3] = sym<Range>(30*xx + 60*xy - 36*x[0] + 30*yy - 36*x[1] + 9, -15*xx - 30*xy + 18*x[0] - 15*yy + 18*x[1] - R(9)/2, 0);
          out[4] = sym<Range>(30*yy - 24*x[1] + 3, -15*yy + 12*x[1] - R(3)/2, 0);
          out[5] = sym<Range>(15*xx/2 - 15*xy - 3*x[0] - 15*yy + 15*x[1] - R(3)/2, -15*xx/4 + 15*xy/2 + 3*x[0]/2 + 15*yy/2 - 15*x[1]/2 + R(3)/4, 0);

          // edge 2
          out[6] = sym<Range>(0, 15*xx - 12*x[0] + R(3)/2, 0);
          out[7] = sym<Range>(0, 15*yy - 12*x[1] + R(3)/2, 0);
          out[8] = sym<Range>(0, 15*xx/4 + 15*xy - 6*x[0] + 15*yy/4 - 6*x[1] + R(3)/2, 0);

          // cell
          out[9] = sym<Range>(-60*xx - 60*xy + 48*x[0], 90*xx + 180*xy - 120*x[0] + 90*yy - 120*x[1] + 36, -60*xy - 60*yy + 48*x[1]);
          out[10] = sym<Range>(60*xx + 60*xy - 48*x[0], -45*xx - 60*xy + 48*x[0] - 15*yy + 24*x[1] - 9, 0);
          out[11] = sym<Range>(0, 15*xx + 60*xy - 24*x[0] + 45*yy - 48*x[1] + 9, -60*xy - 60*yy + 48*x[1]);
          out[12] = sym<Range>(30*xx - 12*x[0], -135*xx - 150*xy + 150*x[0] + 30*x[1] - 24, 60*xy - 12*x[1]);
          out[13] = sym<Range>(-30*xx + 12*x[0], 45*xx + 30*xy - 42*x[0] - 6*x[1] + 6, 0);
          out[14] = sym<Range>(0, -30*xx - 60*xy + 36*x[0] + 12*x[1] - 6, 60*xy - 12*x[1]);
          out[15] = sym<Range>(60*xy - 12*x[0], -150*xy + 30*x[0] - 135*yy + 150*x[1] - 24, 30*yy - 12*x[1]);
          out[16] = sym<Range>(-60*xy + 12*x[0], 60*xy - 12*x[0] + 30*yy - 36*x[1] + 6, 0);
          out[17] = sym<Range>(0, -30*xy + 6*x[0] - 45*yy + 42*x[1] - 6, 30*yy - 12*x[1]);
        }
      }

      /** \brief Evaluate `div(div(phi))` of all shape functions `phi` at a given point
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Second derivative of all shape functions at that point
       */
      void evaluateDivDiv(const typename Traits::DomainType& x,
                          std::vector<typename Traits::DivDivType>& out) const
      {
        out.resize(size());
        if constexpr(k < 2)
          std::fill(out.begin(), out.end(), R(0));
        else if constexpr(k == 2) {
          out[0] = 0; out[1] = 0; out[2] = 30;
          out[3] = 0; out[4] = 0; out[5] = 30;
          out[6] = 0; out[7] = 0; out[8] = 30;
          out[9]  =  120; out[10] = 0;   out[11] = 0;
          out[12] = -240; out[13] = 0;   out[14] = -120;
          out[15] = -240; out[16] = 120; out[17] = 0;
        }
      }

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
    class HellanHermannJohnsonLocalBasis
    {
      using Element = E;
      using Geometry = typename E::Geometry;
      static constexpr int dim = Geometry::mydimension;
      using D = typename Geometry::ctype;

      using ReferenceLocalBasis = HellanHermannJohnsonReferenceLocalBasis<D, R, dim, k>;

      public:
        struct Traits {
          static constexpr int dimDomain = dim;
          static constexpr int dimRange = Geometry::coorddimension;
          using DomainFieldType = D;
          using DomainType = FieldVector<D,dimDomain>;
          using RangeFieldType = R;
          using RangeType = DenseTensor<R,dimRange,dimRange>;
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

          auto J = geometry_->jacobian(x);
          auto Jt = transpose(J);
          auto dx = geometry_->integrationElement(x);

          for (std::size_t i = 0; i < inValues.size(); ++i)
          {
            outValues[i] = tensordot<1>(tensordot<1>(J, inValues[i]), Jt);
            outValues[i] /= Dune::power(dx,2);
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
    class HellanHermannJohnsonLocalInterpolation
    {
      using Element = E;
      using Geometry = typename E::Geometry;
      static constexpr int dim = Geometry::mydimension;
      using D = typename Geometry::ctype;

      using size_type = std::size_t;

      static constexpr unsigned int size()
      {
        return HellanHermannJohnsonLocalCoefficients<dim,k>::size();
      }

    public:

      HellanHermannJohnsonLocalInterpolation()
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
          auto Ji = geometry_.jacobianInverse(xi);
          auto Jit = transpose(Ji);
          auto dx = geometry_.integrationElement(xi);
          return tensordot<1>(tensordot<1>(Ji, f_(xi)), Jit) * (dx*dx);
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
        out.resize(size());
        auto refElem = referenceElement(*geometry_);

        auto local_f = LocalValuedFunction{f, *geometry_};
        auto const& edgeQuadRule = Dune::QuadratureRules<D,dim-1>::rule(refElem.type(0,1), 4);
        auto const& cellQuadRule = Dune::QuadratureRules<D,dim>::rule(refElem.type(), 4);

        // 1. integral moments over the edges
        if constexpr (k == 0) {
          for (int i = 0; i < refElem.size(1); ++i) {
            auto n = refElem.integrationOuterNormal(i); n/= n.two_norm();
            auto geoInCell = refElem.template geometry<1>(i);

            out[i] = C(0);
            D edgeLength = 0;
            for (auto const& [x,w] : edgeQuadRule) {
              auto dx = geoInCell.integrationElement(x) * w;
              edgeLength += dx;

              out[i] += tensordot<1>(tensordot<1>(n,local_f(geoInCell.global(x))),n);
            }
            out[i] *= edgeLength;
          }
        }
        else if constexpr (k == 1) {
          for (int i = 0; i < refElem.size(1); ++i) {
            auto n = refElem.integrationOuterNormal(i); n/= n.two_norm();
            auto geoInCell = refElem.template geometry<1>(i);

            out[2*i] = C(0);
            out[2*i+1] = C(0);
            D edgeLength = 0;
            for (auto const& [x,w] : edgeQuadRule) {
              auto dx = geoInCell.integrationElement(x) * w;
              edgeLength += dx;

              auto nVn = tensordot<1>(tensordot<1>(n,local_f(geoInCell.global(x))),n);
              out[2*i] += (1-x) * nVn;
              out[2*i+1] += x * nVn;
            }
            out[2*i] *= edgeLength;
            out[2*i+1] *= edgeLength;
          }

          std::array B{
            DenseTensor<D,2,2>({{0,1},{1,0}}),
            DenseTensor<D,2,2>({{-2,1},{1,0}}),
            DenseTensor<D,2,2>({{0,-1},{-1,2}})
          };

          for (int i = 0; i < B.size(); ++i) {
            auto geoInCell = refElem.template geometry<0>(0);
            out[6+i] = C(0);
            for (auto const& [x,w] : cellQuadRule) {
              auto dx = geoInCell.integrationElement(x) * w;
              out[6+i] += local_f(geoInCell.global(x)).inner(B[i]);
            }
          }
        }
        else if constexpr (k == 2) {
          for (int i = 0; i < refElem.size(1); ++i) {
            auto n = refElem.integrationOuterNormal(i); n/= n.two_norm();
            auto geoInCell = refElem.template geometry<1>(i);

            out[3*i] = C(0);
            out[3*i+1] = C(0);
            out[3*i+2] = C(0);
            D edgeLength = 0;
            for (auto const& [x,w] : edgeQuadRule) {
              auto dx = geoInCell.integrationElement(x) * w;
              edgeLength += dx;

              auto nVn = tensordot<1>(tensordot<1>(n,local_f(geoInCell.global(x))),n);
              out[3*i] += (2*x*x-3*x+1) * nVn;
              out[3*i+1] += (x*(2*x-1)) * nVn;
              out[3*i+2] += (4*x*(1-x)) * nVn;
            }
            out[3*i] *= edgeLength;
            out[3*i+1] *= edgeLength;
            out[3*i+2] *= edgeLength;
          }

          auto geoInCell = refElem.template geometry<0>(0);
          for (int i = 0; i < 9; ++i)
            out[9+i] = C(0);

          for (auto const& [x,w] : cellQuadRule) {
            auto dx = geoInCell.integrationElement(x) * w;
            auto V = local_f(geoInCell.global(x));

            using T = DenseTensor<D,2,2>;
            out[9]  += V.inner(T({{0,-x[0]-x[1]+1},{-x[0]-x[1]+1,0}}));
            out[10] += V.inner(T({{2*x[0]+2*x[1]-2,-x[0]-x[1]+1},{-x[0]-x[1]+1,0}}));
            out[11] += V.inner(T({{0,-x[0]-x[1]+1},{-x[0]-x[1]+1,-2*x[0]-2*x[1]+2}}));

            out[12] += V.inner(T({{0,x[0]},{x[0],0}}));
            out[13] += V.inner(T({{-2*x[0],x[0]},{x[0],0}}));
            out[14] += V.inner(T({{0,x[0]},{x[0],2*x[0]}}));

            out[15] += V.inner(T({{0,x[1]},{x[1],0}}));
            out[16] += V.inner(T({{-2*x[1],x[1]},{x[1],0}}));
            out[17] += V.inner(T({{0,x[1]},{x[1],2*x[1]}}));
          }
        }
      }

    private:
      std::optional<Geometry> geometry_ = std::nullopt;
    };


    /** \brief Hellan-Hermann-Johnson finite element for simplices, as defined on the reference Element.
     * For more Details, see <dune/functions/functionspacebases/hellanhermannjohnsonbasis.hh>.
     *
     * \tparam D Type used for domain coordinates
     * \tparam R Type used for function values
     * \tparam dim dimension of the reference element
     */
    template<class E, class R, unsigned int k>
    class HellanHermannJohnsonLocalFiniteElement
    {
      using Element = E;
      using Geometry = typename E::Geometry;
      using D = typename Geometry::ctype;
      static constexpr int dim = Geometry::mydimension;

    public:
      HellanHermannJohnsonLocalFiniteElement()
      {
        static_assert(dim==2, "HellanHermannJohnsonLocalFiniteElement only implemented for dim=2");
      }

      /** \brief Export number types, dimensions, etc.
       */
      using size_type = std::size_t;
      using Traits = LocalFiniteElementTraits<
          HellanHermannJohnsonLocalBasis<E, R, k>,
          HellanHermannJohnsonLocalCoefficients<dim, k>,
          HellanHermannJohnsonLocalInterpolation<E, k>>;


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
        return HellanHermannJohnsonLocalCoefficients<dim,k>::size();
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
  //   HellanHermannJohnsonPreBasis
  //   HellanHermannJohnsonNode
  //
  // The pre-basis allows to create the others and is the owner of possible shared
  // state. These components do _not_ depend on the global basis and local view
  // and can be used without a global basis.
  // *****************************************************************************

  template<class GV, unsigned int k, class R>
  class HellanHermannJohnsonNode
    : public LeafBasisNode
  {
  public:
    using size_type = std::size_t;
    using Element = typename GV::template Codim<0>::Entity;
    using FiniteElement = Impl::HellanHermannJohnsonLocalFiniteElement<Element, R, k>;

    HellanHermannJohnsonNode()
      : element_(nullptr)
    {
      this->setSize(finiteElement_.size());
    }

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
    }

    //! The order of the local basis.
    unsigned int order() const { return k; }

  protected:
    FiniteElement finiteElement_;
    Element const* element_;
  };


  /**
   * \brief A pre-basis for a Hellan-Hermann-Johnson
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \tparam GV  The grid view that the FE basis is defined on
   * \tparam R   Range type used for shape function values
   * \note This only works for simplex grids
   */
  template<class GV, unsigned int k, class R>
  class HellanHermannJohnsonPreBasis
    : public LeafPreBasisMapperMixin<GV, Impl::EdgeTwist<typename GV::IndexSet>>
  {
    using Twist = Impl::EdgeTwist<typename GV::IndexSet>;
    using Base = LeafPreBasisMapperMixin<GV, Twist>;
    using Element = typename GV::template Codim<0>::Entity;
    using D = typename GV::ctype;
    static const std::size_t dim = GV::dimension;

    // helper methods to assign each subentity the number of dofs. Used by the LeafPreBasisMapperMixin.
    static constexpr std::size_t hellanHermannJohnsonLayout(Dune::GeometryType type, int gridDim)
    {
      if constexpr(k == 0)
        return type.isLine() ? 1 : 0;
      else if constexpr(k == 1)
        return type.isLine() ? 2 : int(type.dim()) == gridDim ? 3 : 0;
      else if constexpr(k == 2)
        return type.isLine() ? 3 : int(type.dim()) == gridDim ? 9 : 0;
      else
        return 0;
    }

  public:
    //! The grid view that the FE basis is defined on
    using GridView = GV;

    //! Type used for indices and size information
    using size_type = std::size_t;

    //! Template mapping root tree path to type of created tree node
    using Node = HellanHermannJohnsonNode<GridView, k, R>;

  public:

    //! Constructor for a given grid view object
    HellanHermannJohnsonPreBasis(const GV &gv)
      : Base(gv, hellanHermannJohnsonLayout, Twist{gv.indexSet(), hellanHermannJohnsonLayout(GeometryTypes::line,dim)})
    {
      static_assert(dim==2, "HellanHermannJohnsonPreBasis only implemented for dim=2");
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
     * \brief construct a PreBasisFactory for the Hellan-Hermann-Johnson Finite Element
     *
     * \tparam k  The polynomial order of the basis
     * \tparam R  Type of the basis range
     * \return a factory function to create the HellanHermannJohnsonPreBasis
     * \relates HellanHermannJohnsonBasis
     */
    template<unsigned int k, class R = double>
    auto hhj()
    {
      return []<class GV>(GV const &gridView) {
        return HellanHermannJohnsonPreBasis<GV, k, R>(gridView);
      };
    }

    /**
     * \brief construct a PreBasisFactory for the Hellan-Hermann-Johnson Finite Element
     *
     * \tparam k  The polynomial order of the basis
     * \tparam R  Type of the basis range
     * \return a factory function to create the HellanHermannJohnsonPreBasis
     * \relates HellanHermannJohnsonBasis
     */
    template<unsigned int k, class R = double>
    auto hellanHermannJohnson()
    {
      return []<class GV>(GV const &gridView) {
        return HellanHermannJohnsonPreBasis<GV, k, R>(gridView);
      };
    }

  } // end namespace BasisFactory
} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HELLANHERMANNJOHNSONBASIS_HH
