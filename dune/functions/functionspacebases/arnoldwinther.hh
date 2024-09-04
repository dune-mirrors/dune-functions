// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_ARNOLDWINTHER_HH
#define DUNE_LOCALFUNCTIONS_ARNOLDWINTHER_HH

#include <numeric>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune{
  namespace Impl{
    /**
     * \brief Implementation of the conformal Arnold-Winther Local Basis
     * \tparam D Type to represent the field in the domain
     * \tparam R Type to represent the field in the range
     */
    template<class D, class R>
    class ArnoldWintherLocalBasis
    {
      // only implemented for tringles
    public:
      static constexpr unsigned int  dim = 2;
      static constexpr unsigned int coeffSize = 24;
      using Traits = LocalBasisTraits<D, dim, FieldVector<D, dim>, R, dim* dim, FieldMatrix<R, dim, dim>, FieldMatrix<FieldVector<R, dim>, dim, dim>>;

      static constexpr unsigned int size() { return coeffSize; }

      static constexpr unsigned int order() { return 3; }

      /** \brief Evaluate all shape functions at a given point
      *
      *\param[in]  in  The evaluation point
      * \param[out] out Values of all shape functions at that point
      */
      void evaluateFunction(const typename Traits::DomainType& in, std::vector<typename Traits::RangeType>& out) const;

      /** \brief Evaluate Jacobians of all shape functions at a given point
      *
      *\param[in]  in  The evaluation point
      * \param[out] out Jacobians of all shape functions at that point
      */
      void evaluateJacobian(const typename Traits::DomainType& in, std::vector<typename Traits::JacobianType>& out) const;

      /** \brief Evaluate partial derivatives of all shape functions at a given point
      *
      * \param[in] order The partial derivative to be computed, as a multi-index
      * \param[in] in  The evaluation point
      * \param[out] out Jacobians of all shape functions at that point
      */
      void partial(const std::array<unsigned int, dim>& order, const typename Traits::DomainType& in, std::vector<typename Traits::RangeType>& out) const;

    };

    /** \brief Associations of the Arnold-Winther degrees of freedom to subentities of the
     * reference simplex
     */
    class ArnoldWintherLocalCoefficients
    {
      static constexpr unsigned int dim = 2;
    public:
      using size_type = unsigned int;

      ArnoldWintherLocalCoefficients() : localKeys_(size())
      {
        // vertices: 3 DOFs per vertex
        for (size_type i = 0; i < dim + 1; ++i)
        {
          for (size_type j = 0; j < 3; ++j)
            localKeys_[i * 3 + j] = LocalKey{i,dim, j};
        }
        // edges: 4 DOFs per edge
        for (size_type i = 0; i < dim + 1; ++i)
        {
          for (size_type j = 0; j < 4; ++j)
            localKeys_[9 + i * 4 + j] = LocalKey{i,dim - 1, j};
        }

        // element: 3 DOFs
        for (size_type i = 0; i < 3; ++i)
          localKeys_[21 + i] = LocalKey{0,0, i};

      }

      //! number of coefficients
      static constexpr size_type size() { return 24; }
      //! get i'th index
      const LocalKey& localKey(std::size_t i) const { return localKeys_[i]; }

    private:
      std::vector<LocalKey> localKeys_;
    };


    template<class D, class R>
    class ArnoldWintherLocalInterpolation
    {
      using LocalBasis = ArnoldWintherLocalBasis<D, R>;
      using size_type = std::size_t;
      using LocalCoordinate = typename LocalBasis::Traits::DomainType;
      using c_type = typename LocalBasis::Traits::DomainFieldType;
      static constexpr size_type dim = LocalBasis::Traits::dimDomain;

      template<class F>
      using ReturnType = typename std::decay_t<std::remove_cv_t<decltype(std::declval<F>()(std::declval<LocalCoordinate>()))>>;


      template< unsigned int order, class F, class Geometry>
      auto integralMoment(F const& f, Geometry const& geo)const
      {
        auto quad = QuadratureRules<c_type, Geometry::mydimension>::rule(geo.type(), LocalBasis::order() + order);

        ReturnType<F> sum = 0.;

        for (auto const& qp : quad)
        { // maybe this is an overoptimization
          auto QP = geo.global(qp.position());
          if constexpr (order == 0u)
            sum += qp.weight() * f(QP) * geo.integrationElement(qp.position());
          else if constexpr (order == 1u)
          {
            static_assert(Geometry::mydimension == 1, "First order moment only implemented for 1 dimensional facets");
            sum += qp.weight() * c_type(qp.position()) * f(QP) * geo.integrationElement(qp.position());
          }
          else
            DUNE_THROW(NotImplemented, "Higher Order Moments are not implemented");

        }
        return sum;
      }
    public:
      /** \brief Evaluate a given function at the Lagrange nodes
       *
       * \tparam F Type of function to evaluate
       * \tparam C Type used for the values of the function
       * \param[in] ff Function to evaluate
       * \param[out] out Array of function values
       */
      template <typename F, typename C>
      void interpolate(const F& ff, std::vector<C>& out) const
      {
        auto&& f = Impl::makeFunctionWithCallOperator<LocalCoordinate>(ff);

        out.resize(LocalBasis::size());
        auto refElement = Dune::ReferenceElements<double, dim>::simplex();
        auto it = out.begin();

        // point evaluations
        // 9 DOFs in total
        for (auto i = 0; i < refElement.size(dim); ++i)
        {
          auto value = f(refElement.position(i, dim));
          it[0] = value[0][0];
          it[1] = value[0][1];
          it[2] = value[1][1];
          it += 3;
        }

        // integral moment over edges
        // 12 DOFs in total
        for (auto i = 0; i < refElement.size(1); ++i)
        {
          auto average = integralMoment<0>(f, refElement.template geometry<1>(i));
          auto firstMoment = integralMoment<1>(f, refElement.template geometry<1>(i));

          std::size_t lower = (i == 2) ? 1 : 0;
          std::size_t upper = (i == 0) ? 1 : 2;
          auto tangent =
              refElement.position(upper, 2) - refElement.position(lower, 2);
          tangent /= tangent.two_norm();
          std::decay_t<decltype(tangent)> normal = {
              tangent[1], -tangent[0]}; // rotation by -90 degree


          using fRange = ReturnType<typename std::decay_t<decltype(f)>>;
          using protomotedType = typename PromotionTraits<typename FieldTraits<fRange >::field_type, c_type>::PromotedType;

          FieldVector<protomotedType, 2> normalTimesAverage, normalTimesFirstMoment;
          average.mtv(normal, normalTimesAverage);
          firstMoment.mv(normal, normalTimesFirstMoment);
          it[0] = normalTimesAverage * normal;
          it[1] = normalTimesAverage * tangent;
          it[2] = normalTimesFirstMoment * normal;
          it[3] = normalTimesFirstMoment * tangent;
          it += 4;
        }

        // integral moment on element
        // three DOFs in total
        auto average = integralMoment<0>(f, refElement.template geometry<0>(0));
        it[0] = average[0][0];
        it[1] = average[0][1];
        it[2] = average[1][1];


      }
    };

  } // namespace Impl
  /** \brief ArnoldWinther finite element for simplices
  *
  * \tparam D Type used for domain coordinates
  * \tparam R Type used for function values
  */
  template <class D, class R>
  class ArnoldWintherLocalFiniteElement
  {

  public:
    /** \brief Export number types, dimensions, etc.
   */

    using Traits =
      LocalFiniteElementTraits<Impl::ArnoldWintherLocalBasis<D, R>,
      Impl::ArnoldWintherLocalCoefficients,
      Impl::ArnoldWintherLocalInterpolation<D, R>>;

    /** \brief Default constructor
     * TODO remove?
     * \deprecated This explicit implementation only exists to work around a bug
     * in clang 3.8 which disappeared in clang 6
     */
    ArnoldWintherLocalFiniteElement() {}

    /** \brief Returns the local basis, i.e., the set of shape functions
     */
    const typename Traits::LocalBasisType& localBasis() const { return basis_; }

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

    /** \brief The number of shape functions */
    static constexpr std::size_t size()
    {
      return 24;
    }

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type() { return GeometryTypes::simplex(2); }

  private:
    Impl::ArnoldWintherLocalBasis<D, R> basis_;
    Impl::ArnoldWintherLocalCoefficients coefficients_;
    Impl::ArnoldWintherLocalInterpolation<D, R> interpolation_;
  };
} //namespace dune
#include "arnoldwinther.inc.hh"
#endif
