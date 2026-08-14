// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ARNOLDWINTHERBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ARNOLDWINTHERBASIS_HH

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/boundschecking.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/promotiontraits.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/transformedfiniteelementmixin.hh>

/**
 * \file arnoldwintherbasis.hh
 * \brief Conforming Arnold-Winther basis for planar affine triangle grids
 *
 * The Arnold-Winther element discretizes symmetric stresses in
 * \f$H(\operatorname{div})\f$ for two-dimensional elasticity.  This
 * implementation provides the lowest-order conforming element described by
 * Arnold and Winther (2002), using the transformations from Aznaran, Farrell,
 * and Kirby (2021).
 *
 * Only affine triangle grids with `dimension == dimensionworld == 2` are
 * supported. Embedded surface grids require transport between neighboring
 * tangent spaces. Non-affine maps additionally require derivatives of the
 * geometry Jacobian and nonconstant degree-of-freedom transformations, which
 * the generic dune-grid geometry interface does not provide.
 */
namespace Dune::Functions
{
  // forward declaration
  template <class GV, class R> class ArnoldWintherPreBasis;

  /** \brief Global basis for the conforming Arnold-Winther finite element
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \tparam GV Grid view type; must describe an affine planar triangle grid
   * \tparam R Range field type
   */
  template <class GV, class R = double>
  using ArnoldWintherBasis = DefaultGlobalBasis<ArnoldWintherPreBasis<GV, R>>;

  namespace Impl
  {
    /** \brief Edge orientations derived from globally unique vertex identifiers */
    using ArnoldWintherFaceOrientations = Experimental::FaceOrientations<2>;

    /** \brief Construct the edge orientations used by local transformations
     *
     * \param element Grid element whose edge orientations are requested
     * \param idSet Globally unique grid identifier set
     */
    template <class Element, class IdSet>
    ArnoldWintherFaceOrientations arnoldWintherFaceOrientations(const Element& element,
                                                                const IdSet& idSet)
    {
      constexpr int dim = 2;
      const auto& referenceElement = Dune::referenceElement<double, dim>(element.type());
      auto vertexIds = Dune::transformedRangeView(
          referenceElement.subEntities(0, 0, dim),
          [&](auto localVertexIndex) { return idSet.subId(element, localVertexIndex, dim); });
      using namespace Dune::Indices;
      return ArnoldWintherFaceOrientations(element.type(), vertexIds, _1);
    }

    /** \brief Matrix, vector, and third-order tensor types used by Arnold-Winther */
    template <class R, int dim, int dimDomain = dim>
    struct ArnoldWintherTensorTypes
    {
      /** \brief Scalar field type */
      using Scalar = R;

      /** \brief Vector type */
      using Vector = FieldVector<R, dim>;

      /** \brief Matrix type */
      using Matrix = FieldMatrix<R, dim, dim>;

      /** \brief Derivative tensor indexed by domain direction, row, and column */
      using ThreeTensor = std::array<std::array<std::array<R, dim>, dim>, dimDomain>;
    };

    /**
     * \brief Reference local basis of the conforming Arnold-Winther element
     *
     * \tparam D Type to represent the field in the domain
     * \tparam R Type to represent the field in the range
     * \tparam dim Domain dimension; only dimension two is supported
     * \tparam k Symfem element order; only order two is supported
     */
    template <class D, class R, int dim = 2, unsigned int k = 2>
    class ArnoldWintherReferenceLocalBasis
    {
      using Range = ArnoldWintherTensorTypes<R, dim>::Matrix;

    public:
      static_assert(dim == 2, "Arnold-Winther is only implemented in 2D");
      static_assert(k == 2, "Only the lowest-order Arnold-Winther element is implemented");

      /** \brief Number of local basis functions */
      static constexpr unsigned int coeffSize = 24;

      /** \brief Traits describing domains, ranges, and derivatives */
      struct Traits
      {
        /** \brief Domain field type */
        using DomainFieldType = D;

        /** \brief Domain dimension */
        constexpr static int dimDomain = dim;

        /** \brief Domain coordinate type */
        using DomainType = FieldVector<D, dim>;

        /** \brief Range field type */
        using RangeFieldType = R;

        /** \brief Number of entries in the matrix-valued range type */
        constexpr static int dimRange = dim * dim;

        /** \brief Symmetric matrix-valued range type */
        using RangeType = ArnoldWintherTensorTypes<R, dim>::Matrix;

        /** \brief Derivative tensor type */
        using JacobianType = ArnoldWintherTensorTypes<R, dim, dim>::ThreeTensor;

        /** \brief Row-wise divergence type */
        using DivergenceType = ArnoldWintherTensorTypes<R, dim>::Vector;
      };

      /** \brief Number of local basis functions */
      static constexpr unsigned int size()
      {
        return coeffSize;
      }

      /** \brief Polynomial degree of the local basis */
      static constexpr unsigned int order()
      {
        return 3;
      }

      /** \brief Evaluate all shape functions at a given point
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Values of all shape functions at that point
       */
      void evaluateFunction(const typename Traits::DomainType& in,
                            std::vector<typename Traits::RangeType>& out) const;

      /** \brief Evaluate row-wise divergences of all shape functions
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Divergences of all shape functions at that point
       */
      void evaluateDivergence(const typename Traits::DomainType& in,
                              std::vector<typename Traits::DivergenceType>& out) const;

      /** \brief Evaluate Jacobian of all shape functions - not implemented.
       *
       * \param[in]  in  The evaluation point
       * \param[out] out Jacobian of all shape functions at that point
       */
      void evaluateJacobian(const typename Traits::DomainType& in,
                            std::vector<typename Traits::JacobianType>& out) const
      {
        DUNE_THROW(Dune::NotImplemented,
          "ArnoldWintherBasis does to implement the evaluateJacobian method.");
      }

    private:
      template <class RangeType> static RangeType sym(R a00, R a01, R a11)
      {
        return RangeType({{a00, a01}, {a01, a11}});
      }
    };

    /** \brief Associations of the Arnold-Winther degrees of freedom to subentities
     * of the reference simplex
     */
    class ArnoldWintherLocalCoefficients
    {
      static constexpr unsigned int dim = 2;

    public:
      /** \brief Index type */
      using size_type = unsigned int;

      /** \brief Construct the fixed local-key layout */
      ArnoldWintherLocalCoefficients() : localKeys_(size())
      {
        // vertices: 3 DOFs per vertex
        for (size_type i = 0; i < dim + 1; ++i) {
          for (size_type j = 0; j < 3; ++j)
            localKeys_[i * 3 + j] = LocalKey{i, dim, j};
        }

        // edges: 4 DOFs per edge
        for (size_type i = 0; i < dim + 1; ++i) {
          for (size_type j = 0; j < 4; ++j)
            localKeys_[9 + i * 4 + j] = LocalKey{i, dim - 1, j};
        }

        // element: 3 DOFs
        for (size_type i = 0; i < 3; ++i)
          localKeys_[21 + i] = LocalKey{0, 0, i};
      }

      /** \brief Number of coefficients */
      static constexpr size_type size()
      {
        return 24;
      }

      /** \brief Local key of the i-th coefficient */
      const LocalKey& localKey(std::size_t i) const
      {
        return localKeys_[i];
      }

    private:
      std::vector<LocalKey> localKeys_;
    };

    /** \brief Transforms shape function values and derivatives from reference
     * element coordinates to world coordinates using the double contravariant Piola
     * transform
     *
     * See, for example:
     *  Aznaran, Francis & Kirby, Robert & Farrell, Patrick. (2021). Transformations
     *  for Piola-mapped elements.
     */
    struct DoubleContravariantPiolaTransformator
    {
    private:
      template <class ReferenceMatrix, class WorldMatrix, class Jacobian,
                class IntegrationElement>
      static void applyToMatrix(const ReferenceMatrix& referenceValue, WorldMatrix& value,
                                Jacobian const& jacobian, IntegrationElement integrationElement)
      {
        for (std::size_t k = 0; k < jacobian.N(); ++k) {
          for (std::size_t l = 0; l < jacobian.N(); ++l) {
            value[k][l] = 0;
            for (std::size_t i = 0; i < jacobian.M(); ++i)
              for (std::size_t j = 0; j < jacobian.M(); ++j)
                value[k][l] += jacobian[k][i] * referenceValue[i][j] * jacobian[l][j];
            value[k][l] /= integrationElement * integrationElement;
          }
        }
      }

    public:
      /** \brief Double Piola-transform shape-function values to physical tensors
       *
       * \param referenceValues Values on the reference element
       * \param values Transformed values on the physical element
       * \param xi Local evaluation coordinate
       * \param geometry Element geometry
       */
      template <class ReferenceValues, class WorldValues, class LocalCoordinate,
                class Geometry>
      static void applyValues(const ReferenceValues& referenceValues, WorldValues& values,
                              const LocalCoordinate& xi, const Geometry& geometry)
      {
        auto jacobian = geometry.jacobian(xi);
        auto integrationElement = geometry.integrationElement(xi);
        for (std::size_t i = 0; i < referenceValues.size(); ++i)
          applyToMatrix(referenceValues[i], values[i], jacobian, integrationElement);
      }

      /** \brief Transform derivatives of matrix-valued shape functions
       *
       * For the affine geometries supported here, the double Piola matrix is
       * constant, so it is applied independently to every reference derivative.
       *
       * \param referenceJacobians Derivatives on the reference element
       * \param jacobians Transformed derivatives on the physical element
       * \param xi Local evaluation coordinate
       * \param geometry Element geometry
       */
      template <class ReferenceJacobians, class WorldJacobians, class LocalCoordinate,
                class Geometry>
      static void applyJacobians(const ReferenceJacobians& referenceJacobians,
                                 WorldJacobians& jacobians, const LocalCoordinate& xi,
                                 const Geometry& geometry)
      {
        auto jacobian = geometry.jacobian(xi);
        auto integrationElement = geometry.integrationElement(xi);
        for (std::size_t i = 0; i < referenceJacobians.size(); ++i)
          for (std::size_t derivative = 0; derivative < referenceJacobians[i].size(); ++derivative)
            applyToMatrix(referenceJacobians[i][derivative], jacobians[i][derivative], jacobian,
                          integrationElement);
      }

      /** \brief Piola-transform affine reference divergences to world vectors
       *
       * \param[in] referenceDivergences Reference divergence values
       * \param[out] divergences World-dimensional divergence values
       * \param xi Local evaluation coordinate
       * \param geometry Element geometry
       * This uses \f$ div\,\tau = g^{-2}J\widehat{div}\,\hat\tau\f$ with
       * integration element \f$g\f$.  ArnoldWintherLocalFiniteElement::bind()
       * enforces affine geometry before this method can be reached.
       */
      template <class ReferenceDivergences, class WorldDivergences, class LocalCoordinate,
                class Geometry>
      static void applyDivergences(const ReferenceDivergences& referenceDivergences,
                                   WorldDivergences& divergences, const LocalCoordinate& xi,
                                   const Geometry& geometry)
      {
        auto jacobian = geometry.jacobian(xi);
        auto integrationElement2 =
            geometry.integrationElement(xi) * geometry.integrationElement(xi);
        for (std::size_t i = 0; i < referenceDivergences.size(); ++i) {
          divergences[i] = 0;
          for (std::size_t k = 0; k < jacobian.N(); ++k) {
            for (std::size_t j = 0; j < jacobian.M(); ++j)
              divergences[i][k] += jacobian[k][j] * referenceDivergences[i][j];
            divergences[i][k] /= integrationElement2;
          }
        }
      }
    };

    /** \brief Integrate a local-coordinate function over a reference subentity */
    template <class F, class Geometry>
    auto integralMoment(const F& f, const Geometry& geometry, int quadratureOrder)
    {
      using ctype = typename Geometry::ctype;
      using GlobalCoordinate = typename Geometry::GlobalCoordinate;

      const auto quadrature =
          QuadratureRules<ctype, Geometry::mydimension>::rule(geometry.type(), quadratureOrder);

      using ReturnType =
          std::remove_cvref_t<decltype(std::declval<F>()(std::declval<GlobalCoordinate>()))>;
      ReturnType sum = 0;

      for (const auto& quadraturePoint : quadrature) {
        const auto position = geometry.global(quadraturePoint.position());
        sum += quadraturePoint.weight() * f(position) *
               geometry.integrationElement(quadraturePoint.position());
      }
      return sum;
    }

    /** \brief Compute moments against a Lagrange basis on a reference subentity */
    template <class C, unsigned int lagrangeOrder, class F, class Geometry>
    auto lagrangeMoments(const F& f, const Geometry& geometry, int quadratureOrder)
    {
      using ctype = typename Geometry::ctype;
      using LocalCoordinate = typename Geometry::LocalCoordinate;
      using GlobalCoordinate = typename Geometry::GlobalCoordinate;
      using D = LocalCoordinate::field_type;
      static constexpr int dim = Geometry::mydimension;

      using EdgeBasis = Dune::Impl::LagrangeSimplexLocalBasis<D, C, dim, lagrangeOrder>;
      EdgeBasis edgeLagrangeBasis;
      thread_local std::vector<typename EdgeBasis::Traits::RangeType> edgeValues;
      static constexpr std::size_t edgeSize = EdgeBasis::size();

      const auto quadrature =
          QuadratureRules<ctype, Geometry::mydimension>::rule(geometry.type(), quadratureOrder);
      using ReturnType =
          std::remove_cvref_t<decltype(std::declval<F>()(std::declval<GlobalCoordinate>()))>;

      std::array<ReturnType, edgeSize> result{};
      for (const auto& quadraturePoint : quadrature) {
        const auto position = geometry.global(quadraturePoint.position());
        edgeLagrangeBasis.evaluateFunction(quadraturePoint.position(), edgeValues);
        const auto value = f(position) * quadraturePoint.weight() *
                           geometry.integrationElement(quadraturePoint.position());
        for (std::size_t i = 0; i < edgeSize; ++i)
          result[i] += value * edgeValues[i][0];
      }
      return result;
    }

    /** \brief Interpolation into the Arnold-Winther reference finite element
     *
     * The interpolation evaluates tensor values at vertices and integrates tensor
     * moments on edges and in the cell. It does not require derivatives of the
     * interpolated function.
     *
     * \tparam D Domain field type
     * \tparam R Range field type
     */
    template <class D, class R>
    class ArnoldWintherReferenceLocalInterpolation
    {
      using LocalBasis = ArnoldWintherReferenceLocalBasis<D, R>;
      using size_type = std::size_t;
      using LocalCoordinate = typename LocalBasis::Traits::DomainType;
      using ctype = typename LocalBasis::Traits::DomainFieldType;
      static constexpr size_type dim = LocalBasis::Traits::dimDomain;

    public:
      /** \brief Construct with a quadrature order */
      explicit ArnoldWintherReferenceLocalInterpolation(int quadratureOrder = 10)
          : quadratureOrder_(quadratureOrder)
      {
      }

      /** \brief Evaluate a given function at the Lagrange nodes
       *
       * \tparam F Type of function to evaluate
       * \tparam C Type used for the values of the function
       * \param[in] f Function to evaluate
       * \param[out] out Array of function values
       */
      template <typename F, typename C>
      void interpolate(const F& f, std::vector<C>& out) const
      {
        out.resize(LocalBasis::size());
        auto refElement = Dune::ReferenceElements<double, dim>::simplex();
        auto it = out.begin();

        // point evaluations
        // 9 DOFs in total
        for (auto i = 0; i < refElement.size(dim); ++i) {
          auto value = f(refElement.position(i, dim));
          it[0] = value[0][0];
          it[1] = value[0][1];
          it[2] = value[1][1];
          it += 3;
        }

        // integral moment over edges
        // 12 DOFs in total
        for (auto i = 0; i < refElement.size(1); ++i) {
          auto moments =
              lagrangeMoments<C, 1>(f, refElement.template geometry<1>(i), quadratureOrder_);

          const auto lower = refElement.subEntity(i, 1, 0, dim);
          const auto upper = refElement.subEntity(i, 1, 1, dim);
          auto tangent = refElement.position(upper, dim) - refElement.position(lower, dim);
          tangent /= tangent.two_norm();

          // Symfem defines the edge normal as the counterclockwise rotation of
          // its oriented tangent.  This is not always Dune's outward reference
          // normal, so refElement.integrationOuterNormal(i) cannot be used here
          // without an additional, edge-dependent sign correction.
          std::decay_t<decltype(tangent)> normal = {-tangent[1], tangent[0]};

          using fRange = std::decay_t<std::remove_cv_t<decltype(f(refElement.position(i, dim)))>>;
          using PromotedType = typename PromotionTraits<typename FieldTraits<fRange>::field_type,
                                                        ctype>::PromotedType;

          FieldVector<PromotedType, 2> tmp;
          for (auto&& val : moments) {
            val.mtv(normal, tmp);

            // Match Symfem's length-scaled edge-moment convention.
            const auto referenceEdgeLength = refElement.template geometry<1>(i).volume();
            it[0] = dot(tmp, normal) * referenceEdgeLength;
            it[1] = dot(tmp, tangent) * referenceEdgeLength;
            it += 2;
          }
        }

        // integral moment on element
        // three DOFs in total
        auto average = integralMoment(f, refElement.template geometry<0>(0), quadratureOrder_);
        it[0] = average[0][0];
        it[1] = average[0][1];
        it[2] = average[1][1];
      }

      int quadratureOrder_;
    };

    /** \brief Element-bound interpolation into the physical Arnold-Winther element
     *
     * \tparam Element Grid element type
     * \tparam R Range field type
     */
    template <class Element, class R>
    class ArnoldWintherLocalInterpolation
    {
      using size_type = std::size_t;
      using LocalCoordinate = typename Element::Geometry::LocalCoordinate;

      using ctype = typename Element::Geometry::ctype;
      static constexpr size_type dim = Element::Geometry::mydimension;
      static constexpr int size = 24; // number of dofs.

    public:
      /** \brief Construct with a quadrature order */
      explicit ArnoldWintherLocalInterpolation(int quadratureOrder = 10)
          : quadratureOrder_(quadratureOrder)
      {
      }

      /** \brief Bind the interpolation to an element and its edge orientations */
      void bind(const ArnoldWintherFaceOrientations& orientations, const Element& element)
      {
        faceOrientations_ = orientations;
        element_ = &element;
      }

      /** \brief Apply all physical Arnold-Winther degrees of freedom
       *
       * The constant edge directions used here are valid because the owning finite
       * element rejects non-affine geometries during binding.
       *
       * \tparam F Type of function to evaluate
       * \tparam C Type used for the values of the function
       * \param[in] f Function to evaluate
       * \param[out] out Array of function values
       */
      template <typename F, typename C>
      void interpolate(const F& f, std::vector<C>& out) const
      {

        out.resize(size);
        auto it = out.begin();
        auto refElement = referenceElement(*element_);

        // point evaluations
        // 9 DOFs in total
        for (auto i = 0u; i < element_->subEntities(dim); ++i) {
          auto geoInCell = refElement.template geometry<dim>(i);

          auto value = f(geoInCell.center());

          it[0] = value[0][0];
          it[1] = value[0][1];
          it[2] = value[1][1];
          it += 3;
        }

        // integral moment over edges
        // 12 DOFs in total
        static constexpr int momentOrder = 1;
        for (auto i = 0u; i < element_->subEntities(1); ++i) {
          auto edgeGeo = element_->template subEntity<1>(i).geometry();
          auto refEdgeGeo = refElement.template geometry<1>(i);

          auto moments = lagrangeMoments<C, momentOrder>(f, refEdgeGeo, quadratureOrder_);

          const auto lower = refElement.subEntity(i, 1, 0, dim);
          const auto upper = refElement.subEntity(i, 1, 1, dim);
          auto tangent = element_->template subEntity<dim>(upper).geometry().center() -
                         element_->template subEntity<dim>(lower).geometry().center();
          tangent /= tangent.two_norm();

          // Match the oriented-tangent normal convention used by Symfem for the
          // generated reference DOFs; it is not necessarily the outward normal.
          std::decay_t<decltype(tangent)> normal = {-tangent[1], tangent[0]};

          using fRange =
              typename std::decay_t<std::remove_cv_t<decltype(f(std::declval<LocalCoordinate>()))>>;
          using PromotedType = typename PromotionTraits<typename FieldTraits<fRange>::field_type,
                                                        ctype>::PromotedType;

          FieldVector<PromotedType, dim> normalTimesMoment;

          // Symfem scales each edge moment by the edge length.  Since moments
          // above are integrated with the reference-edge geometry, changing
          // variables to the physical edge gives the factor |e|^2/|e_hat|.
          const auto edgeMomentScale = edgeGeo.volume() * edgeGeo.volume() / refEdgeGeo.volume();
          for (std::size_t m = 0; m < momentOrder + 1; ++m) {
            if (faceOrientations_.faceOrientationIndex(i, 1))
              moments[momentOrder - m].mtv(normal, normalTimesMoment);
            else
              moments[m].mtv(normal, normalTimesMoment);
            it[0] = dot(normalTimesMoment, normal) * edgeMomentScale;
            it[1] = dot(normalTimesMoment, tangent) * edgeMomentScale;
            it += 2;
          }
        }

        // integral moment on element
        // three DOFs in total
        auto average = integralMoment(f, refElement.template geometry<0>(0), quadratureOrder_) *
                       element_->geometry().volume() / refElement.template geometry<0>(0).volume();
        it[0] = average[0][0];
        it[1] = average[0][1];
        it[2] = average[1][1];
      }

    private:
      int quadratureOrder_;
      ArnoldWintherFaceOrientations faceOrientations_;
      const Element* element_ = nullptr;
    };

    /**
     * \brief Arnold-Winther-specific block diagonal transformation matrix
     *
     * This models the matrix \f$P\f$ from Aznaran, Farrell, and Kirby. Its fixed
     * blocks correspond to three vertex blocks, three edge blocks, and one cell
     * block. The implementation intentionally only provides multiplication by the
     * transpose, which is the operation required for transforming basis values and
     * derivatives.
     *
     * \tparam T Scalar field type
     */
    template <class T>
    class ArnoldWintherBlockDiagonalMatrix
    {
    public:
      /** \brief Index type */
      using size_type = std::size_t;

      /** \brief Construct an uninitialized transformation */
      ArnoldWintherBlockDiagonalMatrix() = default;

      /** \brief Construct from vertex, edge, and cell transformation blocks */
      ArnoldWintherBlockDiagonalMatrix(const std::array<FieldMatrix<T, 3, 3>, 3>& pointDofs,
                                       const std::array<FieldMatrix<T, 4, 4>, 3>& edgeDofs,
                                       const FieldMatrix<T, 3, 3>& elementDofs)
          : transformPointDofs_(pointDofs), transformEdgeDofs_(edgeDofs),
            transformElementDofs_(elementDofs)
      {
      }

      /** \brief Multiply a block vector by the transposed transformation */
      template <class VectorIn, class VectorOut>
      void mtv(const VectorIn& x, VectorOut& y) const
      {
        DUNE_ASSERT_BOUNDS(static_cast<const void*>(&x) != static_cast<const void*>(&y));
        DUNE_ASSERT_BOUNDS(x.size() == 24);
        DUNE_ASSERT_BOUNDS(y.size() == 24);

        size_type index = 0;
        for (const auto& matrix : transformPointDofs_) {
          applyTransposedBlock(matrix, x, y, index);
          index += matrix.M();
        }
        for (const auto& matrix : transformEdgeDofs_) {
          applyTransposedBlock(matrix, x, y, index);
          index += matrix.M();
        }
        applyTransposedBlock(transformElementDofs_, x, y, index);
      }

    private:
      template <class Value, class Scalar>
      static void assignScaled(Value& out, Scalar factor, const Value& in)
      {
        if constexpr (requires { out = factor * in; })
          out = factor * in;
        else
          for (size_type i = 0; i < out.size(); ++i)
            assignScaled(out[i], factor, in[i]);
      }

      template <class Value, class Scalar>
      static void addScaled(Value& out, Scalar factor, const Value& in)
      {
        if constexpr (requires { out += factor * in; })
          out += factor * in;
        else
          for (size_type i = 0; i < out.size(); ++i)
            addScaled(out[i], factor, in[i]);
      }

      template <class Matrix, class VectorIn, class VectorOut>
      static void applyTransposedBlock(const Matrix& matrix, const VectorIn& x, VectorOut& y,
                                       size_type offset)
      {
        for (size_type i = 0; i < matrix.N(); ++i) {
          assignScaled(y[offset + i], matrix[0][i], x[offset]);
          for (size_type j = 1; j < matrix.M(); ++j)
            addScaled(y[offset + i], matrix[j][i], x[offset + j]);
        }
      }

      std::array<FieldMatrix<T, 3, 3>, 3> transformPointDofs_;
      std::array<FieldMatrix<T, 4, 4>, 3> transformEdgeDofs_;
      FieldMatrix<T, 3, 3> transformElementDofs_;
    };

    /** \brief Element-bound Arnold-Winther local finite element for triangles
     *
     * \tparam Element Grid element type
     * \tparam D Type used for domain coordinates
     * \tparam R Type used for function values
     */
    template <class Element, class D, class R>
    class ArnoldWintherLocalFiniteElement
        : public Impl::TransformedFiniteElementMixin<
              ArnoldWintherLocalFiniteElement<Element, D, R>,
              typename ArnoldWintherReferenceLocalBasis<D, R>::Traits>
    {
      static constexpr int dim = Element::Geometry::mydimension;
      static constexpr int dimWorld = Element::Geometry::coorddimension;
      using ReferenceTraits = typename ArnoldWintherReferenceLocalBasis<D, R>::Traits;
      using This = ArnoldWintherLocalFiniteElement<Element, D, R>;
      using Base = Impl::TransformedFiniteElementMixin<This, ReferenceTraits>;
      friend class Impl::TransformedLocalBasis<This, ReferenceTraits>;

    public:
      /** \brief Local finite-element traits */
      using Traits = LocalFiniteElementTraits<Impl::TransformedLocalBasis<This, ReferenceTraits>,
                                              Impl::ArnoldWintherLocalCoefficients,
                                              Impl::ArnoldWintherLocalInterpolation<Element, R>>;

      /** \brief Construct an unbound local finite element */
      ArnoldWintherLocalFiniteElement() : Base() {}

      /** \brief Returns the assignment of the degrees of freedom to the element
       * subentities
       */
      const typename Traits::LocalCoefficientsType& localCoefficients() const
      {
        return coefficients_;
      }

      /** \brief Return the object evaluating physical degrees of freedom */
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
      static constexpr GeometryType type()
      {
        return GeometryTypes::simplex(2);
      }

      /**
       * \brief Bind the local finite element and assemble its transformation
       *
       * \param orientations Edge orientations relative to global vertex IDs
       * \param element Element to bind to
       */
      void bind(const ArnoldWintherFaceOrientations& orientations, Element const& element)
      {
        if constexpr (dim != dimWorld)
          DUNE_THROW(Dune::NotImplemented, "Arnold-Winther requires dimension=dimensionworld=2");
        else {
          if (not element.geometry().affine())
            DUNE_THROW(Dune::NotImplemented,
                       "Arnold-Winther requires an affine geometry: the generic "
                       "dune-grid Geometry interface does not expose the Jacobian "
                       "derivatives required by the non-affine double-Piola "
                       "divergence and DOF transformations");
          faceOrientations_ = orientations;
          element_ = &element;
          interpolation_.bind(orientations, element);
          fillMatrix(element.geometry());
        }
      }

    protected:
      /** \brief Return the reference local basis used by the transformation mixin */
      Impl::ArnoldWintherReferenceLocalBasis<D, R> const& referenceLocalBasis() const
      {
        return basis_;
      }

      /** \brief Apply the basis and double-Piola transformations
       *
       * The input and output containers must provide random access. The value type
       * selects the transformation for values, derivatives, or divergences.
       */
      template <class InputValues, class OutputValues>
      void transform(InputValues const& inValues, OutputValues& outValues) const
      {
        using InputValue = typename InputValues::value_type;
        std::vector<InputValue> transformedReferenceValues(size());
        mat_.mtv(inValues, transformedReferenceValues);
        // bind() rejects non-affine geometries, so the Piola map is constant and
        // can be evaluated at any reference point.
        const auto x = Dune::referenceElement<double, 2>(GeometryTypes::simplex(2)).position(0, 0);

        if constexpr (std::is_same_v<InputValue, typename ReferenceTraits::RangeType>)
          DoubleContravariantPiolaTransformator::applyValues(transformedReferenceValues, outValues,
                                                             x, element_->geometry());
        else if constexpr (std::is_same_v<InputValue, typename ReferenceTraits::DivergenceType>)
          DoubleContravariantPiolaTransformator::applyDivergences(
              transformedReferenceValues, outValues, x, element_->geometry());
        else if constexpr (std::is_same_v<InputValue, typename ReferenceTraits::JacobianType>)
          DoubleContravariantPiolaTransformator::applyJacobians(transformedReferenceValues,
                                                                outValues, x, element_->geometry());
        else
          static_assert(Dune::AlwaysFalse<InputValue>::value,
                        "Unsupported Arnold-Winther transformed value type");
      }

    private:
      template <class Geometry>
      void fillMatrix(Geometry const& geometry)
      {
        // Per-edge blocks of P=Q^{-1}.  They convert the four physical edge
        // moments (two Lagrange moments, each with nn and nt components) to the
        // corresponding reference functionals, including edge reorientation.
        std::array<Dune::FieldMatrix<R, 4, 4>, 3> W_k;
        // Symmetric-tensor block of P=Q^{-1}.  Before inversion W represents
        // vec_s(J tau J^T).  The integration element g specializes P to the
        // vertex blocks (g^2 W^{-1}) and cell block (g W^{-1}).
        Dune::FieldMatrix<R, 3, 3> W;
        // The geometry is affine, so any point in the reference triangle is valid.
        auto x = Dune::referenceElement<double, 2>(GeometryTypes::simplex(2)).position(0, 0);
        const auto jacobian = geometry.jacobian(x);
        const auto integrationElement = geometry.integrationElement(x);
        const auto jacobianDeterminant = jacobian.determinant();

        // By default, edges point from the vertex with the smaller index
        // to the vertex with the larger index. Note that the alpha and beta are
        // invariant of orientation, since the normals/tangents appear twice in
        // their definitions.

        // Get local tangents and construct the edge transformation blocks.
        auto refElement = Dune::referenceElement<double, 2>(geometry.type());
        for (std::size_t i = 0; i < 3; ++i) {
          const auto lower = refElement.subEntity(i, 1, 0, 2);
          const auto upper = refElement.subEntity(i, 1, 1, 2);
          auto tangent = refElement.position(upper, 2) - refElement.position(lower, 2);

          tangent /= tangent.two_norm();

          Dune::FieldMatrix<R, 2, 2> referenceG = {{-tangent[1], tangent[0]},
                                                   {tangent[0], tangent[1]}};
          auto tmp = tangent, tmp2 = tangent;
          jacobian.mv(tangent, tmp);
          jacobian.mtv(tmp, tmp2);
          referenceG.mtv(tmp2, tmp);
          // This computes referenceG^T J^T J tangent.
          const auto alpha = tmp[0] / jacobianDeterminant;
          const auto beta = tmp[1] / jacobianDeterminant;

          // Store the inverse edge push-forward block.  Symfem's edge moments are
          // scaled by their edge length, so the |e_hat|/|e| factor in the
          // unscaled push-forward cancels and no edge-length ratio appears here.
          W_k[i] = 0;
          if (faceOrientations_.faceOrientationIndex(i, 1)) {
            W_k[i][2][0] = 1.;
            W_k[i][3][0] = -alpha / beta;
            W_k[i][3][1] = 1. / beta;
            W_k[i][0][2] = 1.;
            W_k[i][1][2] = -alpha / beta;
            W_k[i][1][3] = 1. / beta;

          } else {
            W_k[i][0][0] = 1.;
            W_k[i][1][0] = -alpha / beta;
            W_k[i][1][1] = 1. / beta;
            W_k[i][2][2] = 1.;
            W_k[i][3][2] = -alpha / beta;
            W_k[i][3][3] = 1. / beta;
          }
        }
        // Assemble and invert the symmetric-tensor component block.
        W[0][0] = jacobian[0][0] * jacobian[0][0];
        W[0][1] = 2. * jacobian[0][0] * jacobian[0][1];
        W[0][2] = jacobian[0][1] * jacobian[0][1];
        W[1][0] = jacobian[0][0] * jacobian[1][0];
        W[1][1] = jacobian[0][0] * jacobian[1][1] + jacobian[0][1] * jacobian[1][0];
        W[1][2] = jacobian[0][1] * jacobian[1][1];
        W[2][0] = jacobian[1][0] * jacobian[1][0];
        W[2][1] = 2. * jacobian[1][0] * jacobian[1][1];
        W[2][2] = jacobian[1][1] * jacobian[1][1];
        W.invert();
        W *= integrationElement * integrationElement;
        mat_ = ArnoldWintherBlockDiagonalMatrix<R>{
            std::array<Dune::FieldMatrix<R, 3, 3>, 3>{W, W, W}, W_k, W / integrationElement};
      }

    private:
      Impl::ArnoldWintherReferenceLocalBasis<D, R> basis_;
      Traits::LocalCoefficientsType coefficients_;
      Traits::LocalInterpolationType interpolation_;
      // Matrix P from Aznaran, Farrell, and Kirby.
      ArnoldWintherBlockDiagonalMatrix<R> mat_;
      ArnoldWintherFaceOrientations faceOrientations_;
      const Element* element_ = nullptr;
    };

  } // namespace Impl

  // forward declaration
  template <class GV, class R> class ArnoldWintherNode;

  /** \brief Pre-basis for the conforming Arnold-Winther global basis
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \tparam GV Grid view type
   * \tparam R Range field type
   */
  template <class GV, class R>
  class ArnoldWintherPreBasis : public LeafPreBasisMapperMixin<GV>
  {
    static constexpr int dim = GV::dimension;
    static_assert(dim == 2, "ArnoldWinther PreBasis only implemented for 2d simplices");
    using Base = LeafPreBasisMapperMixin<GV>;

    // Assign the fixed number of DOFs to each subentity type.
    static constexpr auto arnoldWintherMapperLayout(Dune::GeometryType type, int)
    {
      if (type.isVertex())
        return 3;
      if (type.isLine())
        return 4;
      if (type.isTriangle())
        return 3;
      return 0;
    }

  public:
    /** \brief Grid view type */
    using GridView = GV;

    /** \brief Type used for indices and sizes */
    using size_type = std::size_t;

    /** \brief Local basis-tree node type */
    using Node = ArnoldWintherNode<GridView, R>;

    /** \brief Construct for a grid view */
    explicit ArnoldWintherPreBasis(const GV& gv) : Base(gv, arnoldWintherMapperLayout)
    {
      if constexpr (GV::dimension != GV::dimensionworld)
        DUNE_THROW(Dune::NotImplemented,
          "Arnold-Winther requires dimension=dimensionworld=2.");
    }

    /** \brief Update the stored grid view after a grid change */
    void update(const GridView& gv)
    {
      Base::update(gv);
    }

    /** \brief Create an unbound local basis-tree node */
    Node makeNode() const
    {
      return Node{this->gridView()};
    }

    /** \brief Return the globally consistent edge orientations of an element */
    template <class Element>
    auto faceOrientations(const Element& element) const
    {
      return Impl::arnoldWintherFaceOrientations(element, this->gridView().grid().globalIdSet());
    }
  };

  /** \brief Leaf basis node for the conforming Arnold-Winther basis
   *
   * \tparam GV Grid view type
   * \tparam R Range field type
   */
  template <class GV, class R>
  class ArnoldWintherNode : public LeafBasisNode
  {
  public:
    /** \brief Type used for indices and sizes */
    using size_type = std::size_t;

    /** \brief Codimension-zero entity type */
    using Element = typename GV::template Codim<0>::Entity;

    /** \brief Element-bound local finite-element type */
    using FiniteElement = Impl::ArnoldWintherLocalFiniteElement<Element, typename GV::ctype, R>;

    /** \brief Construct an unbound node for a grid view */
    explicit ArnoldWintherNode(const GV& gridView) : gridView_(&gridView) {}

    /** \brief Return the bound element */
    const Element& element() const
    {
      return *element_;
    }

    /** \brief Return the LocalFiniteElement for the element we are bound to
     *
     * The LocalFiniteElement implements the corresponding interfaces of the
     * dune-localfunctions module
     */
    const FiniteElement& finiteElement() const
    {
      return finiteElement_;
    }

    /** \brief Bind the node and local finite element to an element */
    void bind(const Element& element)
    {
      if (not element.type().isSimplex())
        DUNE_THROW(Dune::NotImplemented,
          "ArnoldWintherBasis can only be bound to simplex elements");
      element_ = &element;
      finiteElement_.bind(
          Impl::arnoldWintherFaceOrientations(element, gridView_->grid().globalIdSet()), *element_);
      this->setSize(finiteElement_.size());
    }

    /** \brief Polynomial order of the local finite element */
    unsigned int order() const
    {
      return 3;
    }

  protected:
    FiniteElement finiteElement_;
    const Element* element_ = nullptr;
    const GV* gridView_;
  };

  namespace BasisFactory
  {
    /**
     * \brief Create a pre-basis factory that can create ArnoldWinther pre-basis
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     * \tparam Range Number type used for shape-function values
     */
    template <typename Range = double>
    auto arnoldWinther()
    {
      return [](auto const& gridView) {
        return ArnoldWintherPreBasis<std::decay_t<decltype(gridView)>, Range>(gridView);
      };
    }
  } // namespace BasisFactory
} // namespace Dune::Functions

#include <dune/functions/functionspacebases/arnoldwintherbasis_inc.hh>

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ARNOLDWINTHERBASIS_HH
