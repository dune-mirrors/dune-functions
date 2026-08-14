#ifndef DUNE_C1ELEMENTS_ARNOLDWINTHER_HH
#define DUNE_C1ELEMENTS_ARNOLDWINTHER_HH

#include <array>
#include <numeric>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
// #include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/transformedfiniteelementmixin.hh>

namespace Dune {
namespace Functions {

/**
* \brief Implementation of the conforming Arnold-Winther element
*   This is a finite element used to discretize the stress for two-dimensional
* elasticity, originally proposed in "Arnold, D. N., & Winther, R. (2002).
* Mixed finite elements for elasticity." As such, its shape functions take
* values in the space of symmetric 2x2 matrices, whose (rowwise) divergence is
* in L2.
* This comes with some complications in the Dune framework, in particular,
* there is no data structure for 3-tensors, which is the JacobianType of this
* finite element.  This therefore omits the evaluateJacobian method and only
* implements evaluateDivergence directly.  The transformation is based on
* "Aznaran, Francis & Kirby, Robert & Farrell, Patrick.
* (2021). Transformations for Piola-mapped elements."
*
* Only grids with dimension=dimensionworld=2 are supported.  A conforming
* Arnold-Winther space on an embedded surface would have to identify traction
* vectors and vertex tensor values in different tangent spaces.  This requires
* an explicit tangent-space transport on the non-coplanar neighboring elements
* occurring in ordinary surface meshes.
*
* Non-affine geometries are rejected.  Their double-Piola divergence contains
* derivatives of the geometry Jacobian and integration element, which the
* generic dune-grid Geometry interface does not expose.  Moreover, their DOF
* push-forward is no longer represented by the constant blocks used here.
*/
namespace Impl {
using ArnoldWintherFaceOrientations = Experimental::FaceOrientations<2>;

template <class Element, class IdSet>
ArnoldWintherFaceOrientations
arnoldWintherFaceOrientations(const Element &element, const IdSet &idSet) {
  constexpr int dim = 2;
  const auto &referenceElement = Dune::referenceElement<double, dim>(element.type());
  auto vertexIds = Dune::transformedRangeView(
      referenceElement.subEntities(0, 0, dim), [&](auto localVertexIndex) {
        return idSet.subId(element, localVertexIndex, dim);
      });
  using namespace Dune::Indices;
  return ArnoldWintherFaceOrientations(element.type(), vertexIds, _1);
}

template <class R, int dim, int dimDomain = dim>
struct ArnoldWintherTensorTypes {
  using Scalar = R;

  using Vector = FieldVector<R, dim>;

  using Matrix = FieldMatrix<R, dim, dim>;

  using ThreeTensor =
      std::array<std::array<std::array<R, dim>, dim>, dimDomain>;
};

/**
 * \brief Implementation of the conformal Arnold-Winther Local Basis
 * \tparam D Type to represent the field in the domain
 * \tparam R Type to represent the field in the range
 */
template <class D, class R, int dim = 2, unsigned int k = 2>
 class ArnoldWintherReferenceLocalBasis {
  // only implemented for triangles
  // only lowest order implemented
  using Range = ArnoldWintherTensorTypes<R, dim>::Matrix;
public:
  static_assert(dim == 2, "AW Element only implemented in 2d");
  static constexpr unsigned int coeffSize = 24;
  // using Traits = LocalBasisTraits<D, dim, typename ArnoldWinterTensorTypes<D,
  // dim>::Vector, R, dim* dim, typename ArnoldWinterTensorTypes<R,
  // dim>::Matrix, typename ArnoldWinterTensorTypes<R, dim>::Vector>;

  struct Traits {
    //! \brief Export type for domain field
    using DomainFieldType = D;

    //! \brief dimension of the domain
    constexpr static int dimDomain = dim;

    //! \brief domain type
    using DomainType = Dune::FieldVector<D, dim>;

    //! \brief Export type for range field
    using RangeFieldType = R;

    //! \brief dimension of the range
    // TODO discuss what this should be. For now we take it as the entries in
    // the range matrix
    constexpr static int dimRange = dim * dim;

    //! \brief range type
    using RangeType = ArnoldWintherTensorTypes<R, dim>::Matrix;

    /** \brief Type to represent derivative
     */
    using JacobianType =
        ArnoldWintherTensorTypes<R, dim, dim>::ThreeTensor;

    /** \brief Type to represent the rowwise divergence
     */
    using DivergenceType = ArnoldWintherTensorTypes<R, dim>::Vector;
  };

  static constexpr unsigned int size() { return coeffSize; }

  static constexpr unsigned int order() { return 3; }

  /** \brief Evaluate all shape functions at a given point
   *
   *\param[in]  in  The evaluation point
   * \param[out] out Values of all shape functions at that point
   */
  void evaluateFunction(const typename Traits::DomainType &in,
                        std::vector<typename Traits::RangeType> &out) const;

  /** \brief Evaluate Jacobians of all shape functions at a given point
   *
   *\param[in]  in  The evaluation point
   * \param[out] out Jacobians of all shape functions at that point
   */
  void evaluateJacobian(const typename Traits::DomainType &in,
                        std::vector<typename Traits::JacobianType> &out) const;

  /** \brief Evaluate Jacobians of all shape functions at a given point
   *
   *\param[in]  in  The evaluation point
   * \param[out] out Jacobians of all shape functions at that point
   */
  void
  evaluateDivergence(const typename Traits::DomainType &in,
                     std::vector<typename Traits::DivergenceType> &out) const;

  /** \brief Evaluate partial derivatives of all shape functions at a given
   * point
   *
   * \param[in] order The partial derivative to be computed, as a multi-index
   * \param[in] in  The evaluation point
   * \param[out] out Jacobians of all shape functions at that point
   */
  void partial(const std::array<unsigned int, dim> &order,
               const typename Traits::DomainType &in,
               std::vector<typename Traits::RangeType> &out) const;
  private:
    template <class Range>
    static Range sym(R a00, R a01, R a11)
    {
      return Range({{a00,a01},{a01,a11}});
    }
};

/** \brief Associations of the Arnold-Winther degrees of freedom to subentities
 * of the reference simplex
 */
class ArnoldWintherLocalCoefficients {
  static constexpr unsigned int dim = 2;

public:
  using size_type = unsigned int;

  ArnoldWintherLocalCoefficients() : localKeys_(size()) {
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

  //! number of coefficients
  static constexpr size_type size() { return 24; }
  //! get i'th index
  const LocalKey &localKey(std::size_t i) const { return localKeys_[i]; }

private:
  std::vector<LocalKey> localKeys_;
};

/** \brief Transforms shape function values and derivatives from reference
 * element coordinates to world coordinates using the double contravariant Piola
 * transform
 *
 * \TODO Maybe make this a private class, since the treatment of divergences is
 * tailored to AW
 *
 * See for example:
 *  Aznaran, Francis & Kirby, Robert & Farrell, Patrick. (2021). Transformations
 * for Piola-mapped elements.
 */
struct DoubleContravariantPiolaTransformator {
private:
  template <class ReferenceMatrix, class WorldMatrix, class Jacobian,
            class IntegrationElement>
  static void applyToMatrix(const ReferenceMatrix& referenceValue,
                            WorldMatrix& value, Jacobian const &jacobian,
                            IntegrationElement integrationElement) {
    for (std::size_t k = 0; k < jacobian.N(); ++k) {
      for (std::size_t l = 0; l < jacobian.N(); ++l) {
        value[k][l] = 0;
        for (std::size_t i = 0; i < jacobian.M(); ++i)
          for (std::size_t j = 0; j < jacobian.M(); ++j)
            value[k][l] += jacobian[k][i] * referenceValue[i][j]
                           * jacobian[l][j];
        value[k][l] /= integrationElement * integrationElement;
      }
    }
  }

public:
  /** \brief Double Piola-transform shape-function values to world tensors. */
  template <class ReferenceValues, class WorldValues, class LocalCoordinate,
            class Geometry>
  static void applyValues(const ReferenceValues& referenceValues,
                          WorldValues& values, const LocalCoordinate &xi,
                          const Geometry &geometry) {
    auto jacobian = geometry.jacobian(xi);
    auto integrationElement = geometry.integrationElement(xi);
    for (std::size_t i = 0; i < referenceValues.size(); ++i)
      applyToMatrix(referenceValues[i], values[i], jacobian,
                    integrationElement);
  }

  /** \brief Transform derivatives of matrix-valued shape functions.
   *
   * For the affine geometries supported here, the double Piola matrix is
   * constant, so it is applied independently to every reference derivative.
   */
  template <class ReferenceJacobians, class WorldJacobians,
            class LocalCoordinate, class Geometry>
  static void applyJacobians(const ReferenceJacobians& referenceJacobians,
                             WorldJacobians& jacobians,
                             const LocalCoordinate &xi,
                             const Geometry &geometry) {
    auto jacobian = geometry.jacobian(xi);
    auto integrationElement = geometry.integrationElement(xi);
    for (std::size_t i = 0; i < referenceJacobians.size(); ++i)
      for (std::size_t derivative = 0;
           derivative < referenceJacobians[i].size(); ++derivative)
        applyToMatrix(referenceJacobians[i][derivative],
                      jacobians[i][derivative], jacobian,
                      integrationElement);
  }

  /** \brief Piola-transform affine reference divergences to world vectors
   *
   * \param[in] referenceDivergences Reference divergence values
   * \param[out] divergences World-dimensional divergence values
   * This uses \f$ div\,\tau = g^{-2}J\widehat{div}\,\hat\tau\f$ with
   * integration element \f$g\f$.  ArnoldWintherLocalFiniteElement::bind()
   * enforces affine geometry before this method can be reached.
   */
  template <class ReferenceDivergences, class WorldDivergences,
            class LocalCoordinate, class Geometry>
  static void applyDivergences(const ReferenceDivergences& referenceDivergences,
                               WorldDivergences& divergences,
                               const LocalCoordinate &xi,
                               const Geometry &geometry) {
    auto jacobian = geometry.jacobian(xi);
    auto integrationElement2 =
        geometry.integrationElement(xi) * geometry.integrationElement(xi);
    for (std::size_t i = 0; i < referenceDivergences.size(); ++i) {
      divergences[i] = 0;
      for (std::size_t k = 0; k < jacobian.N(); ++k) {
        for (std::size_t j = 0; j < jacobian.M(); ++j)
          divergences[i][k] +=
              jacobian[k][j] * referenceDivergences[i][j];
        divergences[i][k] /= integrationElement2;
      }
    }
  }
};


template <unsigned int momentOrder, class F, class Geometry>
static auto integralMoment(F const &f, Geometry const &geo, int quadOrder) {
  using ctype = typename Geometry::ctype;
  using LocalCoordinate = typename Geometry::LocalCoordinate;
  using GlobalCoordinate = typename Geometry::GlobalCoordinate;

  auto quad = QuadratureRules<ctype, Geometry::mydimension>::rule(
      geo.type(), quadOrder);

  typename std::decay_t<std::remove_cv_t<decltype(std::declval<F>()(
        std::declval<GlobalCoordinate>()))>> sum = 0.;

  for (auto const &qp : quad) {
    auto QP = geo.global(qp.position());
    if constexpr (momentOrder == 0u)
      sum += qp.weight() * f(QP) * geo.integrationElement(qp.position());
    else
    {
      static_assert(
          Geometry::mydimension == 1,
          "Higher order moment only implememted for 1 dimensional facets");
      sum += qp.weight() * Dune::power(ctype(qp.position()), momentOrder) * f(QP) *
              geo.integrationElement(qp.position());
    }
  }
  return sum;
}

template <class C, unsigned int lagrangeOrder, class F, class Geometry>
static auto LagrangeMoment(F const &f, Geometry const &geo, int quadOrder) {
  using ctype = typename Geometry::ctype;
  using LocalCoordinate = typename Geometry::LocalCoordinate;
  using GlobalCoordinate = typename Geometry::GlobalCoordinate;
  using D = LocalCoordinate::field_type;
  static constexpr int dim = Geometry::mydimension;

  Dune::Impl::LagrangeSimplexLocalBasis<D,C, dim, lagrangeOrder> edgeLagrangebasis;
  thread_local std::vector< typename Dune::Impl::LagrangeSimplexLocalBasis<D,C, dim, lagrangeOrder>::Traits::RangeType> edgeValues;
  static constexpr std::size_t edgeSize = edgeLagrangebasis.size();

  auto quad = QuadratureRules<ctype, Geometry::mydimension>::rule(
      geo.type(), quadOrder);


  using ReturnType = std::remove_cvref_t<decltype(std::declval<F>()(
        std::declval<GlobalCoordinate>()))>;

  std::array<ReturnType, edgeSize> result;

  for (std::size_t i = 0; i < edgeSize; ++i)
    result[i] = 0.;

  for (auto const &qp : quad) {
    auto QP = geo.global(qp.position());
    edgeLagrangebasis.evaluateFunction(qp.position(), edgeValues);
    auto value = f(QP)*qp.weight()*geo.integrationElement(qp.position());
    for (std::size_t i = 0; i < edgeSize; ++i){
      result[i] += value*edgeValues[i][0];
    }
  }
  return result;
}

/**\brief The Arnold-Winther degrees of freedom on the reference Trianlge
 * \TODO add functionalDescriptors
 * \TODO Actually, all we need now is the global interpolation, see
 * cubichermitebasis.hh
 */
template <class D, class R>
class ArnoldWintherReferenceLocalInterpolation {
  using LocalBasis = ArnoldWintherReferenceLocalBasis<D, R>;
  using size_type = std::size_t;
  using LocalCoordinate = typename LocalBasis::Traits::DomainType;
  using c_type = typename LocalBasis::Traits::DomainFieldType;
  static constexpr size_type dim = LocalBasis::Traits::dimDomain;

public:
  ArnoldWintherReferenceLocalInterpolation(int quadOrder = 10)
  : quadratureOrder(quadOrder){}

  /** \brief Evaluate a given function at the Lagrange nodes
   *
   * \tparam F Type of function to evaluate
   * \tparam C Type used for the values of the function
   * \param[in] ff Function to evaluate
   * \param[out] out Array of function values
   */
  template <typename F, typename C>
  void interpolate(const F &f, std::vector<C> &out) const {

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
      auto moments = LagrangeMoment<C,1>(f, refElement.template geometry<1>(i), quadratureOrder);

      const auto lower = refElement.subEntity(i, 1, 0, dim);
      const auto upper = refElement.subEntity(i, 1, 1, dim);
      auto tangent =
          refElement.position(upper, dim) - refElement.position(lower, dim);
      tangent /= tangent.two_norm();
      // Symfem defines the edge normal as the counterclockwise rotation of
      // its oriented tangent.  This is not always Dune's outward reference
      // normal, so refElement.integrationOuterNormal(i) cannot be used here
      // without an additional, edge-dependent sign correction.
      std::decay_t<decltype(tangent)> normal = {
          -tangent[1], tangent[0]};

      using fRange = std::decay_t<std::remove_cv_t<decltype(f(refElement.position(i, dim)))>>;
      using protomotedType =
          typename PromotionTraits<typename FieldTraits<fRange>::field_type,
                                   c_type>::PromotedType;

      FieldVector<protomotedType, 2> tmp;
      for (auto&& val : moments){
        val.mtv(normal, tmp);

        // Match Symfem's length-scaled edge-moment convention.
        const auto referenceEdgeLength =
            refElement.template geometry<1>(i).volume();
        it[0] = dot(tmp, normal) * referenceEdgeLength;
        it[1] = dot(tmp, tangent) * referenceEdgeLength;
        it += 2;
      }

    }

    // integral moment on element
    // three DOFs in total
    auto average = integralMoment<0>(f, refElement.template geometry<0>(0), quadratureOrder);
    it[0] = average[0][0];
    it[1] = average[0][1];
    it[2] = average[1][1];
  }

  int quadratureOrder ;
};


 /** \brief The Arnold-Winther degrees of freedom on the reference Trianlge
 * \TODO add functionalDescriptors
 */
template <class Element, class R>
class ArnoldWintherLocalInterpolation {
  using size_type = std::size_t;
  using LocalCoordinate = typename Element::Geometry::LocalCoordinate;

  using ctype = typename Element::Geometry::ctype;
  static constexpr size_type dim = Element::Geometry::mydimension;
  static constexpr int size = 24; // number of dofs.
public:
  ArnoldWintherLocalInterpolation(int quadOrder = 10)
  : quadratureOrder(quadOrder)
  {}

  void bind(const ArnoldWintherFaceOrientations &orientations,
            Element const& e){
    faceOrientations_ = orientations;
    element = &e;
  }

  /** \brief Evaluate a given function at the Lagrange nodes
   * \TODO this implemenation assumes constant normals/tangents. Generalize to curved grids
   * \tparam F Type of function to evaluate
   * \tparam C Type used for the values of the function
   * \param[in] ff Function to evaluate
   * \param[out] out Array of function values
   */
  template <typename F, typename C>
  void interpolate(const F &f, std::vector<C> &out) const {

    out.resize(size);
    auto it = out.begin();
    auto refElement = referenceElement(*element);
    // point evaluations
    // 9 DOFs in total
    for (auto i = 0u; i < element->subEntities(dim); ++i) {
      auto geoInCell = refElement.template geometry<dim>(i);

      auto value = f(geoInCell.center());

      it[0] = value[0][0];
      it[1] = value[0][1];
      it[2] = value[1][1];
      it += 3;
    }

    // integral moment over edges
    // 12 DOFs in total
    static constexpr int momentOrder =1;
    for (auto i = 0u; i < element->subEntities(1); ++i) {
      auto edgeGeo = (*element).template subEntity<1>(i).geometry();
      auto refEdgeGeo = refElement.template geometry<1>(i);

      auto moments = LagrangeMoment<C, momentOrder>(f, refEdgeGeo, quadratureOrder);

      const auto lower = refElement.subEntity(i, 1, 0, dim);
      const auto upper = refElement.subEntity(i, 1, 1, dim);
      auto tangent = (*element).template subEntity<dim>(upper).geometry().center() - (*element).template subEntity<dim>(lower).geometry().center();
      tangent /= tangent.two_norm();
      // Match the oriented-tangent normal convention used by Symfem for the
      // generated reference DOFs; it is not necessarily the outward normal.
      std::decay_t<decltype(tangent)> normal = {-tangent[1], tangent[0]};

      using fRange = typename std::decay_t<std::remove_cv_t<decltype(f(std::declval<LocalCoordinate>()))>>;
      using protomotedType =
          typename PromotionTraits<typename FieldTraits<fRange>::field_type,
                                   ctype>::PromotedType;

      FieldVector<protomotedType, dim> normalTimesMoment;
      // Symfem scales each edge moment by the edge length.  Since moments
      // above are integrated with the reference-edge geometry, changing
      // variables to the physical edge gives the factor |e|^2/|e_hat|.
      const auto edgeMomentScale =
          edgeGeo.volume() * edgeGeo.volume() / refEdgeGeo.volume();
      for (std::size_t m = 0; m < momentOrder + 1; ++m){
        if (faceOrientations_.faceOrientationIndex(i, 1))
          moments[momentOrder -m].mtv(normal, normalTimesMoment);
        else
          moments[m].mtv(normal, normalTimesMoment);
        it[0] = dot(normalTimesMoment, normal) * edgeMomentScale;
        it[1] = dot(normalTimesMoment, tangent) * edgeMomentScale;
        it += 2;
      }

    }

    // integral moment on element
    // three DOFs in total
    auto average = integralMoment<0>(f, refElement.template geometry<0>(0), quadratureOrder)*(*element).geometry().volume()/refElement.template geometry<0>(0).volume();
    it[0] = average[0][0];
    it[1] = average[0][1];
    it[2] = average[1][1];
  }

private:
  int quadratureOrder;
  ArnoldWintherFaceOrientations faceOrientations_;
  Element const* element = nullptr;

};

// \TODO make this fullfill Dune interfaces
// \TODO make generic in matrix Types and sizes
// \TODO maybe make this private class
/**
 * \brief Block Diagonal Matrix with hardcoded dimensions that fit the Arnold
 * Winther FE transformation. It models the $P$ Matrix in the above mentioned
 * paper. The only operation needed for this purpose is the multiplication of
 * its transpose with a vector and access to its inverse. This Vector however,
 * has values which are Matrices or 3-tensors (for now, Matrices of Vectors).
 *
 * \tparam T FieldType
 */
template <class T>
class ArnoldWintherBlockDiagonalMatrix {
  using This = ArnoldWintherBlockDiagonalMatrix<T>;

  std::array<Dune::FieldMatrix<T, 3, 3>, 3> transformPointDofs_;
  std::array<Dune::FieldMatrix<T, 4, 4>, 3> transformEdgeDofs_;
  Dune::FieldMatrix<T, 3, 3> transformElementDofs_;

public:
  using value_type = T;
  using field_type = T;
  using size_type = std::size_t;

  ArnoldWintherBlockDiagonalMatrix() = default;
  ArnoldWintherBlockDiagonalMatrix(
      std::array<FieldMatrix<T, 3, 3>, 3> const &pointDofs,
      std::array<FieldMatrix<T, 4, 4>, 3> const &edgeDofs,
      FieldMatrix<T, 3, 3> const &elementDofs)
      : transformPointDofs_(pointDofs), transformEdgeDofs_(edgeDofs),
        transformElementDofs_(elementDofs) {}

  template <class VectorIn, class VectorOut>
  void mv(VectorIn const &x, VectorOut &y) const {
    DUNE_ASSERT_BOUNDS((void *)(&x) != (void *)(&y));
    DUNE_ASSERT_BOUNDS(x.size() == 24);
    DUNE_ASSERT_BOUNDS(y.size() == 24);

    size_type index = 0;
    for (auto const &mat : transformPointDofs_) {
      applyBlock(mat, x, y, index, false);
      index += mat.M();
    }
    for (auto const &mat : transformEdgeDofs_) {
      applyBlock(mat, x, y, index, false);
      index += mat.M();
    }
    applyBlock(transformElementDofs_, x, y, index, false);
  }

  template <class VectorIn, class VectorOut>
  void mtv(VectorIn const &x, VectorOut &y) const {
    DUNE_ASSERT_BOUNDS((void *)(&x) != (void *)(&y));
    DUNE_ASSERT_BOUNDS(x.size() == 24);
    DUNE_ASSERT_BOUNDS(y.size() == 24);

    size_type index = 0;
    for (auto const &mat : transformPointDofs_) {
      applyBlock(mat, x, y, index, true);
      index += mat.M();
    }
    for (auto const &mat : transformEdgeDofs_) {
      applyBlock(mat, x, y, index, true);
      index += mat.M();
    }
    applyBlock(transformElementDofs_, x, y, index, true);
  }

  This getInverse() {
    This inverse = *this;
    inverse.invert();
    return inverse;
  }

private:
  template <class Value, class Scalar>
  static void assignScaled(Value &out, Scalar factor, Value const &in) {
    if constexpr (requires { out = factor * in; })
      out = factor * in;
    else
      for (size_type i = 0; i < out.size(); ++i)
        assignScaled(out[i], factor, in[i]);
  }

  template <class Value, class Scalar>
  static void addScaled(Value &out, Scalar factor, Value const &in) {
    if constexpr (requires { out += factor * in; })
      out += factor * in;
    else
      for (size_type i = 0; i < out.size(); ++i)
        addScaled(out[i], factor, in[i]);
  }

  template <class Matrix, class VectorIn, class VectorOut>
  static void applyBlock(Matrix const &matrix, VectorIn const &x, VectorOut &y,
                         size_type offset, bool transpose) {
    for (size_type i = 0; i < matrix.N(); ++i) {
      auto coefficient = [&](size_type j) -> auto const & {
        return transpose ? matrix[j][i] : matrix[i][j];
      };
      assignScaled(y[offset + i], coefficient(0), x[offset]);
      for (size_type j = 1; j < matrix.M(); ++j)
        addScaled(y[offset + i], coefficient(j), x[offset + j]);
    }
  }

  void invert() {
    for (auto &mat : transformPointDofs_)
      mat.invert();
    for (auto &mat : transformEdgeDofs_)
      mat.invert();
    transformElementDofs_.invert();
  }
};

/** \brief ArnoldWinther finite element for simplices
 *
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
  using Base = Impl::TransformedFiniteElementMixin<
      This, ReferenceTraits>;
  friend class Impl::TransformedLocalBasis<
      This, ReferenceTraits>;

public:
  /** \brief Export number types, dimensions, etc.
   */

  using Traits = LocalFiniteElementTraits<
      Impl::TransformedLocalBasis<This, ReferenceTraits>,
      Impl::ArnoldWintherLocalCoefficients,
      Impl::ArnoldWintherLocalInterpolation<Element, R>>;

  ArnoldWintherLocalFiniteElement() : Base() {}

  /** \brief Returns the assignment of the degrees of freedom to the element
   * subentities
   */
  const typename Traits::LocalCoefficientsType &localCoefficients() const {
    return coefficients_;
  }

  /** \brief Returns object that evaluates degrees of freedom
   */
  const typename Traits::LocalInterpolationType &localInterpolation() const {
    return interpolation_;
  }

  /** \brief The number of shape functions */
  static constexpr std::size_t size() { return 24; }

  /** \brief The reference element that the local finite element is defined on
   */
  static constexpr GeometryType type() { return GeometryTypes::simplex(2); }

  /**
   * \brief binds the transformation to an element and its elementinformation
   *        Fills the transformation Matrix.
   *
   * \tparam Element
   * \param orientations Face orientations relative to global vertex IDs
   * \param element
   */
  void bind(const ArnoldWintherFaceOrientations &orientations,
            Element const &element) {
    if constexpr (dim != dimWorld)
      DUNE_THROW(Dune::NotImplemented,
                 "Arnold-Winther requires dimension=dimensionworld=2");
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
  // CRTP Interface
  /** \brief Returns the local basis, i.e., the set of shape functions
   */
  Impl::ArnoldWintherReferenceLocalBasis<D, R> const &
  referenceLocalBasis() const {
    return basis_;
  }

  /** Apply the transformation. Note that we do distinguish for
   * Vector/Matrix Type via the DoublePiolas function overload,
   * We assume random access containers.
  */
  template <class InputValues, class OutputValues>
  void transform(InputValues const &inValues, OutputValues &outValues) const {
    using InputValue = typename InputValues::value_type;
    std::vector<InputValue> transformedReferenceValues(size());
    mat_.mtv(inValues, transformedReferenceValues);
    // bind() rejects non-affine geometries, so the Piola map is constant and
    // can be evaluated at any reference point.
    const auto x = Dune::referenceElement<double, 2>(GeometryTypes::simplex(2))
                       .position(0, 0);

    if constexpr (std::is_same_v<
                      InputValue, typename ReferenceTraits::RangeType>)
      DoubleContravariantPiolaTransformator::applyValues(
          transformedReferenceValues, outValues, x, element_->geometry());
    else if constexpr (std::is_same_v<
                           InputValue, typename ReferenceTraits::DivergenceType>)
      DoubleContravariantPiolaTransformator::applyDivergences(
          transformedReferenceValues, outValues, x, element_->geometry());
    else if constexpr (std::is_same_v<
                           InputValue, typename ReferenceTraits::JacobianType>)
      DoubleContravariantPiolaTransformator::applyJacobians(
          transformedReferenceValues, outValues, x, element_->geometry());
    else
      static_assert(Dune::AlwaysFalse<InputValue>::value,
                    "Unsupported Arnold-Winther transformed value type");
  }

private:
  template <class Geometry>
  void fillMatrix(Geometry const &geometry) {
    // Per-edge blocks of P=Q^{-1}.  They convert the four physical edge
    // moments (two Lagrange moments, each with nn and nt components) to the
    // corresponding reference functionals, including edge reorientation.
    std::array<Dune::FieldMatrix<R, 4, 4>, 3> W_k;
    // Symmetric-tensor block of P=Q^{-1}.  Before inversion W represents
    // vec_s(J tau J^T).  The integration element g specializes P to the
    // vertex blocks (g^2 W^{-1}) and cell block (g W^{-1}).
    Dune::FieldMatrix<R, 3, 3> W;
    // The geometry is affine, so any point in the reference triangle is valid.
    auto x = Dune::referenceElement<double, 2>(GeometryTypes::simplex(2))
                 .position(0, 0);
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
      auto tangent =
          refElement.position(upper, 2) - refElement.position(lower, 2);

      tangent /= tangent.two_norm();

      Dune::FieldMatrix<R, 2, 2> referenceG =
          {{-tangent[1], tangent[0]}, {tangent[0], tangent[1]}};
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
    // Fill W  (not yet inverted)
    // TODO this should be improved to handle DiagonalMatrices as well. Since we
    // only have simplices, I think this case currently cannot arise tho.
    // first W tilde
    W[0][0] = jacobian[0][0] * jacobian[0][0];
    W[0][1] = 2. * jacobian[0][0] * jacobian[0][1];
    W[0][2] = jacobian[0][1] * jacobian[0][1];
    W[1][0] = jacobian[0][0] * jacobian[1][0];
    W[1][1] = jacobian[0][0] * jacobian[1][1] +
              jacobian[0][1] * jacobian[1][0];
    W[1][2] = jacobian[0][1] * jacobian[1][1];
    W[2][0] = jacobian[1][0] * jacobian[1][0];
    W[2][1] = 2. * jacobian[1][0] * jacobian[1][1];
    W[2][2] = jacobian[1][1] * jacobian[1][1];
    W.invert();
    // now we have the inverted W breve
    W *= integrationElement * integrationElement;
    // fill matrix
    mat_ = ArnoldWintherBlockDiagonalMatrix<R>{
        std::array<Dune::FieldMatrix<R, 3, 3>, 3>{W, W, W}, W_k,
        W / integrationElement};
  }

private:
  Impl::ArnoldWintherReferenceLocalBasis<D, R> basis_;
  Traits::LocalCoefficientsType coefficients_;
  Traits::LocalInterpolationType interpolation_;
  // Blockmatrix This is the matrix P from the paper mentioned above
  ArnoldWintherBlockDiagonalMatrix<R> mat_;
  ArnoldWintherFaceOrientations faceOrientations_;
  const Element *element_;
};

} // namespace Impl

template <class GV, class R> class ArnoldWintherNode;

template <class GV, typename R>
class ArnoldWintherPreBasis
: public LeafPreBasisMapperMixin<GV>
{
  static const int dim = GV::dimension;
  static_assert(dim == 2,
                "ArnoldWinther PreBasis only implemented for 2d simplices");
  using Base = LeafPreBasisMapperMixin<GV>;

  // helper methods to assign each subentity the number of dofs. Used by the
  // LeafPreBasisMapperMixin.
  static constexpr auto arnoldWintherMapperLayout(Dune::GeometryType type,
                                                  int gridDim) {
    assert(gridDim == 2);
    if (type.isVertex())
      return 3; // three evaluation dof per vertex
    if (type.isLine())
      return 4;
    if ((type.isTriangle()))
      return 3;
    else
      return 0;
  }

public:
  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  using Node = ArnoldWintherNode<GridView, R>;

  //! Constructor for a given grid view object
  ArnoldWintherPreBasis(const GV &gv)
      : Base(gv, arnoldWintherMapperLayout)
  {
    if constexpr (GV::dimension != GV::dimensionworld)
      DUNE_THROW(Dune::NotImplemented,
                 "Arnold-Winther requires dimension=dimensionworld=2.");
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update(GridView const &gv) {
    Base::update(gv);
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const { return Node{this->gridView()}; }

  template <class Element>
  auto faceOrientations(const Element &element) const {
    return Impl::arnoldWintherFaceOrientations(
        element, this->gridView().grid().globalIdSet());
  }
};

template <class GV, class R>
class ArnoldWintherNode : public LeafBasisNode {
public:
  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;

public:
  using FiniteElement =
      Impl::ArnoldWintherLocalFiniteElement<Element, typename GV::ctype, R>;

  ArnoldWintherNode(GV const &gridView)
      : gridView_(&gridView) {
    // this->setSize(finiteElement_.size());
  }

  ~ArnoldWintherNode() {}

  //! Return current element, throw if unbound
  const Element &element() const { return *element_; }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the
   * dune-localfunctions module
   */
  const FiniteElement &finiteElement() const { return finiteElement_; }

  //! Bind to element.
  void bind(const Element &e) {
    if (not e.type().isSimplex())
      DUNE_THROW(Dune::NotImplemented,
                 "ArnoldWintherBasis can only be bound to simplex elements");
    element_ = &e;
    finiteElement_.bind(
        Impl::arnoldWintherFaceOrientations(
            e, gridView_->grid().globalIdSet()),
        *element_);
    this->setSize(finiteElement_.size());
  }

  unsigned int order() const { return 3; }

protected:
  FiniteElement finiteElement_;
  Element const *element_;
  GV const *gridView_;
};

namespace BasisFactory {
/**
 * \brief Create a pre-basis factory that can create ArnoldWinther pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam Range Numbertype used for shape function values
 *
 *
 */

template <typename Range = double> auto arnoldWinther() {
  return [](auto const &gridView) {
    return ArnoldWintherPreBasis<std::decay_t<decltype(gridView)>, Range>(
        gridView);
  };
}
} // namespace BasisFactory
} // namespace Functions
} // namespace Dune

#include <dune/functions/functionspacebases/arnoldwintherbasis_inc.hh>

#endif
