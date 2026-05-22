#ifndef DUNE_C1ELEMENTS_ARNOLDWINTHER_HH
#define DUNE_C1ELEMENTS_ARNOLDWINTHER_HH

#include <dune/common/math.hh>
#include <numeric>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/scalarvectorview.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

#include <dune/functions/common/densevectorview.hh>
#include <dune/functions/common/multidot.hh>
#include <dune/functions/common/mapperutilities.hh>

#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/transformedfiniteelementmixin.hh>

namespace Dune {
namespace Functions {

/**
* \brief Implementation of the conforming Arnold-Winther element
*   This is a finite element used to discretize the stress for two-dimensional
* elasticity, originally proposed in "Arnold, D. N., & Winther, R. (2002).
* Mixed finite elements for elasticity." As such, its shapefunctions take
values
* in the space of symmetric 2x2 matrices, whose (rowwise) divergence is in L2.
* This comes with some complications in the Dune framework, in particular,
* there is no datastructure for 3-tensors, which is the JacobianType of this
finite element.
* This therefore omits the evaluateJacobian method and only implements
evaluateDivergence directly
*   The implementation, in particular the transformation accounting for the
* non-affinity is based on "Aznaran, Francis & Kirby, Robert & Farrell,
Patrick.
* (2021). Transformations for Piola-mapped elements."
*/
// \TODO Rework this to incorporate some stuff for nonaffine mappings
namespace Impl {
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
  /** \brief Double Piola-transform a set of shape-function values
   *
   * \param[in,out] values The values to be Piola-transformed
   */
  template <typename R, int dim, typename LocalCoordinate, typename Geometry>
  static auto
  apply(std::vector<typename Impl::ArnoldWintherTensorTypes<R, dim>::Matrix>
            &values,
        const LocalCoordinate &xi, const Geometry &geometry) {
    auto jacobian = geometry.jacobian(xi);
    auto JT = geometry.jacobianTransposed(xi);

    auto integrationElement = geometry.integrationElement(xi);
    assert(values[0].N() == values[0].M());
    assert(values[0].N() == jacobian.N());
    assert(jacobian.M() == jacobian.N());

    for (auto &value : values) {
      // value = jacobian * value * transpose(jacobian);
      value = multiDot(value, JT, JT);
      value /= (integrationElement * integrationElement);
    }

    return;
  }

  /** \brief Piola-transform a set of shape-function derivatives
   *
   * \param[in,out] gradients The shape function derivatives to be
   * Piola-transformed
   * \TODO This transfroms the divergence. There is no implementation of
   * gradients We want the physical divergence of tau \f$ div \tau =
   * \frac{1}{(det J)^2}J\hat{div}\tau  \f$ \bug The current implementation
   * works only for affine geometries. The Piola transformation for non-affine
   * geometries requires second derivatives of the geometry, which we don't get
   *   from the dune-grid Geometry interface.
   */
  template <typename R, int dim, typename LocalCoordinate, typename Geometry>
  static auto
  apply(std::vector<typename Impl::ArnoldWintherTensorTypes<R, dim>::Vector>
            &divergences,
        const LocalCoordinate &xi, const Geometry &geometry) {
    auto jacobian = geometry.jacobian(xi);
    auto integrationElement2 =
        geometry.integrationElement(xi) * geometry.integrationElement(xi);
    assert(dim == jacobian.N());
    assert(jacobian.M() ==
           jacobian.N()); // \TODO think this through for (linear) surface
                          // meshs. Maybe this works

    for (auto &value : divergences) {
      auto tmp = value;
      value = 0;
      for (std::size_t k = 0; k < dim; k++) {
        for (auto &&[jacobian_k_i, i] : sparseRange(jacobian[k]))
          value[k] += jacobian_k_i * tmp[i];
        value[k] /= integrationElement2;
      }
    }
  }

  /** \brief Wrapper around a callable that applies the inverse Piola
   * transform
   *
   * The LocalInterpolation implementations in dune-localfunctions expect
   * local-valued functions, but the ones dune-functions expect
   * global-valued ones.  Therefore, we need to stuff the inverse Piola
   * transform between dune-functions and dune-localfunctions, and this is
   * what this class does.
   */
  template <class Function, class LocalCoordinate, class Element>
  class LocalValuedFunction {
    const Function &f_;
    const Element &element_;

  public:
    LocalValuedFunction(const Function &f, const Element &element)
        : f_(f), element_(element) {}

    auto operator()(const LocalCoordinate &xi) const {
      auto globalValue = f_(xi);

      // Apply the inverse Piola transform
      auto jacobianInverse = element_.geometry().jacobianInverse(xi);
      auto integrationElement = element_.geometry().integrationElement(xi);

      globalValue = jacobianInverse * globalValue * transpose(jacobianInverse);

      globalValue *= integrationElement * integrationElement;

      return globalValue;
    }
  };
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

      std::size_t lower = (i == 2) ? 1 : 0;
      std::size_t upper = (i == 0) ? 1 : 2;
      auto tangent =
          refElement.position(upper, 2) - refElement.position(lower, 2);
      tangent /= tangent.two_norm();
      std::decay_t<decltype(tangent)> normal = {
          -tangent[1], tangent[0]};

      using fRange = std::decay_t<std::remove_cv_t<decltype(f(refElement.position(i, dim)))>>;
      using protomotedType =
          typename PromotionTraits<typename FieldTraits<fRange>::field_type,
                                   c_type>::PromotedType;

      FieldVector<protomotedType, 2> tmp;
      for (auto&& val : moments){
        val.mtv(normal, tmp);

        it[0] = dot(tmp, normal)*refElement.template geometry<1>(i).volume();
        it[1] = dot(tmp, tangent)*refElement.template geometry<1>(i).volume();
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
  using GlobalCoordinate = typename Element::Geometry::GlobalCoordinate;

  using ctype = typename Element::Geometry::ctype;
  static constexpr size_type dim = Element::Geometry::mydimension;
  static constexpr size_type dimWorld = Element::Geometry::coordimension;
  static constexpr int size = 24; // number of dofs. TODO generalize this!
public:
  ArnoldWintherLocalInterpolation(int quadOrder = 10)
  : quadratureOrder(quadOrder)
  {}

  void bind(std::bitset<3> data, Element const& e){
    edgeOrientation_ = data;
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
    auto lf = DoubleContravariantPiolaTransformator::LocalValuedFunction<F, typename Element::Geometry::LocalCoordinate, Element>(f, *element);
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

      std::size_t lower = (i == 2) ? 1 : 0;
      std::size_t upper = (i == 0) ? 1 : 2;
      auto tangent = (*element).template subEntity<dim>(upper).geometry().center() - (*element).template subEntity<dim>(lower).geometry().center();
      tangent /= tangent.two_norm();
      std::decay_t<decltype(tangent)> normal = {
          -tangent[1], tangent[0]};

      using fRange = typename std::decay_t<std::remove_cv_t<decltype(f(std::declval<GlobalCoordinate>()))>>;
      using protomotedType =
          typename PromotionTraits<typename FieldTraits<fRange>::field_type,
                                   ctype>::PromotedType;

      FieldVector<protomotedType, 2> normalTimesMoment;
      for (std::size_t m = 0; m < momentOrder + 1; ++m){
        if (edgeOrientation_[i])
          moments[momentOrder -m].mtv(normal, normalTimesMoment);
        else
          moments[m].mtv(normal, normalTimesMoment);
        it[0] = dot(normalTimesMoment, normal)*edgeGeo.volume();
        it[1] = dot(normalTimesMoment, tangent)*edgeGeo.volume();
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
  std::bitset<3> edgeOrientation_;
  Element const* element = nullptr;

};

// Class offering reading and writing access to a Vector starting from an index
// \TODO check that this actually works also in debug mode or drop this
// \TODO maybe make this private class of BlockDiagonalMatrix
template <class Vector>
class VectorSlice {
public:
  using value_type = typename Vector::value_type;
  using size_type = typename Vector::size_type;

  VectorSlice() = delete;
  VectorSlice(Vector &vec, size_type index)
      : VectorSlice(vec, index, vec.size()) {}
  VectorSlice(Vector &vec, size_type index, size_type end)
      : vec_(vec), i_(index), end_(end) {}
  VectorSlice &operator=(value_type scalar) {
    for (size_type i = 0u; i < size(); ++i)
      vec_[i_ + i] = scalar;
    return *this;
  }
  value_type const &operator[](size_type const &index) const {
    return vec_[i_ + index];
  }
  value_type &operator[](size_type const &index) { return vec_[i_ + index]; }

  size_type N() const { return end_ - i_; }

  size_type size() const { return end_ - i_; }

private:
  Vector &vec_;
  size_type i_;
  size_type end_;
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
    auto &&xx = Dune::Impl::asVector(x);
    auto &&yy = Dune::Impl::asVector(y);
    DUNE_ASSERT_BOUNDS((void *)(&x) != (void *)(&y));
    DUNE_ASSERT_BOUNDS(xx.N() == M());
    DUNE_ASSERT_BOUNDS(yy.N() == N());

    size_type index = 0;
    for (auto const &mat : transformPointDofs_) {
      auto res = VectorSlice(y, index, index + 3);
      res = 0;
      mat.mv(VectorSlice(x, index, index + 3), res);
      index += mat.M();
    }
    for (auto const &mat : transformEdgeDofs_) {
      auto res = VectorSlice(y, index, index + 4);
      res = 0;
      mat.mv(VectorSlice(x, index, index + 4), res);
      index += mat.M();
    }

    auto res = VectorSlice(y, index, index + 3);
    res = 0;
    transformElementDofs_.mv(VectorSlice(x, index, index + 3), res);
  }

  template <class VectorIn, class VectorOut>
  void mtv(VectorIn const &x, VectorOut &y) const {
    auto &&xx = Dune::Impl::asVector(x);
    auto &&yy = Dune::Impl::asVector(y);
    DUNE_ASSERT_BOUNDS((void *)(&x) != (void *)(&y));
    DUNE_ASSERT_BOUNDS(xx.N() == M());
    DUNE_ASSERT_BOUNDS(yy.N() == N());

    // using y_field_type = typename FieldTraits<VectorOut>::field_type;
    size_type index = 0;
    for (auto const &mat : transformPointDofs_) {
      auto res = VectorSlice(y, index, index + 3);
      res = 0;
      mat.mtv(VectorSlice(x, index, index + 3), res);
      index += mat.M();
    }
    for (auto const &mat : transformEdgeDofs_) {
      auto res = VectorSlice(y, index, index + 4);
      res = 0;
      mat.mtv(VectorSlice(x, index, index + 4), res);
      index += mat.M();
    }

    auto res = VectorSlice(y, index, index + 3);
    res = 0;
    transformElementDofs_.mtv(VectorSlice(x, index, index + 3), res);
  }

  This getInverse() {
    This inverse = *this;
    inverse.invert();
    return inverse;
  }

private:
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
  using Base = Impl::TransformedFiniteElementMixin<
      ArnoldWintherLocalFiniteElement<Element, D, R>,
      typename ArnoldWintherReferenceLocalBasis<D, R>::Traits>;
  friend class Impl::TransformedLocalBasis<
      ArnoldWintherLocalFiniteElement<Element, D, R>,
      typename ArnoldWintherReferenceLocalBasis<D, R>::Traits>;

public:
  /** \brief Export number types, dimensions, etc.
   */

  using Traits = LocalFiniteElementTraits<
      Impl::ArnoldWintherReferenceLocalBasis<D, R>,
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
   * \param data     Edge orientations
   * \param element
   */
  void bind(std::bitset<3> data, Element const &element) {
    edgeOrientation_ = data;
    element_ = &element;
    interpolation_.bind(data, element);
    fillMatrix(element.geometry()); // barycenter, because we need some value.
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

  // Here we cannot directly use
  // mat_.mtv(inValues, outValues);
  // because mv expects the DenseVector interface.
  auto inValuesDenseVector = Impl::DenseVectorView(inValues);
  auto outValuesDenseVector = Impl::DenseVectorView(outValues);
  mat_.mtv(inValuesDenseVector, outValuesDenseVector);


  auto x = Dune::referenceElement<double, 2>(GeometryTypes::simplex(2))
                .position(0, 0);
  DoubleContravariantPiolaTransformator::apply<R, 2>(outValues, x,
                                                    element_->geometry());
  }

private:
  template <class Geometry>
  void fillMatrix(Geometry const &geometry) {
    std::array<R, 3> alpha;
    std::array<R, 3> beta;

    std::array<Dune::FieldVector<R, 2>, 3> referenceTangents; // normalized

    std::array<R, 3> referenceEdgeLength;
    std::array<R, 3> globalEdgeLength;

    std::array<Dune::FieldMatrix<R, 2, 2>, 3> referenceG;

    std::array<Dune::FieldMatrix<R, 4, 4>, 3> W_k;
    Dune::FieldMatrix<R, 3, 3> W;
    // \TODO For linear triangles any point inside is fine, for curved ones one
    // would need to chose them according to the dofs
    auto x = Dune::referenceElement<double, 2>(GeometryTypes::simplex(2))
                 .position(0, 0);
    auto jacobianTransposed = geometry.jacobianTransposed(x);
    auto jacobianDeterminant = jacobianTransposed.determinant(); //geometry.integrationElement(x);

    // By default, edges point from the vertex with the smaller index
    // to the vertex with the larger index. Note that the alpha and beta are
    // invariant of orientation, since the normals/tangents appear twice in
    // their definitions.

    // get local and global Tangents
    auto refElement = Dune::referenceElement<double, 2>(geometry.type());
    for (std::size_t i = 0; i < 3; ++i) {
      std::size_t lower = (i == 2) ? 1 : 0;
      std::size_t upper = (i == 0) ? 1 : 2;
      auto tangent =
          refElement.position(upper, 2) - refElement.position(lower, 2);

      referenceEdgeLength[i] = tangent.two_norm();
      tangent /= referenceEdgeLength[i];

      auto globalEdge = geometry.global(refElement.position(upper, 2)) -
                        geometry.global(refElement.position(lower, 2));

      globalEdgeLength[i] = globalEdge.two_norm();

      referenceG[i] = {{-tangent[1], tangent[0]}, {tangent[0], tangent[1]}};
      auto tmp = tangent, tmp2 = tangent;
      jacobianTransposed.mtv(tangent, tmp);
      jacobianTransposed.mv(tmp, tmp2);
      referenceG[i].mtv(tmp2, tmp);
      // auto tmp = transpose(referenceG[i]) * jacobianTransposed * (transpose(jacobianTransposed) * tangent);
      alpha[i] = tmp[0] / jacobianDeterminant;
      beta[i] = tmp[1] / jacobianDeterminant;



      // already inverted W_k
      // By using Lagrange moments these matrices are orientation invariant.
      W_k[i] = 0;
      if (edgeOrientation_[i]){

        W_k[i][2][0] = 1.;
        W_k[i][3][0] = -alpha[i] / beta[i];
        W_k[i][3][1] = 1. / beta[i];
        W_k[i][0][2] = 1.;
        W_k[i][1][2] = -alpha[i] / beta[i];
        W_k[i][1][3] = 1. / beta[i];

      }
      else{
        W_k[i][0][0] = 1.;
        W_k[i][1][0] = -alpha[i] / beta[i];
        W_k[i][1][1] = 1. / beta[i];
        W_k[i][2][2] = 1.;
        W_k[i][3][2] = -alpha[i] / beta[i];
        W_k[i][3][3] = 1. / beta[i];
      }


      W_k[i] *= globalEdgeLength[i] / referenceEdgeLength[i];
    }
    // Fill W  (not yet inverted)
    // TODO this should be improved to handle DiagonalMatrices as well. Since we
    // only have simplices, I think this case currently cannot arise tho.
    // first W tilde
    W[0][0] = jacobianTransposed[0][0] * jacobianTransposed[0][0];
    W[0][1] = 2. * jacobianTransposed[0][0] * jacobianTransposed[1][0];
    W[0][2] = jacobianTransposed[1][0] * jacobianTransposed[1][0];
    W[1][0] = jacobianTransposed[0][0] * jacobianTransposed[0][1];
    W[1][1] = jacobianTransposed[0][0] * jacobianTransposed[1][1] +
              jacobianTransposed[1][0] * jacobianTransposed[0][1];
    W[1][2] = jacobianTransposed[1][0] * jacobianTransposed[1][1];
    W[2][0] = jacobianTransposed[0][1] * jacobianTransposed[0][1];
    W[2][1] = 2. * jacobianTransposed[0][1] * jacobianTransposed[1][1];
    W[2][2] = jacobianTransposed[1][1] * jacobianTransposed[1][1];
    W.invert();
    // now we have the inverted W breve
    W *= jacobianDeterminant * jacobianDeterminant;
    // fill matrix
    mat_ = ArnoldWintherBlockDiagonalMatrix<R>{
        std::array<Dune::FieldMatrix<R, 3, 3>, 3>{W, W, W}, W_k,
        W / jacobianDeterminant};
  }

private:
  Impl::ArnoldWintherReferenceLocalBasis<D, R> basis_;
  Traits::LocalCoefficientsType coefficients_;
  Traits::LocalInterpolationType interpolation_;
  // Blockmatrix This is the matrix P from the paper mentioned above
  ArnoldWintherBlockDiagonalMatrix<R> mat_;
  std::bitset<3> edgeOrientation_;
  const Element *element_;
};

} // namespace Impl

template <class GV, class R> class ArnoldWintherNode;

template <class GV, typename R>
class ArnoldWintherPreBasis
: public LeafPreBasisMapperMixin<GV>//, Impl::ModuloEdgeTwist<typename GV::IndexSet>>
{
  static const int dim = GV::dimension;
  static_assert(dim == 2,
                "ArnoldWinther PreBasis only implemented for 2d simplices");
  using Base = LeafPreBasisMapperMixin<GV>;//, Impl::ModuloEdgeTwist<typename GV::IndexSet>>;

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

  using SubEntityMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;

public:
  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  using Node = ArnoldWintherNode<GridView, R>;

  //! Constructor for a given grid view object
  ArnoldWintherPreBasis(const GV &gv)
      : Base(gv, arnoldWintherMapperLayout),//, Impl::ModuloEdgeTwist{gv.indexSet(), arnoldWintherMapperLayout(GeometryTypes::line, dim), 2}),
        mapper_({gv, mcmgElementLayout()})
  {
    data_ = Impl::computeEdgeOrientations(mapper_);
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update(GridView const &gv) {
    Base::update(gv);
    mapper_.update(this->gridView());
    data_ = Impl::computeEdgeOrientations(mapper_);
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const { return Node{mapper_, data_}; }

private:
  SubEntityMapper mapper_;
public:
  std::vector<std::bitset<3>> data_;
};

template <class GV, class R>
class ArnoldWintherNode : public LeafBasisNode {
public:
  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;

public:
  using FiniteElement =
      Impl::ArnoldWintherLocalFiniteElement<Element, typename GV::ctype, R>;

  ArnoldWintherNode(Mapper const &m, std::vector<std::bitset<3>> const &data)
      : mapper_(&m), data_(&data) {
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
    finiteElement_.bind((*data_)[mapper_->index(e)], *element_);
    this->setSize(finiteElement_.size());
  }

  unsigned int order() const { return 3; }

protected:
  FiniteElement finiteElement_;
  Element const *element_;
  Mapper const *mapper_;
  std::vector<std::bitset<3>> const *data_;
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
template <class T>
struct FieldTraits<typename Functions::Impl::VectorSlice<T>> {
  typedef typename FieldTraits<T>::field_type field_type;
  typedef typename FieldTraits<T>::real_type real_type;
};
} // namespace Dune

#include <dune/functions/functionspacebases/arnoldwintherbasis_inc.hh>

#endif