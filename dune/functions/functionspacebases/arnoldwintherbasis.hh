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

#include <dune/functions/common/densevectorview.hh>

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
template <class D, class R> class ArnoldWintherReferenceLocalBasis {
  // only implemented for triangles
  // only lowest order implemented

public:
  static constexpr unsigned int dim = 2;
  // TODO generalize this to surface maps
  static constexpr unsigned int dimAmbient = 2;

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
    using DomainType = Dune::FieldVector<D, dimDomain>;

    //! \brief Export type for range field
    using RangeFieldType = R;

    //! \brief dimension of the range
    // TODO discuss what this should be. For now we take it as the entries in
    // the range matrix
    constexpr static int dimRange = dimAmbient * dimAmbient;

    //! \brief range type
    using RangeType = ArnoldWintherTensorTypes<R, dimAmbient>::Matrix;

    /** \brief Type to represent derivative
     */
    using JacobianType =
        ArnoldWintherTensorTypes<R, dimAmbient, dim>::ThreeTensor;

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
  ArnoldWintherReferenceLocalInterpolation(int quadOrder = 6)
  : quadratureOrder(quadOrder){}

  /** \brief Evaluate a given function at the Lagrange nodes
   *
   * \tparam F Type of function to evaluate
   * \tparam C Type used for the values of the function
   * \param[in] ff Function to evaluate
   * \param[out] out Array of function values
   */
  template <typename F, typename C>
  void interpolate(const F &ff, std::vector<C> &out) const {
    auto &&f = Dune::Impl::makeFunctionWithCallOperator<LocalCoordinate>(ff);

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
      auto average = integralMoment<0>(f, refElement.template geometry<1>(i), quadratureOrder);
      auto firstMoment =
          integralMoment<1>(f, refElement.template geometry<1>(i), quadratureOrder+1);

      std::size_t lower = (i == 2) ? 1 : 0;
      std::size_t upper = (i == 0) ? 1 : 2;
      auto tangent =
          refElement.position(upper, 2) - refElement.position(lower, 2);
      tangent /= tangent.two_norm();
      std::decay_t<decltype(tangent)> normal = {
          tangent[1], -tangent[0]}; // rotation by -90 degree

      using fRange = std::decay_t<std::remove_cv_t<decltype(f(refElement.position(i, dim)))>>;
      using protomotedType =
          typename PromotionTraits<typename FieldTraits<fRange>::field_type,
                                   c_type>::PromotedType;

      FieldVector<protomotedType, 2> normalTimesAverage, normalTimesFirstMoment;
      average.mtv(normal, normalTimesAverage);
      firstMoment.mv(normal, normalTimesFirstMoment);
      it[0] = dot(normalTimesAverage, normal);
      it[1] = dot(normalTimesAverage, tangent);
      it[2] = dot(normalTimesFirstMoment, normal);
      it[3] = dot(normalTimesFirstMoment, tangent);
      it += 4;
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

// generated with sympy from symfem library
template <class D, class R>
void ArnoldWintherReferenceLocalBasis<D, R>::evaluateFunction(
    const typename Traits::DomainType &in,
    std::vector<typename Traits::RangeType> &out) const {
  out.resize(size());
  auto iter = out.begin(); // generated with sympy from symfem library
  auto const &x = in[0], y = in[1];
  double diag;
  // 1th basis function
  diag = x * (-6 * x * y + y * (6 - 9 * y));

  *(iter++) = {
      {x * x * (2 * x + 9 * y - 3) + y * (y * (18 - 10 * y) - 9) + 1, diag},
      {diag, 6 * x * y * y + y * y * (3 * y - 3)}};

  // 2th basis function
  diag = x * (x * (-10 * x - 45 * y + 18) + y * (48 - 45 * y) - 9) +
         y * (y * (18 - 10 * y) - 9) + 1;

  *(iter++) = {{x * (x * (15 * x + 45 * y - 24) + y * (30 * y - 36) + 9), diag},
               {diag, x * (30 * x * y + y * (45 * y - 36)) +
                          y * (y * (15 * y - 24) + 9)}};

  // 3th basis function
  diag = x * (-9 * x * y + y * (6 - 6 * y));

  *(iter++) = {{x * x * (3 * x + 6 * y - 3), diag},
               {diag, x * (x * (18 - 10 * x) + 9 * y * y - 9) +
                          y * y * (2 * y - 3) + 1}};

  // 4th basis function
  diag = x * (6 * x * y + y * (9 * y - 6));

  *(iter++) = {{x * x * (-2 * x - 9 * y + 3), diag},
               {diag, -6 * x * y * y + y * y * (3 - 3 * y)}};

  // 5th basis function
  diag = x * (x * (10 * x + 27 * y - 12) + y * (18 * y - 18) + 3);

  *(iter++) = {{x * x * (-9 * x - 18 * y + 9), diag},
               {diag, x * (-30 * x * y + y * (24 - 27 * y)) +
                          y * (y * (9 - 6 * y) - 3)}};

  // 6th basis function
  diag = 0;

  *(iter++) = {{0, diag}, {diag, x * (x * (10 * x - 12) + 3)}};

  // 7th basis function
  diag = 0;

  *(iter++) = {{y * (y * (10 * y - 12) + 3), diag}, {diag, 0}};

  // 8th basis function
  diag = x * (18 * x * y + y * (27 * y - 18)) + y * (y * (10 * y - 12) + 3);

  *(iter++) = {{x * (x * (-6 * x - 27 * y + 9) + y * (24 - 30 * y) - 3), diag},
               {diag, -18 * x * y * y + y * y * (9 - 9 * y)}};

  // 9th basis function
  diag = x * (9 * x * y + y * (6 * y - 6));

  *(iter++) = {{x * x * (-3 * x - 6 * y + 3), diag},
               {diag, -9 * x * y * y + y * y * (3 - 2 * y)}};

  // 10th basis function
  diag = x * (54 * x * y + y * (36 * y - 36));

  *(iter++) = {{x * x * (-18 * x - 36 * y + 18), diag},
               {diag, x * (x * (60 * x - 96) + y * (24 - 54 * y) + 36) +
                          y * (y * (30 - 12 * y) - 18)}};

  // 11th basis function
  diag = x * (x * (-60 * x - 180 * y + 96) + y * (132 - 90 * y) - 36);

  *(iter++) = {{x * (x * (60 * x + 90 * y - 78) - 24 * y + 18), diag},
               {diag, x * (180 * x * y + y * (180 * y - 192)) +
                          y * (y * (30 * y - 66) + 36)}};

  // 12th basis function
  diag = x * (-108 * x * y + y * (72 - 72 * y));

  *(iter++) = {{x * x * (36 * x + 72 * y - 36), diag},
               {diag, x * (x * (180 - 120 * x) + y * (108 * y - 48) - 60) +
                          y * (y * (24 * y - 48) + 24)}};

  // 13th basis function
  diag = x * (x * (120 * x + 360 * y - 180) + y * (180 * y - 240) + 60);

  *(iter++) = {{x * (x * (-120 * x - 180 * y + 144) + 24 * y - 24), diag},
               {diag, x * (-360 * x * y + y * (360 - 360 * y)) +
                          y * (y * (120 - 60 * y) - 60)}};

  // 14th basis function
  diag = x * (36 * x * y + y * (54 * y - 36));

  *(iter++) = {{x * (x * (-12 * x - 54 * y + 30) + 24 * y - 18) +
                    y * (y * (60 * y - 96) + 36),
                diag},
               {diag, -36 * x * y * y + y * y * (18 - 18 * y)}};

  // 15th basis function
  diag = x * (90 * x * y + y * (180 * y - 132)) + y * (y * (60 * y - 96) + 36);

  *(iter++) = {
      {x * (x * (-30 * x - 180 * y + 66) + y * (192 - 180 * y) - 36), diag},
      {diag, x * y * (24 - 90 * y) + y * (y * (78 - 60 * y) - 18)}};

  // 16th basis function
  diag = x * (-72 * x * y + y * (72 - 108 * y));

  *(iter++) = {{x * (x * (24 * x + 108 * y - 48) - 48 * y + 24) +
                    y * (y * (180 - 120 * y) - 60),
                diag},
               {diag, 72 * x * y * y + y * y * (36 * y - 36)}};

  // 17th basis function
  diag =
      x * (-180 * x * y + y * (240 - 360 * y)) + y * (y * (180 - 120 * y) - 60);

  *(iter++) = {
      {x * (x * (60 * x + 360 * y - 120) + y * (360 * y - 360) + 60), diag},
      {diag, x * y * (180 * y - 24) + y * (y * (120 * y - 144) + 24)}};

  // 18th basis function
  diag = x * (9 * M_SQRT2 * x * y - 9 * M_SQRT2 * y * y);

  *(iter++) = {{x * (x * (-3 * M_SQRT2 * x + 9 * M_SQRT2 * y + 6 * M_SQRT2) -
                     3 * M_SQRT2),
                diag},
               {diag, x * y * (-9 * M_SQRT2 * y + 12 * M_SQRT2) +
                          y * (3 * M_SQRT2 * y * y - 3 * M_SQRT2)}};

  // 19th basis function
  diag = x * (45 * M_SQRT2 * x * y + y * (45 * M_SQRT2 * y - 36 * M_SQRT2));

  *(iter++) = {
      {x * (x * (-15 * M_SQRT2 * x - 45 * M_SQRT2 * y + 12 * M_SQRT2) +
            3 * M_SQRT2),
       diag},
      {diag, x * y * (-45 * M_SQRT2 * y + 12 * M_SQRT2) +
                 y * (y * (-15 * M_SQRT2 * y + 18 * M_SQRT2) - 3 * M_SQRT2)}};

  // 20th basis function
  diag = x * (-18 * M_SQRT2 * x * y + 18 * M_SQRT2 * y * y);

  *(iter++) = {{x * (x * (6 * M_SQRT2 * x - 18 * M_SQRT2 * y - 6 * M_SQRT2) +
                     12 * M_SQRT2 * y),
                diag},
               {diag, x * y * (18 * M_SQRT2 * y - 12 * M_SQRT2) +
                          y * y * (-6 * M_SQRT2 * y + 6 * M_SQRT2)}};

  // 21th basis function
  diag = x * (-90 * M_SQRT2 * x * y + y * (-90 * M_SQRT2 * y + 72 * M_SQRT2));

  *(iter++) = {{x * (x * (30 * M_SQRT2 * x + 90 * M_SQRT2 * y - 30 * M_SQRT2) -
                     12 * M_SQRT2 * y),
                diag},
               {diag, x * y * (90 * M_SQRT2 * y - 12 * M_SQRT2) +
                          y * y * (30 * M_SQRT2 * y - 30 * M_SQRT2)}};

  // 22th basis function
  diag = 0;

  *(iter++) = {{x * (-24 * x - 24 * y + 24), diag}, {diag, 0}};

  // 23th basis function
  diag = 24 * x * y;

  *(iter++) = {{x * (-24 * x - 48 * y + 24), diag},
               {diag, -48 * x * y + y * (24 - 24 * y)}};

  // 24th basis function
  diag = 0;

  *(iter++) = {{0, diag}, {diag, -24 * x * y + y * (24 - 24 * y)}};
}
// generated with sympy from symfem library
template <class D, class R>
void ArnoldWintherReferenceLocalBasis<D, R>::evaluateDivergence(
    const typename Traits::DomainType &in,
    std::vector<typename Traits::DivergenceType> &out) const {
  out.resize(size());
  auto iter = out.begin(); // generated with sympy from symfem library
  auto const &x = in[0], y = in[1];
  double diag;
  // 1th basis function

  *(iter++) = {0, 0};

  // 2th basis function

  *(iter++) = {0, 0};

  // 3th basis function

  *(iter++) = {0, 0};

  // 4th basis function

  *(iter++) = {0, 0};

  // 5th basis function

  *(iter++) = {0, 0};

  // 6th basis function

  *(iter++) = {0, 0};

  // 7th basis function

  *(iter++) = {0, 0};

  // 8th basis function

  *(iter++) = {0, 0};

  // 9th basis function

  *(iter++) = {0, 0};

  // 10th basis function

  *(iter++) = {0, 24 * x + 24 * y - 18};

  // 11th basis function

  *(iter++) = {-24 * x - 24 * y + 18, 0};

  // 12th basis function

  *(iter++) = {0, -48 * x - 24 * y + 24};

  // 13th basis function

  *(iter++) = {48 * x + 24 * y - 24, 0};

  // 14th basis function

  *(iter++) = {24 * x + 24 * y - 18, 0};

  // 15th basis function

  *(iter++) = {0, 24 * x + 24 * y - 18};

  // 16th basis function

  *(iter++) = {-24 * x - 48 * y + 24, 0};

  // 17th basis function

  *(iter++) = {0, -24 * x - 48 * y + 24};

  // 18th basis function

  *(iter++) = {12 * M_SQRT2 * x - 3 * M_SQRT2, 12 * M_SQRT2 * x - 3 * M_SQRT2};

  // 19th basis function

  *(iter++) = {-12 * M_SQRT2 * x + 3 * M_SQRT2, 12 * M_SQRT2 * x - 3 * M_SQRT2};

  // 20th basis function

  *(iter++) = {-12 * M_SQRT2 * x + 12 * M_SQRT2 * y,
               -12 * M_SQRT2 * x + 12 * M_SQRT2 * y};

  // 21th basis function

  *(iter++) = {12 * M_SQRT2 * x - 12 * M_SQRT2 * y,
               -12 * M_SQRT2 * x + 12 * M_SQRT2 * y};

  // 22th basis function

  *(iter++) = {-48 * x - 24 * y + 24, 0};

  // 23th basis function

  *(iter++) = {-24 * x - 48 * y + 24, -48 * x - 24 * y + 24};

  // 24th basis function

  *(iter++) = {0, -24 * x - 48 * y + 24};
}
// generated with sympy from symfem library
template <class D, class R>
void ArnoldWintherReferenceLocalBasis<D, R>::evaluateJacobian(
    const typename Traits::DomainType &in,
    std::vector<typename Traits::JacobianType> &out) const {
  out.resize(size());
  auto iter = out.begin(); // generated with sympy from symfem library
  auto const &x = in[0], y = in[1];
  FieldVector<double, 2> diag;
  // 1th basis function
  diag = {-12 * x * y + y * (6 - 9 * y), x * (-6 * x - 18 * y + 6)};

  *(iter++) = {
      {{x * (6 * x + 18 * y - 6), 9 * x * x + y * (36 - 30 * y) - 9}, diag},
      {diag, {6 * y * y, 12 * x * y + y * (9 * y - 6)}}};

  // 2th basis function
  diag = {x * (-30 * x - 90 * y + 36) + y * (48 - 45 * y) - 9,
          x * (-45 * x - 90 * y + 48) + y * (36 - 30 * y) - 9};

  *(iter++) = {{{x * (45 * x + 90 * y - 48) + y * (30 * y - 36) + 9,
                 x * (45 * x + 60 * y - 36)},
                diag},
               {diag,
                {60 * x * y + y * (45 * y - 36),
                 x * (30 * x + 90 * y - 36) + y * (45 * y - 48) + 9}}};

  // 3th basis function
  diag = {-18 * x * y + y * (6 - 6 * y), x * (-9 * x - 12 * y + 6)};

  *(iter++) = {
      {{x * (9 * x + 12 * y - 6), 6 * x * x}, diag},
      {diag,
       {x * (36 - 30 * x) + 9 * y * y - 9, 18 * x * y + y * (6 * y - 6)}}};

  // 4th basis function
  diag = {12 * x * y + y * (9 * y - 6), x * (6 * x + 18 * y - 6)};

  *(iter++) = {{{x * (-6 * x - 18 * y + 6), -9 * x * x}, diag},
               {diag, {-6 * y * y, -12 * x * y + y * (6 - 9 * y)}}};

  // 5th basis function
  diag = {x * (30 * x + 54 * y - 24) + y * (18 * y - 18) + 3,
          x * (27 * x + 36 * y - 18)};

  *(iter++) = {{{x * (-27 * x - 36 * y + 18), -18 * x * x}, diag},
               {diag,
                {-60 * x * y + y * (24 - 27 * y),
                 x * (-30 * x - 54 * y + 24) + y * (18 - 18 * y) - 3}}};

  // 6th basis function
  diag = {0, 0};

  *(iter++) = {{{0, 0}, diag}, {diag, {x * (30 * x - 24) + 3, 0}}};

  // 7th basis function
  diag = {0, 0};

  *(iter++) = {{{0, y * (30 * y - 24) + 3}, diag}, {diag, {0, 0}}};

  // 8th basis function
  diag = {36 * x * y + y * (27 * y - 18),
          x * (18 * x + 54 * y - 18) + y * (30 * y - 24) + 3};

  *(iter++) = {{{x * (-18 * x - 54 * y + 18) + y * (24 - 30 * y) - 3,
                 x * (-27 * x - 60 * y + 24)},
                diag},
               {diag, {-18 * y * y, -36 * x * y + y * (18 - 27 * y)}}};

  // 9th basis function
  diag = {18 * x * y + y * (6 * y - 6), x * (9 * x + 12 * y - 6)};

  *(iter++) = {{{x * (-9 * x - 12 * y + 6), -6 * x * x}, diag},
               {diag, {-9 * y * y, -18 * x * y + y * (6 - 6 * y)}}};

  // 10th basis function
  diag = {108 * x * y + y * (36 * y - 36), x * (54 * x + 72 * y - 36)};

  *(iter++) = {{{x * (-54 * x - 72 * y + 36), -36 * x * x}, diag},
               {diag,
                {x * (180 * x - 192) + y * (24 - 54 * y) + 36,
                 x * (24 - 108 * y) + y * (60 - 36 * y) - 18}}};

  // 11th basis function
  diag = {x * (-180 * x - 360 * y + 192) + y * (132 - 90 * y) - 36,
          x * (-180 * x - 180 * y + 132)};

  *(iter++) = {
      {{x * (180 * x + 180 * y - 156) - 24 * y + 18, x * (90 * x - 24)}, diag},
      {diag,
       {360 * x * y + y * (180 * y - 192),
        x * (180 * x + 360 * y - 192) + y * (90 * y - 132) + 36}}};

  // 12th basis function
  diag = {-216 * x * y + y * (72 - 72 * y), x * (-108 * x - 144 * y + 72)};

  *(iter++) = {{{x * (108 * x + 144 * y - 72), 72 * x * x}, diag},
               {diag,
                {x * (360 - 360 * x) + y * (108 * y - 48) - 60,
                 x * (216 * y - 48) + y * (72 * y - 96) + 24}}};

  // 13th basis function
  diag = {x * (360 * x + 720 * y - 360) + y * (180 * y - 240) + 60,
          x * (360 * x + 360 * y - 240)};

  *(iter++) = {
      {{x * (-360 * x - 360 * y + 288) + 24 * y - 24, x * (24 - 180 * x)},
       diag},
      {diag,
       {-720 * x * y + y * (360 - 360 * y),
        x * (-360 * x - 720 * y + 360) + y * (240 - 180 * y) - 60}}};

  // 14th basis function
  diag = {72 * x * y + y * (54 * y - 36), x * (36 * x + 108 * y - 36)};

  *(iter++) = {{{x * (-36 * x - 108 * y + 60) + 24 * y - 18,
                 x * (24 - 54 * x) + y * (180 * y - 192) + 36},
                diag},
               {diag, {-36 * y * y, -72 * x * y + y * (36 - 54 * y)}}};

  // 15th basis function
  diag = {180 * x * y + y * (180 * y - 132),
          x * (90 * x + 360 * y - 132) + y * (180 * y - 192) + 36};

  *(iter++) = {
      {{x * (-90 * x - 360 * y + 132) + y * (192 - 180 * y) - 36,
        x * (-180 * x - 360 * y + 192)},
       diag},
      {diag,
       {y * (24 - 90 * y), x * (24 - 180 * y) + y * (156 - 180 * y) - 18}}};

  // 16th basis function
  diag = {-144 * x * y + y * (72 - 108 * y), x * (-72 * x - 216 * y + 72)};

  *(iter++) = {{{x * (72 * x + 216 * y - 96) - 48 * y + 24,
                 x * (108 * x - 48) + y * (360 - 360 * y) - 60},
                diag},
               {diag, {72 * y * y, 144 * x * y + y * (108 * y - 72)}}};

  // 17th basis function
  diag = {-360 * x * y + y * (240 - 360 * y),
          x * (-180 * x - 720 * y + 240) + y * (360 - 360 * y) - 60};

  *(iter++) = {
      {{x * (180 * x + 720 * y - 240) + y * (360 * y - 360) + 60,
        x * (360 * x + 720 * y - 360)},
       diag},
      {diag,
       {y * (180 * y - 24), x * (360 * y - 24) + y * (360 * y - 288) + 24}}};

  // 18th basis function
  diag = {18 * M_SQRT2 * x * y - 9 * M_SQRT2 * y * y,
          x * (9 * M_SQRT2 * x - 18 * M_SQRT2 * y)};

  *(iter++) = {
      {{x * (-9 * M_SQRT2 * x + 18 * M_SQRT2 * y + 12 * M_SQRT2) - 3 * M_SQRT2,
        9 * M_SQRT2 * x * x},
       diag},
      {diag,
       {y * (-9 * M_SQRT2 * y + 12 * M_SQRT2),
        x * (-18 * M_SQRT2 * y + 12 * M_SQRT2) + 9 * M_SQRT2 * y * y -
            3 * M_SQRT2}}};

  // 19th basis function
  diag = {90 * M_SQRT2 * x * y + y * (45 * M_SQRT2 * y - 36 * M_SQRT2),
          x * (45 * M_SQRT2 * x + 90 * M_SQRT2 * y - 36 * M_SQRT2)};

  *(iter++) = {
      {{x * (-45 * M_SQRT2 * x - 90 * M_SQRT2 * y + 24 * M_SQRT2) + 3 * M_SQRT2,
        -45 * M_SQRT2 * x * x},
       diag},
      {diag,
       {y * (-45 * M_SQRT2 * y + 12 * M_SQRT2),
        x * (-90 * M_SQRT2 * y + 12 * M_SQRT2) +
            y * (-45 * M_SQRT2 * y + 36 * M_SQRT2) - 3 * M_SQRT2}}};

  // 20th basis function
  diag = {-36 * M_SQRT2 * x * y + 18 * M_SQRT2 * y * y,
          x * (-18 * M_SQRT2 * x + 36 * M_SQRT2 * y)};

  *(iter++) = {{{x * (18 * M_SQRT2 * x - 36 * M_SQRT2 * y - 12 * M_SQRT2) +
                     12 * M_SQRT2 * y,
                 x * (-18 * M_SQRT2 * x + 12 * M_SQRT2)},
                diag},
               {diag,
                {y * (18 * M_SQRT2 * y - 12 * M_SQRT2),
                 x * (36 * M_SQRT2 * y - 12 * M_SQRT2) +
                     y * (-18 * M_SQRT2 * y + 12 * M_SQRT2)}}};

  // 21th basis function
  diag = {-180 * M_SQRT2 * x * y + y * (-90 * M_SQRT2 * y + 72 * M_SQRT2),
          x * (-90 * M_SQRT2 * x - 180 * M_SQRT2 * y + 72 * M_SQRT2)};

  *(iter++) = {{{x * (90 * M_SQRT2 * x + 180 * M_SQRT2 * y - 60 * M_SQRT2) -
                     12 * M_SQRT2 * y,
                 x * (90 * M_SQRT2 * x - 12 * M_SQRT2)},
                diag},
               {diag,
                {y * (90 * M_SQRT2 * y - 12 * M_SQRT2),
                 x * (180 * M_SQRT2 * y - 12 * M_SQRT2) +
                     y * (90 * M_SQRT2 * y - 60 * M_SQRT2)}}};

  // 22th basis function
  diag = {0, 0};

  *(iter++) = {{{-48 * x - 24 * y + 24, -24 * x}, diag}, {diag, {0, 0}}};

  // 23th basis function
  diag = {24 * y, 24 * x};

  *(iter++) = {{{-48 * x - 48 * y + 24, -48 * x}, diag},
               {diag, {-48 * y, -48 * x - 48 * y + 24}}};

  // 24th basis function
  diag = {0, 0};

  *(iter++) = {{{0, 0}, diag}, {diag, {-24 * y, -24 * x - 48 * y + 24}}};
}
// generated with sympy from symfem library
template <class D, class R>
void ArnoldWintherReferenceLocalBasis<D, R>::partial(
    const std::array<unsigned int, dim> &order,
    const typename Traits::DomainType &in,
    std::vector<typename Traits::RangeType> &out) const {
  out.resize(size());
  auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
  if (totalOrder == 0) {
    evaluateFunction(in, out);
  } else if (totalOrder == 1) {
    if (order[0] == 1) {
      auto iter = out.begin();
      // generated with sympy from symfem library
      auto const &x = in[0], y = in[1];
      double diag;
      // 1th basis function
      diag = -12 * x * y + y * (6 - 9 * y);

      *(iter++) = {{x * (6 * x + 18 * y - 6), diag}, {diag, 6 * y * y}};

      // 2th basis function
      diag = x * (-30 * x - 90 * y + 36) + y * (48 - 45 * y) - 9;

      *(iter++) = {{x * (45 * x + 90 * y - 48) + y * (30 * y - 36) + 9, diag},
                   {diag, 60 * x * y + y * (45 * y - 36)}};

      // 3th basis function
      diag = -18 * x * y + y * (6 - 6 * y);

      *(iter++) = {{x * (9 * x + 12 * y - 6), diag},
                   {diag, x * (36 - 30 * x) + 9 * y * y - 9}};

      // 4th basis function
      diag = 12 * x * y + y * (9 * y - 6);

      *(iter++) = {{x * (-6 * x - 18 * y + 6), diag}, {diag, -6 * y * y}};

      // 5th basis function
      diag = x * (30 * x + 54 * y - 24) + y * (18 * y - 18) + 3;

      *(iter++) = {{x * (-27 * x - 36 * y + 18), diag},
                   {diag, -60 * x * y + y * (24 - 27 * y)}};

      // 6th basis function
      diag = 0;

      *(iter++) = {{0, diag}, {diag, x * (30 * x - 24) + 3}};

      // 7th basis function
      diag = 0;

      *(iter++) = {{0, diag}, {diag, 0}};

      // 8th basis function
      diag = 36 * x * y + y * (27 * y - 18);

      *(iter++) = {{x * (-18 * x - 54 * y + 18) + y * (24 - 30 * y) - 3, diag},
                   {diag, -18 * y * y}};

      // 9th basis function
      diag = 18 * x * y + y * (6 * y - 6);

      *(iter++) = {{x * (-9 * x - 12 * y + 6), diag}, {diag, -9 * y * y}};

      // 10th basis function
      diag = 108 * x * y + y * (36 * y - 36);

      *(iter++) = {{x * (-54 * x - 72 * y + 36), diag},
                   {diag, x * (180 * x - 192) + y * (24 - 54 * y) + 36}};

      // 11th basis function
      diag = x * (-180 * x - 360 * y + 192) + y * (132 - 90 * y) - 36;

      *(iter++) = {{x * (180 * x + 180 * y - 156) - 24 * y + 18, diag},
                   {diag, 360 * x * y + y * (180 * y - 192)}};

      // 12th basis function
      diag = -216 * x * y + y * (72 - 72 * y);

      *(iter++) = {{x * (108 * x + 144 * y - 72), diag},
                   {diag, x * (360 - 360 * x) + y * (108 * y - 48) - 60}};

      // 13th basis function
      diag = x * (360 * x + 720 * y - 360) + y * (180 * y - 240) + 60;

      *(iter++) = {{x * (-360 * x - 360 * y + 288) + 24 * y - 24, diag},
                   {diag, -720 * x * y + y * (360 - 360 * y)}};

      // 14th basis function
      diag = 72 * x * y + y * (54 * y - 36);

      *(iter++) = {{x * (-36 * x - 108 * y + 60) + 24 * y - 18, diag},
                   {diag, -36 * y * y}};

      // 15th basis function
      diag = 180 * x * y + y * (180 * y - 132);

      *(iter++) = {
          {x * (-90 * x - 360 * y + 132) + y * (192 - 180 * y) - 36, diag},
          {diag, y * (24 - 90 * y)}};

      // 16th basis function
      diag = -144 * x * y + y * (72 - 108 * y);

      *(iter++) = {{x * (72 * x + 216 * y - 96) - 48 * y + 24, diag},
                   {diag, 72 * y * y}};

      // 17th basis function
      diag = -360 * x * y + y * (240 - 360 * y);

      *(iter++) = {
          {x * (180 * x + 720 * y - 240) + y * (360 * y - 360) + 60, diag},
          {diag, y * (180 * y - 24)}};

      // 18th basis function
      diag = 18 * M_SQRT2 * x * y - 9 * M_SQRT2 * y * y;

      *(iter++) = {{x * (-9 * M_SQRT2 * x + 18 * M_SQRT2 * y + 12 * M_SQRT2) -
                        3 * M_SQRT2,
                    diag},
                   {diag, y * (-9 * M_SQRT2 * y + 12 * M_SQRT2)}};

      // 19th basis function
      diag = 90 * M_SQRT2 * x * y + y * (45 * M_SQRT2 * y - 36 * M_SQRT2);

      *(iter++) = {{x * (-45 * M_SQRT2 * x - 90 * M_SQRT2 * y + 24 * M_SQRT2) +
                        3 * M_SQRT2,
                    diag},
                   {diag, y * (-45 * M_SQRT2 * y + 12 * M_SQRT2)}};

      // 20th basis function
      diag = -36 * M_SQRT2 * x * y + 18 * M_SQRT2 * y * y;

      *(iter++) = {{x * (18 * M_SQRT2 * x - 36 * M_SQRT2 * y - 12 * M_SQRT2) +
                        12 * M_SQRT2 * y,
                    diag},
                   {diag, y * (18 * M_SQRT2 * y - 12 * M_SQRT2)}};

      // 21th basis function
      diag = -180 * M_SQRT2 * x * y + y * (-90 * M_SQRT2 * y + 72 * M_SQRT2);

      *(iter++) = {{x * (90 * M_SQRT2 * x + 180 * M_SQRT2 * y - 60 * M_SQRT2) -
                        12 * M_SQRT2 * y,
                    diag},
                   {diag, y * (90 * M_SQRT2 * y - 12 * M_SQRT2)}};

      // 22th basis function
      diag = 0;

      *(iter++) = {{-48 * x - 24 * y + 24, diag}, {diag, 0}};

      // 23th basis function
      diag = 24 * y;

      *(iter++) = {{-48 * x - 48 * y + 24, diag}, {diag, -48 * y}};

      // 24th basis function
      diag = 0;

      *(iter++) = {{0, diag}, {diag, -24 * y}};
    } else if (order[1] == 1) {
      auto iter = out.begin();
      // generated with sympy from symfem library
      auto const &x = in[0], y = in[1];
      double diag;
      // 1th basis function
      diag = x * (-6 * x - 18 * y + 6);

      *(iter++) = {{9 * x * x + y * (36 - 30 * y) - 9, diag},
                   {diag, 12 * x * y + y * (9 * y - 6)}};

      // 2th basis function
      diag = x * (-45 * x - 90 * y + 48) + y * (36 - 30 * y) - 9;

      *(iter++) = {{x * (45 * x + 60 * y - 36), diag},
                   {diag, x * (30 * x + 90 * y - 36) + y * (45 * y - 48) + 9}};

      // 3th basis function
      diag = x * (-9 * x - 12 * y + 6);

      *(iter++) = {{6 * x * x, diag}, {diag, 18 * x * y + y * (6 * y - 6)}};

      // 4th basis function
      diag = x * (6 * x + 18 * y - 6);

      *(iter++) = {{-9 * x * x, diag}, {diag, -12 * x * y + y * (6 - 9 * y)}};

      // 5th basis function
      diag = x * (27 * x + 36 * y - 18);

      *(iter++) = {{-18 * x * x, diag},
                   {diag, x * (-30 * x - 54 * y + 24) + y * (18 - 18 * y) - 3}};

      // 6th basis function
      diag = 0;

      *(iter++) = {{0, diag}, {diag, 0}};

      // 7th basis function
      diag = 0;

      *(iter++) = {{y * (30 * y - 24) + 3, diag}, {diag, 0}};

      // 8th basis function
      diag = x * (18 * x + 54 * y - 18) + y * (30 * y - 24) + 3;

      *(iter++) = {{x * (-27 * x - 60 * y + 24), diag},
                   {diag, -36 * x * y + y * (18 - 27 * y)}};

      // 9th basis function
      diag = x * (9 * x + 12 * y - 6);

      *(iter++) = {{-6 * x * x, diag}, {diag, -18 * x * y + y * (6 - 6 * y)}};

      // 10th basis function
      diag = x * (54 * x + 72 * y - 36);

      *(iter++) = {{-36 * x * x, diag},
                   {diag, x * (24 - 108 * y) + y * (60 - 36 * y) - 18}};

      // 11th basis function
      diag = x * (-180 * x - 180 * y + 132);

      *(iter++) = {
          {x * (90 * x - 24), diag},
          {diag, x * (180 * x + 360 * y - 192) + y * (90 * y - 132) + 36}};

      // 12th basis function
      diag = x * (-108 * x - 144 * y + 72);

      *(iter++) = {{72 * x * x, diag},
                   {diag, x * (216 * y - 48) + y * (72 * y - 96) + 24}};

      // 13th basis function
      diag = x * (360 * x + 360 * y - 240);

      *(iter++) = {
          {x * (24 - 180 * x), diag},
          {diag, x * (-360 * x - 720 * y + 360) + y * (240 - 180 * y) - 60}};

      // 14th basis function
      diag = x * (36 * x + 108 * y - 36);

      *(iter++) = {{x * (24 - 54 * x) + y * (180 * y - 192) + 36, diag},
                   {diag, -72 * x * y + y * (36 - 54 * y)}};

      // 15th basis function
      diag = x * (90 * x + 360 * y - 132) + y * (180 * y - 192) + 36;

      *(iter++) = {{x * (-180 * x - 360 * y + 192), diag},
                   {diag, x * (24 - 180 * y) + y * (156 - 180 * y) - 18}};

      // 16th basis function
      diag = x * (-72 * x - 216 * y + 72);

      *(iter++) = {{x * (108 * x - 48) + y * (360 - 360 * y) - 60, diag},
                   {diag, 144 * x * y + y * (108 * y - 72)}};

      // 17th basis function
      diag = x * (-180 * x - 720 * y + 240) + y * (360 - 360 * y) - 60;

      *(iter++) = {{x * (360 * x + 720 * y - 360), diag},
                   {diag, x * (360 * y - 24) + y * (360 * y - 288) + 24}};

      // 18th basis function
      diag = x * (9 * M_SQRT2 * x - 18 * M_SQRT2 * y);

      *(iter++) = {{9 * M_SQRT2 * x * x, diag},
                   {diag, x * (-18 * M_SQRT2 * y + 12 * M_SQRT2) +
                              9 * M_SQRT2 * y * y - 3 * M_SQRT2}};

      // 19th basis function
      diag = x * (45 * M_SQRT2 * x + 90 * M_SQRT2 * y - 36 * M_SQRT2);

      *(iter++) = {{-45 * M_SQRT2 * x * x, diag},
                   {diag, x * (-90 * M_SQRT2 * y + 12 * M_SQRT2) +
                              y * (-45 * M_SQRT2 * y + 36 * M_SQRT2) -
                              3 * M_SQRT2}};

      // 20th basis function
      diag = x * (-18 * M_SQRT2 * x + 36 * M_SQRT2 * y);

      *(iter++) = {{x * (-18 * M_SQRT2 * x + 12 * M_SQRT2), diag},
                   {diag, x * (36 * M_SQRT2 * y - 12 * M_SQRT2) +
                              y * (-18 * M_SQRT2 * y + 12 * M_SQRT2)}};

      // 21th basis function
      diag = x * (-90 * M_SQRT2 * x - 180 * M_SQRT2 * y + 72 * M_SQRT2);

      *(iter++) = {{x * (90 * M_SQRT2 * x - 12 * M_SQRT2), diag},
                   {diag, x * (180 * M_SQRT2 * y - 12 * M_SQRT2) +
                              y * (90 * M_SQRT2 * y - 60 * M_SQRT2)}};

      // 22th basis function
      diag = 0;

      *(iter++) = {{-24 * x, diag}, {diag, 0}};

      // 23th basis function
      diag = 24 * x;

      *(iter++) = {{-48 * x, diag}, {diag, -48 * x - 48 * y + 24}};

      // 24th basis function
      diag = 0;

      *(iter++) = {{0, diag}, {diag, -24 * x - 48 * y + 24}};
    } else
      DUNE_THROW(NotImplemented, "Invalid partial derivative");
  } else
    DUNE_THROW(NotImplemented, "Higher order derivatives are not implemented");
}


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
  void interpolate(const F &ff, std::vector<C> &out) const {
    auto &&f = Dune::Impl::makeFunctionWithCallOperator<LocalCoordinate>(ff);

    out.resize(size);
    auto it = out.begin();
    auto refElement = referenceElement(*element);

    // point evaluations
    // 9 DOFs in total
    for (auto i = 0u; i < element->subEntities(dim); ++i) {
      // auto value = f((*element).template subEntity<dim>(i).geometry().center());
      auto value = f(element->geometry().global(refElement.position(i,2)));

      it[0] = value[0][0];
      it[1] = value[0][1];
      it[2] = value[1][1];
      it += 3;
    }

    // integral moment over edges
    // 12 DOFs in total
    for (auto i = 0u; i < element->subEntities(1); ++i) {
      auto average = integralMoment<0>(f, (*element).template subEntity<1>(i).geometry(), quadratureOrder);
      auto firstMoment =
          integralMoment<1>(f, (*element).template subEntity<1>(i).geometry(), quadratureOrder + 1);

      std::size_t lower = (i == 2) ? 1 : 0;
      std::size_t upper = (i == 0) ? 1 : 2;
      auto tangent = (*element).template subEntity<dim>(upper).geometry().center() - (*element).template subEntity<dim>(lower).geometry().center();
      tangent /= tangent.two_norm();
      std::decay_t<decltype(tangent)> normal = {
          tangent[1], -tangent[0]}; // rotation by -90 degree

      using fRange = typename std::decay_t<std::remove_cv_t<decltype(f(std::declval<GlobalCoordinate>()))>>;
      using protomotedType =
          typename PromotionTraits<typename FieldTraits<fRange>::field_type,
                                   ctype>::PromotedType;

      FieldVector<protomotedType, 2> normalTimesAverage, normalTimesFirstMoment;
      average.mtv(normal, normalTimesAverage);
      firstMoment.mv(normal, normalTimesFirstMoment);
      it[0] = dot(normalTimesAverage, normal);
      it[1] = dot(normalTimesAverage, tangent);
      it[2] = dot(normalTimesFirstMoment, normal);
      it[3] = dot(normalTimesFirstMoment, tangent);
      // if edgeOrientation is wrong, then x -> 1-x therefore first moment -> zero moment - first moment
      // Note that because both normals and the tangent switches sign, these cancel
      // if (edgeOrientation_[i]){
      //   it[2] = it[0] - it[2];
      //   it[3] = it[1] - it[3];
      // }
      it += 4;
    }

    // integral moment on element
    // three DOFs in total
    auto average = integralMoment<0>(f, (*element).geometry(), quadratureOrder);
    it[0] = average[0][0];
    it[1] = average[0][1];
    it[2] = average[1][1];
  }

private:
  int quadratureOrder;
  std::bitset<3> edgeOrientation_;
  Element const* element = nullptr;

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
    auto integrationElement = geometry.integrationElement(xi);
    assert(values[0].N() == values[0].M());
    assert(values[0].N() == jacobian.N());
    assert(jacobian.M() == jacobian.N());

    for (auto &value : values) {
      value = jacobian * value * transpose(jacobian);
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

      referenceG[i] = {{tangent[1], tangent[0]}, {-tangent[0], tangent[1]}};
      auto tmp = tangent, tmp2 = tangent;
      jacobianTransposed.mtv(tangent, tmp);
      jacobianTransposed.mv(tmp, tmp2);
      referenceG[i].mtv(tmp2, tmp);
      // auto tmp = transpose(referenceG[i]) * jacobianTransposed * (transpose(jacobianTransposed) * tangent);
      alpha[i] = tmp[0] / jacobianDeterminant;
      beta[i] = tmp[1] / jacobianDeterminant;

      // already inverted W_k
      W_k[i] = 0;
      W_k[i][0][0] = 1.;
      W_k[i][1][0] = -alpha[i] / beta[i];
      W_k[i][1][1] = 1. / beta[i];
      // These entries transform the moments of order 1. They are orientation
      // dependent (due to the direction change in parametrization)
      if (!edgeOrientation_[i]) {
        W_k[i][2][2] = 1.;
        W_k[i][3][2] = -alpha[i] / beta[i];
        W_k[i][3][3] = 1. / beta[i];
      } else {
        // First order moment over a wronly oriented edge equals the zero order
        // moment minus the first order moment
        W_k[i][2][0] = 1.;
        W_k[i][2][2] = -1.;
        W_k[i][3][0] = -alpha[i] / beta[i];
        W_k[i][3][1] = 1. / beta[i];
        W_k[i][3][2] = alpha[i] / beta[i];
        W_k[i][3][3] = -1. / beta[i];
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
  /**
   * \brief Apply the transformation to some Vector of Shapevalues, Jacobians or
   * Hessians
   *
   * \tparam Values Vector
   * \param values
   */
  template <class Values> void apply(Values &values) const {
    Values tmp = values;
    mat_.mv(tmp, values);
  }

  /**
   * \brief Apply the Inverse transformation to some Vector of Shapevalues,
   * Jacobians or Hessians
   *
   * \tparam Values Vector
   * \param values
   */
  template <class Values> void applyInverse(Values &values) const {
    Values tmp = values;
    mat_.getInverse().mv(tmp, values);
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
class ArnoldWintherPreBasis : public LeafPreBasisMapperMixin<GV> {
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
      : Base(gv, arnoldWintherMapperLayout),
        mapper_({gv, mcmgElementLayout()}) {
    updateState(gv);
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update(GridView const &gv) {
    Base::update(gv);
    updateState(gv);
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const { return Node{mapper_, data_}; }

protected:
  // \TODO replace with external functionality once argyris is merged
  void updateState(GridView const &gridView)
  {
    data_.resize(gridView.size(0));
    // compute orientation for all elements
    std::bitset<3> orientation = 0;
    auto const &idSet = gridView.grid().globalIdSet();
    auto const& indexSet = gridView.indexSet();

    for (const auto &element : elements(gridView)) {
      const auto &refElement = referenceElement(element);
      auto elementIndex = mapper_.index(element);

      orientation = 0;

      for (std::size_t i = 0; i < element.subEntities(dim - 1); i++) {
        // Local vertex indices within the element
        auto localV0 = refElement.subEntity(i, dim - 1, 0, dim);
        auto localV1 = refElement.subEntity(i, dim - 1, 1, dim);

        // Global vertex indices within the grid
        auto globalV0 = idSet.subId(element, localV0, dim);
        auto globalV1 = idSet.subId(element, localV1, dim);
        // The edge is flipped if the local ordering disagrees with global
        // ordering
        if(indexSet.subIndex(element.template subEntity<1>(i),0,2)>indexSet.subIndex(element.template subEntity<1>(i),1,2)){
        // if ((localV0 < localV1 && globalV0 > globalV1) ||
            // (localV0 > localV1 && globalV0 < globalV1)) {
          orientation[i] = 1;
        }
      }
      data_[elementIndex] = orientation;
    }
  }

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
    this->setSize(finiteElement_.size());
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
#endif