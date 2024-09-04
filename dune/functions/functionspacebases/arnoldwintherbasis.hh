// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ARNOLDWINTHERBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ARNOLDWINTHERBASIS_HH

#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/scalarvectorview.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/functions/functionspacebases/arnoldwinther.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/doublepiola.hh>
#include <dune/functions/functionspacebases/lineartransformedlocalfiniteelement.hh>
#include <dune/functions/functionspacebases/vectorfloatcmp.hh>

#include <dune/localfunctions/tensormatvec.hh>

/**
 * \brief Implementation of the conforming Arnold-Winther element
 * This is a finite element used to discretize the stress for two-dimensional
 * elasticity, originally proposed in "Arnold, D. N., & Winther, R. (2002).
 * Mixed finite elements for elasticity." As such, its shapefunctions take
 values
 * in the space of symmetric 2x2 matrices, whose (rowwise) divergence is in L2.

 * This comes with some complications in the Dune framework, in particular, there is
 * no datastructure for 3-tensors, which is the JacobianType of this finite
 element.
 * Hence, we must add certain functionality, implemented directly below, which
 * might have sideeffects. Use this file with caution.
 *   The implementation, in particular the transformation accounting for the
 * non-affinity is based on "Aznaran, Francis & Kirby, Robert & Farrell, Patrick.
 * (2021). Transformations for Piola-mapped elements."
 */
// TODO Rework this to incorporate some stuff for nonaffine mappings
// In fact, I don't think this element is well defined for nonaffine mappings
namespace Dune {

// We need to stretch the FMatrix / FVector interface a bit to allow the
// necessary tensor operations. This could be replaced by introducing a
// dependency on dune-tensor
template <class K, class L, int n, int m>
auto operator*(FieldMatrix<K, n, m> const &mat, FieldVector<L, m> vec) {
  FieldVector<typename PromotionTraits<K, L>::PromotedType, n> result;
  mat.mv(vec, result);
  return result;
}

/** \brief Compute type of the result of an arithmetic operation involving two
 * different number types.
 */
template <typename T1, typename T2, int n>
struct PromotionTraits<FieldVector<T1, n>, T2> {
  typedef FieldVector<typename PromotionTraits<T1, T2>::PromotedType, n>
      PromotedType;
};

template <typename T1, typename T2, int n>
struct PromotionTraits<T2, FieldVector<T1, n>> {
  typedef FieldVector<typename PromotionTraits<T1, T2>::PromotedType, n>
      PromotedType;
};

namespace Functions {
namespace Impl {

/**
 * \brief Class offering reading and writing access to a Vector starting from an index
 * \note This does not implement a vector interface, but only some minimal features used in the code below.
 */
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

  VectorSlice& operator=(value_type scalar) {
    for (size_type i = 0u; i < size(); ++i)
      vec_[i_ + i] = scalar;
    return *this;
  }

  value_type const& operator[](size_type const &index) const {
    return vec_[i_ + index];
  }

  value_type& operator[](size_type const &index) { return vec_[i_ + index]; }

  size_type N() const { return end_ - i_; }

  size_type size() const { return end_ - i_; }

private:
  Vector& vec_;
  size_type i_;
  size_type end_;
};

template <class T>
struct FieldTraits<typename Functions::Impl::VectorSlice<T>> {
  typedef typename FieldTraits<T>::field_type field_type;
  typedef typename FieldTraits<T>::real_type real_type;
};

//  TODO make this fulfill Dune interfaces
// TODO make generic in matrix Types and sizes
/**
 * \brief Block Diagonal Matrix with hardcoded dimensions that fit the Arnold
 * Winther fe transformation. It models the $P$ Matrix in the above mentioned
 * paper. The only operation needed for this purpose is the multiplication of
 * its transpose with a vector and access to its inverse. This Vector however,
 * has values which are Matrices or order 3-tensors for the jacobian case (for now, Matrices of Vectors).
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

  // implementation that takes 3-tensors as arguments. Unfortunately, this does
  // not work by default
  template <class KIn, class KOut, int dim>
  void mv(std::vector<FieldMatrix<FieldVector<KIn, dim>, dim, dim>> const &x,
          std::vector<FieldMatrix<FieldVector<KOut, dim>, dim, dim>> &y) const {
    auto &&xx = Dune::Impl::asVector(x);
    auto &&yy = Dune::Impl::asVector(y);
    DUNE_ASSERT_BOUNDS((void *)(&x) != (void *)(&y));
    DUNE_ASSERT_BOUNDS(xx.N() == M());
    DUNE_ASSERT_BOUNDS(yy.N() == N());
    size_type index = 0;
    for (auto const &mat : transformPointDofs_) {
      for (auto i = 0u; i < mat.N(); ++i)
        for (auto k = 0u; k < y[index + i].N(); ++k)
          for (auto l = 0u; l < y[index + i].M(); ++l) {
            y[index + i][k][l] = 0;
            for (auto j = 0u; j < mat.M(); ++j)
              y[index + i][k][l] += mat[i][j] * x[index + j][k][l];
          }

      index += mat.M();
    }
    for (auto const &mat : transformEdgeDofs_) {
      for (auto i = 0u; i < mat.N(); ++i)
        for (auto k = 0u; k < y[index + i].N(); ++k)
          for (auto l = 0u; l < y[index + i].M(); ++l) {
            y[index + i][k][l] = 0;
            for (auto j = 0u; j < mat.M(); ++j)
              y[index + i][k][l] += mat[i][j] * x[index + j][k][l];
          }
      index += mat.M();
    }

    for (auto i = 0u; i < transformElementDofs_.N(); ++i)
      for (auto k = 0u; k < y[index + i].N(); ++k)
        for (auto l = 0u; l < y[index + i].M(); ++l) {
          y[index + i][k][l] = 0;
          for (auto j = 0u; j < transformElementDofs_.M(); ++j)
            y[index + i][k][l] +=
                transformElementDofs_[i][j] * x[index + j][k][l];
        }
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

  // implementation that takes 3-tensors as arguments. Unfortunately, this does
  // not work by default
  template <class KIn, class KOut, int dim>
  void
  mtv(std::vector<FieldMatrix<FieldVector<KIn, dim>, dim, dim>> const &x,
      std::vector<FieldMatrix<FieldVector<KOut, dim>, dim, dim>> &y) const {
    auto &&xx = Dune::Impl::asVector(x);
    auto &&yy = Dune::Impl::asVector(y);
    DUNE_ASSERT_BOUNDS((void *)(&x) != (void *)(&y));
    DUNE_ASSERT_BOUNDS(xx.N() == M());
    DUNE_ASSERT_BOUNDS(yy.N() == N());
    size_type index = 0;
    for (auto const &mat : transformPointDofs_) {
      for (auto i = 0u; i < mat.N(); ++i)
        for (auto k = 0u; k < y[index + i].N(); ++k)
          for (auto l = 0u; l < y[index + i].M(); ++l) {
            y[index + i][k][l] = 0;
            for (auto j = 0u; j < mat.M(); ++j)
              y[index + i][k][l] += mat[j][i] * x[index + j][k][l];
          }

      index += mat.M();
    }
    for (auto const &mat : transformEdgeDofs_) {
      for (auto i = 0u; i < mat.N(); ++i)
        for (auto k = 0u; k < y[index + i].N(); ++k)
          for (auto l = 0u; l < y[index + i].M(); ++l) {
            y[index + i][k][l] = 0;
            for (auto j = 0u; j < mat.M(); ++j)
              y[index + i][k][l] += mat[j][i] * x[index + j][k][l];
          }
      index += mat.M();
    }

    for (auto i = 0u; i < transformElementDofs_.N(); ++i)
      for (auto k = 0u; k < y[index + i].N(); ++k)
        for (auto l = 0u; l < y[index + i].M(); ++l) {
          y[index + i][k][l] = 0;
          for (auto j = 0u; j < transformElementDofs_.M(); ++j)
            y[index + i][k][l] +=
                transformElementDofs_[j][i] * x[index + j][k][l];
        }
  }

  static constexpr unsigned int N() { return 24; }
  static constexpr unsigned int M() { return 24; }
  static constexpr unsigned int rows() { return 24; }
  static constexpr unsigned int cols() { return 24; }

  This getInverse() const {
    This tmp = *this;
    tmp.invert();
    return tmp;
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

template<class GV>
using EdgeOrientationGlobalState = std::tuple<Dune::MultipleCodimMultipleGeomTypeMapper<GV>, std::vector<bool>;

/**
 * \brief compute EdgeOrientation for each edge in the GridView. Currently, for the order we use the vertexindex
 * in the sequential case and a global Lexicographic coordinate comparison for distributed grids.
 */
template<class GV>
void computeGlobalState(EdgeOrientationGlobalState<GV>& state, const GV &gv) {

    auto& [elementMapper_, elementInformation_ ] = state;
    elementMapper_.update(gv);
    elementInformation_.resize(elementMapper_.size(1));
    // compute orientation for all elements
    unsigned short orientation;
    bool sequentialSetup = (gv.comm().size() == 1);
    auto const &indexSet = gv.indexSet();
    for (const auto &element : elements(gv)) {
      const auto &refElement = referenceElement(element);
      auto elementIndex = elementMapper_.index(element);

      for (std::size_t i = 0; i < element.subEntities(dim - 1); i++) {
        // Local vertex indices within the element
        auto localV0 = refElement.subEntity(i, dim - 1, 0, dim);
        auto localV1 = refElement.subEntity(i, dim - 1, 1, dim);
        if (sequentialSetup) {
          // Global vertex indices within the grid
          auto globalV0 = indexSet.subIndex(element, localV0, dim);
          auto globalV1 = indexSet.subIndex(element, localV1, dim);

          if ((localV0 < localV1 && globalV0 > globalV1) ||
              (localV0 > localV1 && globalV0 < globalV1))
            elementInformation_[elementMapper_.subIndex(element, i, dim-1)] |= (1 << i);
        }
        else {

          // sort lexicographically by coordinate
          // this ensures consistent orientation also for distributed grids

          auto globalV0 =
              element.template subEntity<dim>(localV0).geometry().corner(0);
          auto globalV1 =
              element.template subEntity<dim>(localV1).geometry().corner(0);

          if ((localV0 < localV1 && vectorGreater(globalV0, globalV1)) ||
              (localV0 > localV1 && vectorLess(globalV0, globalV1)))
            elementInformation_[elementMapper_.subIndex(element, i, dim-1)] |= (1 << i);

        }
      }
    }
  }


/** \brief Linear transformation that maps the reference basis onto the
 * pull-backs of physical nodal basis for the ArnoldWinther element
 * \tparam R
 * RangeFieldType of finite element
 */
template <class R, class GV>
class ArnoldWintherTransformator {
public:
  static constexpr int dim = 2;
  /** The Globalstate we use to find the local State for each element*/
  using GlobalState = typename EdgeOrientationGlobalState<GV>;
  /**
   * \brief binds the transformation to an element and its elementinformation
   *        Fills the transformation Matrix.
   *
   * \tparam Element
   * \param element
   */
  template <class Element>
  void bind(Element const &element) {
    element_ = &element
    fillMatrix(
        Dune::referenceElement<double, 2>(GeometryTypes::simplex(2)).position(0, 0),// barycenter, because we need some value.
        element_.geometry());
  } // TODO actually we don't and it would be more correct to use the jacobians at the dof locations

  /**
   * \brief Apply the transformation to some Vector of Shapevalues, Jacobians or
   * Hessians
   *
   * \tparam Values Vector
   * \param values
   */
  template <class Values>
  void apply(Values &values) const {
    Values tmp = values;
    mat_.mtv(tmp, values);
  }

  /**
   * \brief Apply the Inverse transformation to some Vector of Shapevalues,
   * Jacobians or Hessians
   *
   * \tparam Values Vector
   * \param values
   */
  template <class Values>
  void applyInverse(Values &values) const {
    Values tmp = values;
    mat_.getInverse().mv(tmp, values);
  }

  LocalState const& localState() const{
    return localState_;
  }

private:
  template <class LocalCoordinate, class Geometry>
  void fillMatrix(LocalCoordinate const &x, Geometry const &geometry) {
    std::array<R, 3> alpha;
    std::array<R, 3> beta;

    std::array<Dune::FieldVector<R, 2>, 3> referenceTangents; // normalized

    std::array<R, 3> referenceEdgeLength;
    std::array<R, 3> globalEdgeLength;

    std::array<Dune::FieldMatrix<R, 2, 2>, 3> referenceG;

    std::array<Dune::FieldMatrix<R, 4, 4>, 3> W_k;
    Dune::FieldMatrix<R, 3, 3> W;

    auto jacobianTransposed = geometry.jacobianTransposed(x);
    auto jacobianDeterminant = geometry.integrationElement(x);

    // By default, edges point from the vertex with the smaller index
    // to the vertex with the larger index. Note that the alpha and beta are invariant
    // of orientation, since the normals/tangents appear twice in their definitions.

    // get local and global Tangents
    auto refElement = Dune::referenceElement<double, 2>(geometry.type());
    for (std::size_t i = 0; i < 3; ++i) {
      std::size_t lower = (i == 2) ? 1 : 0;
      std::size_t upper = (i == 0) ? 1 : 2;
      auto tangent = refElement.position(upper, 2) - refElement.position(lower, 2);

      referenceEdgeLength[i] = tangent.two_norm();
      tangent/= referenceEdgeLength[i];

      auto globalEdge = geometry.global(refElement.position(upper, 2)) -
                        geometry.global(refElement.position(lower, 2));

      globalEdgeLength[i] = globalEdge.two_norm();

      referenceG[i] = {{tangent[1], tangent[0]}, {-tangent[0], tangent[1]}};
      auto tmp = transpose(referenceG[i]) * jacobianTransposed *
                 transpose(jacobianTransposed) * tangent;
      alpha[i] = tmp[0] / jacobianDeterminant;
      beta[i] = tmp[1] / jacobianDeterminant;

      // already inverted W_k
      W_k[i] = 0;
      W_k[i][0][0] = 1.;
      W_k[i][1][0] = -alpha[i] / beta[i];
      W_k[i][1][1] = 1. / beta[i];
      // These entries transform the moments of order 1. They are orientation
      // dependent (due to the direction change in paramtrization)
      if (!std::apply([&e = *element_](auto &&mapper, auto &&data) { return data[mapper.subIndex(e, i, dim - 1)]; },*globalState_))
      {
        W_k[i][2][2] = 1.;
        W_k[i][3][2] = -alpha[i] / beta[i];
        W_k[i][3][3] = 1. / beta[i];
      } else {
        // First order moment over a wronly oriented edge equals the zero order
        // moment minus the first order moment
        W_k[i][2][0] = 1.;
        W_k[i][2][2] = -1.;
        W_k[i][3][0] = -alpha[i] / beta[i];
        W_k[i][3][1] = 1./ beta[i];
        W_k[i][3][2] = alpha[i] / beta[i];
        W_k[i][3][3] = -1. / beta[i];
      }
      W_k[i] *= globalEdgeLength[i] / referenceEdgeLength[i];
    }
    // Fill W (not yet inverted)
    // TODO this should be improved to handle DiagonalMatrices as well. Since we
    // only have simplices, I think this case currently cannot arise tho.

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
    W *= jacobianDeterminant * jacobianDeterminant;
    // fill matrix
    mat_ = ArnoldWintherBlockDiagonalMatrix<R>{
        std::array<Dune::FieldMatrix<R, 3, 3>, 3>{W, W, W}, W_k,
        W / jacobianDeterminant};
  }

  // Blockmatrix This is the matrix P from the paper mentioned above
  ArnoldWintherBlockDiagonalMatrix<R> mat_;
  GlobalState const *globalState_;
  Element const* element_;
}; // Transformator


  template<class Element>
  struct DoubleContravariantPiolaTransformator {

    void bind(Element const& e){
      element_ = &e;
    }

    /** \brief Double Piola-transform a set of shape-function values
      *
      * \param[in,out] values The values to be Piola-transformed
      */
    template <typename Values, typename LocalCoordinate, typename Geometry>
    static auto apply(Values &values) {
      auto xi = Dune::referenceElement<double, 2>(GeometryTypes::simplex(2)).position(0, 0);
      auto const& geometry = e->geometry();
      auto jacobian = geometry.jacobian(xi);
      auto integrationElement = geometry.integrationElement(xi);
      assert(values[0].N() == values[0].M());
      assert(values[0].N() == jacobian.N());
      assert(jacobian.M() == jacobian.N());

      // value = jacobian * value * transpose(jacobian);
      for (auto &value : values) {
        auto tmp = value;
        value = 0;
        for (std::size_t k = 0; k < jacobian.N(); k++)
          for (std::size_t l = 0; l < jacobian.N(); l++)
            for (auto&& [jacobian_k_i, i] : sparseRange(jacobian[k]))
              for (auto&& [jacobian_l_j, j] : sparseRange(jacobian[l]))
                value[k][l] += jacobian_k_i * tmp[i][j] * jacobian_l_j;

        value = value/(integrationElement * integrationElement);
      }
      return;
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
            auto value = f_(xi);

            // Apply the inverse Piola transform
            auto jacobian = element_.geometry().jacobianInverse(xi);
            auto integrationElement =
                element_.geometry().integrationElement(xi);

            auto tmp = value;
            value = 0;
            for (std::size_t k = 0; k < jacobian.N(); k++)
              for (std::size_t l = 0; l < jacobian.N(); l++)
                for (auto&& [jacobian_k_i, i] : sparseRange(jacobian[k]))
                  for (auto&& [jacobian_l_j, j] : sparseRange(jacobian[l]))
                    value[k][l] += jacobian_k_i * tmp[i][j] * jacobian_l_j;

            value = value/(integrationElement * integrationElement);

            return value;
          }
      };
  };

} // namespace Impl

template <class GV, class R>
class ArnoldWintherNode;

template <class GV, typename R>
class ArnoldWintherPreBasis {
  static const int dim = GV::dimension;
  static_assert(dim == 2,
                "ArnoldWinther PreBasis only implemented for 2d simplices");

public:
  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  using Node = ArnoldWintherNode<GridView, R>;

  static constexpr size_type maxMultiIndexSize = 1;
  static constexpr size_type minMultiIndexSize = 1;
  static constexpr size_type multiIndexBufferSize = 1;

  //! Constructor for a given grid view object
  ArnoldWintherPreBasis(const GV &gv)
      : gridView_(gv) {
        computeGlobalState(globalState_, gridView_);
      }

  //! Initialize the global indices
  void initializeIndices() {}

  //! Obtain the grid view that the basis is defined on
  const GridView &gridView() const { return gridView_; }

  //! Update the stored grid view, to be called if the grid has changed
  void update(const GridView &gv) {
    gridView_ = gv;
    computeGlobalState(globalState_, gridView_);
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const { return Node{&globalState_}; }

  //! Same as size(prefix) with empty prefix
  size_type size() const {
    return 3 * gridView_.size(2) + 4 * gridView_.size(1) +
           3 * gridView_.size(0);
  }

  //! Return number of possible values for next position in multi index
  template <class SizePrefix>
  size_type size(const SizePrefix prefix) const {
    // this basically means this is a leaf node with dimrange 1 right?
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? size() : 0;
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const { return size(); }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const { return 24; }

  template <typename It>
  It indices(const Node &node, It it) const {
    const auto &gridIndexSet = gridView().indexSet();
    const auto &element = node.element();

    // throw if Element is not simplex
    if (not(element.type().isSimplex()))
      DUNE_THROW(Dune::NotImplemented,
                 "ArnoldWinther Basis only implemented for simplex elements");
    for (size_type i = 0, end = node.finiteElement().size(); i < end;
         ++it, ++i) {
      Dune::LocalKey localKey =
          node.finiteElement().localCoefficients().localKey(i);
      // TODO probably unnecessary
      if (!gridIndexSet.contains(element))
        DUNE_THROW(Dune::RangeError, "Element is not in gridIndexSet!");
      if (localKey.codim() == 0)
        *it = {{(size_type)3 *
                    gridIndexSet.subIndex(element, localKey.subEntity(), 0) +
                localKey.index()}};
      else if (localKey.codim() == 1)
        *it = {{(size_type)3 * gridView().size(0) +
                4 * gridIndexSet.subIndex(element, localKey.subEntity(), 1) +
                localKey.index()}};
      else if (localKey.codim() == 2)
        *it = {{(size_type)3 * gridView().size(0) + 4 * gridView().size(1) +
                3 * gridIndexSet.subIndex(element, localKey.subEntity(), 2) +
                localKey.index()}};
    }
    return it;
  }

protected:
  GridView gridView_;
  Impl::EdgeOrientationGlobalState globalState_;
};

template <class GV, class R>
class ArnoldWintherNode : public LeafBasisNode {
public:
  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;

private:
  static constexpr unsigned int dim = GV::dimension;
  using LocalValuedFE = ArnoldWintherLocalFiniteElement<typename GV::ctype, R>;
  using GlobalValuedFE = Impl::GlobalValuedLocalFiniteElement<
      Impl::DoubleContravariantPiolaTransformator, LocalValuedFE, Element>;
  using GlobalState = Impl::EdgeOrientationGlobalState<GV>;

public:
  using FiniteElement = Impl::LinearTransformedLocalFiniteElement<
      Impl::ArnoldWintherTransformator<R>, GlobalValuedFE, Element>;

  ArnoldWintherNode(GlobalState const& state)
      : globalState_(state)
      , finiteElement_(state) {
    this->setSize(finiteElement_->size());
  }

  ~ArnoldWintherNode(){}

  //! Return current element, throw if unbound
  const Element &element() const { return element_; }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the
   * dune-localfunctions module
   */
  const FiniteElement &finiteElement() const { return *finiteElement_; }

  //! Bind to element.
  void bind(const Element &e) {
    if (not e.type().isSimplex())
      DUNE_THROW(Dune::NotImplemented,
                 "ArnoldWintherBasis can only be bound to simplex elements");
    element_ = e;
    finiteElement_.bind(element_);
  }

  unsigned int order() const { return 3; }

protected:
  const GlobalState *globalState;
  FiniteElement finiteElement_;
  Element element_;
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

template <typename Range = double>
auto arnoldWinther() {
  return [](auto const &gridView) {
    return ArnoldWintherPreBasis<std::decay_t<decltype(gridView)>, Range>(
        gridView);
  };
}
} // namespace BasisFactory

} // namespace Functions


} // namespace Dune

#endif
