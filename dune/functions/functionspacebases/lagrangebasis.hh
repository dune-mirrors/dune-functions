// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH

#include <type_traits>
#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/leafprebasismixin.hh>


namespace Dune {
namespace Functions {

// *****************************************************************************
// This is the reusable part of the LagrangeBasis. It contains
//
//   LagrangePreBasis
//   LagrangeNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename R=double>
class LagrangeNode;

template<typename GV, int k, typename R=double>
class LagrangePreBasis;



/**
 * \brief A pre-basis for a PQ-lagrange bases with given order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam k   The polynomial order of ansatz functions; -1 means 'order determined at run-time'
 * \tparam R   Range type used for shape function values
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * \warning For pyramid elements in 3d, the shape functions are different for
 *    run-time and compile-time order. While the former are defined using the
 *    Duffy-transformation the latter are continuous and piecewise polynomial
 *    with discontinuous gradients.
 */
template<typename GV, int k, typename R>
class LagrangePreBasis :
  public LeafPreBasisMixin< LagrangePreBasis<GV,k,R> >
{
  static const int dim = GV::dimension;
  static const bool useDynamicOrder = (k<0);

  static Dune::MCMGLayout dofLayout(unsigned int order)
  {
    return [order](Dune::GeometryType type, int dim) -> std::size_t {
      if (order==0)
        return type.dim() == dim;
      if (type.isSimplex())
        return Dune::binomial(order-1, type.dim());
      if (type.isCube())
        return Dune::power(order-1, type.dim());
      if (type.isPyramid() and (order>2))
        return (order-2)*(order-1)*(2*order-3)/6;
      if (type.isPrism() and (order>1))
        return (order-1)*(order-1)*(order-2)/2;
      return 0;
    };
  }

  // The following methods compute orientations of line/triangle/quadrilateral
  // in the following sense: Given a range of global indices (or IDs) of the vertices
  // the computed orientation is a number that indicates how the geometry has
  // to be reoriented to match the reference orientation.
  //
  // For a 1d line this is a simple bit indicating if it has to be flipped
  // corresponding to a reflection along the x-axis.
  // For a 2d triangle/quadrilateral the bits indicate the refelctions needed
  // to get the desired orientation:
  // Bit 0: Reflect along the x-axis
  // Bit 1: Reflect along the y-axis
  // Bit 2: Reflect across the diagonal
  // Notice that the order of these operations does matter and they
  // have to be executed in order with the significance of the bit.

  template<class VertexIndices>
  auto computeLineOrientation(const VertexIndices& vertexIndices) const
  {
    return vertexIndices[0] > vertexIndices[1];
  }

  template<class VertexIndices>
  auto computeTriangleOrientation(const VertexIndices& vertexIndices) const
  {
    auto flipVertices = std::bitset<3>();
    flipVertices[0] = vertexIndices[0] > vertexIndices[1];
    if (flipVertices[0])
    {
      flipVertices[1] = vertexIndices[1] > vertexIndices[2];
      flipVertices[2] = vertexIndices[0] > vertexIndices[2];
    }
    else
    {
      flipVertices[1] = vertexIndices[0] > vertexIndices[2];
      flipVertices[2] = vertexIndices[1] > vertexIndices[2];
    }
    return flipVertices.to_ulong();
  }

  template<class VertexIndices>
  auto computeQuadrilateralOrientation(const VertexIndices& vertexIndices) const
  {
    std::size_t i_min = 0;
    for(auto i: Dune::range(1, 4))
      if (vertexIndices[i] < vertexIndices[i_min])
        i_min = i;
    auto flipVertices = std::bitset<3>(i_min);
    if (i_min==0)
      flipVertices[2] = vertexIndices[1] > vertexIndices[2];
    else if (i_min==1)
      flipVertices[2] = vertexIndices[0] > vertexIndices[3];
    else if (i_min==2)
      flipVertices[2] = vertexIndices[0] < vertexIndices[3];
    else if (i_min==3)
      flipVertices[2] = vertexIndices[1] < vertexIndices[2];
    return flipVertices.to_ulong();
  }

  // The following methods flip DOFs indices according to the orientation
  // as computed by the preceding methods. For a line we simple need
  // to reverse the order-1 DOFs. For triangle/quadrilateral we us
  // precomputed tables.
  //
  // Given a point p in the subentity identified by its local index i
  // wrt the reference element orientation, these methods compute the
  // local index j of this point with respect to the global orientation.
  // To find j we can apply the global to local transformation to p.
  // Then the local index of the obtained point wrt the reference element
  // orientation is j.
  auto flipLineDOF(auto index, auto orientation) const
  {
    if (orientation==1)
      index = order()-2-index;
    return index;
  }

  auto flipTriangleDOF(auto index, auto orientation) const
  {
    return triangleFlipTable_[orientation][index];
  }

  auto flipQuadrilateralDOF(auto index, auto orientation) const
  {
    return quadrilateralFlipTable_[orientation][index];
  }

  // The following methods compute flip tables for DOFs associated to a triangle/quadrilateral.
  // Given an orientation as computed by the preceding methods, the corresponding row of
  // the table gives the renumbering of the DOFs associated with the subentity such
  // that the orientation matches with the global indices and is thus unique.

  static auto computeTriangleFlipTable(int order)
  {
    auto flipTable = std::array<std::vector<std::size_t>, 8>();
    if (order<3)
      return flipTable;
    // Number of DOFs on the base line
    std::size_t m = order-2;
    // Total number of DOFs
    std::size_t n = m*(m+1)/2;
    // In order to reflect we need to flip pairs of vertices.
    // To do this we first map flat DOF indices to barycentric
    // multi-indices wrt. the vertices and summing up to m-1.
    // Then the DOF reorientation of a vertex flip simply
    // requires to flip the corresponding components of the
    // barycentric multi-index.
    using BarycentricIndex = std::array<std::size_t, 3>;
    auto barycentricIndices = std::vector<BarycentricIndex>();
    for(auto i : Dune::range(m))
      for(auto j : Dune::range(m-i))
        barycentricIndices.push_back(BarycentricIndex{m-1-(i+j), j, i});
    auto flatToBarycentric = [&](auto i) { return barycentricIndices[i]; };
    auto barycentricToFlat = [&](auto b) { return b[1] + m*b[2] - b[2]*(b[2]-1)/2; };
    // Loop of all 8 orinetations
    for(auto j : Dune::range(8))
    {
      flipTable[j].resize(n);
      for(auto i : Dune::range(n))
      {
        auto b = flatToBarycentric(i);
        // Bit 0: Reflect along x-axis -> flip vertex 0 and 1
        if (j & (1 << 0))
          std::swap(b[0], b[1]);
        // Bit 1: Reflect along y-axis -> flip vertex 0 and 2
        if (j & (1 << 1))
          std::swap(b[0], b[2]);
        // Bit 2: Reflect across diagonal -> flip vertex 1 and 2
        if (j & (1 << 2))
          std::swap(b[1], b[2]);
        flipTable[j][i] = barycentricToFlat(b);
      }
    }
    return flipTable;
  }

  static auto computeQuadrilateralFlipTable(int order)
  {
    auto flipTable = std::array<std::vector<std::size_t>, 8>();
    if (order<2)
      return flipTable;
    // Number of DOFs on the base line
    std::size_t m = order-1;
    // Total number of DOFs
    std::size_t n = m*m;
    // In order to reflect we map the flat DOF index into a cartesian
    // multi-index wrt. the x- and y-axis.
    using CartesianIndex = std::array<std::size_t, 3>;
    auto flatToCartesian = [&](auto i) { return CartesianIndex{i % m, i / m}; };
    auto cartesianToFlat = [&](auto c) { return c[0] + m*c[1]; };
    for(auto j : Dune::range(8))
    {
      flipTable[j].resize(n);
      for(auto i : Dune::range(n))
      {
        auto c = flatToCartesian(i);
        // Bit 0: Reflect along x-axis -> reverse x-index
        if (j & (1 << 0))
          c[0] = m-1-c[0];
        // Bit 1: Reflect along y-axis -> reverse y-index
        if (j & (1 << 1))
          c[1] = m-1-c[1];
        // Bit 2: Reflect across diagonal -> flip x- and y-index
        if (j & (1 << 2))
          std::swap(c[0], c[1]);
        flipTable[j][i] = cartesianToFlat(c);
      }
    }
    return flipTable;
  }


public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  using Node = LagrangeNode<GV, k, R>;

  //! Constructor for a given grid view object with compile-time order
  LagrangePreBasis(const GridView& gv)
  : LagrangePreBasis(gv, std::numeric_limits<unsigned int>::max())
  {}

  //! Constructor for a given grid view object and run-time order
  LagrangePreBasis(const GridView& gv, unsigned int order) :
    gridView_(gv), order_(order)
    , mapper_(gridView_, dofLayout(useDynamicOrder ? order : k))
  {
    if (!useDynamicOrder && order!=std::numeric_limits<unsigned int>::max())
      DUNE_THROW(RangeError, "Template argument k has to be -1 when supplying a run-time order!");
  }

  //! Initialize the global indices
  void initializeIndices()
  {
    triangleFlipTable_ = computeTriangleFlipTable(order());
    quadrilateralFlipTable_ = computeQuadrilateralFlipTable(order());
  }

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
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{order_};
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return mapper_.size();
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    // That cast to unsigned int is necessary because GV::dimension is an enum,
    // which is not recognized by the power method as an integer type...
    return power(order()+1, (unsigned int)GV::dimension);
  }

//  auto computeFaceOrientations(const typename Node::Element& element) const
//  {
//    const auto& indexSet = gridView().indexSet();
//    const auto& re = Dune::referenceElement<double,dim>(element.type());

//    auto globalVertexIndicesOfFace = [&](auto faceIndex, auto faceCodim) {
//      return Dune::transformedRangeView(re.subEntities(faceIndex, faceCodim, dim), [&](auto localVertexIndex) {
//        return indexSet.subIndex(element, localVertexIndex, dim);
//      });
//    };

//    if constexpr (dim==2)
//    {
//      auto orientations = std::array<std::array<std::size_t, 6>, 1>();
//      for(auto i : Dune::range(re.size(1)))
//      {
//          orientations[0][i] = computeLineOrientation(globalVertexIndicesOfFace(i, 1));
//      }
//      return orientations;
//    }
//    if constexpr (dim==3)
//    {
//      auto orientations = std::array<std::array<std::size_t, 6>, 1>();
//      for(auto i : Dune::range(re.size(1)))
//      {
//        if (re.type(i, 1).isTriangle())
//          orientations[i] = computeTriangleOrientation(globalVertexIndicesOfFace(i, 1));
//      }
//      return orientations;
//    }
//  }

  template<class Node, class It>
  It indices(const Node& node, It it) const
  {
    const auto& element = node.element();
    const auto& localCoefficients = node.finiteElement().localCoefficients();
    const auto& indexSet = gridView().indexSet();
    const auto& re = Dune::referenceElement<double,dim>(element.type());
    auto globalVertexIndicesOfFace = [&](auto faceIndex, auto faceCodim) {
      return Dune::transformedRangeView(re.subEntities(faceIndex, faceCodim, dim), [&](auto localVertexIndex) {
        return indexSet.subIndex(element, localVertexIndex, dim);
      });
    };

    // Instead of computing the orientation of the subentity for each
    // affected DOF, we may want to precompute them once per subentity
    // in advance
//    auto faceOrientations = computeFaceOrientations(node.element());

    for(auto localIndex : Dune::range(localCoefficients.size()))
    {
      auto localKey = localCoefficients.localKey(localIndex);
      auto globalBaseIndex = mapper_.subIndex(element, localKey.subEntity(), localKey.codim());
      auto globalIndex = globalBaseIndex + localKey.index();
      if constexpr(dim == 2)
      {
        if ((localKey.codim() == 1) and (order()>2))
        {
          auto orientation = computeLineOrientation(globalVertexIndicesOfFace(localKey.subEntity(), localKey.codim()));
          globalIndex = globalBaseIndex + flipLineDOF(localKey.index(), orientation);
        }
      }
      if constexpr(dim == 3)
      {
        if ((localKey.codim() == 2) and (order()>2))
        {
          auto orientation = computeLineOrientation(globalVertexIndicesOfFace(localKey.subEntity(), localKey.codim()));
          globalIndex = globalBaseIndex + flipLineDOF(localKey.index(), orientation);
        }
        else if (re.type(localKey.subEntity(), localKey.codim()).isTriangle() and (order()>3))
        {
          auto orientation = computeTriangleOrientation(globalVertexIndicesOfFace(localKey.subEntity(), localKey.codim()));
          globalIndex = globalBaseIndex + flipTriangleDOF(localKey.index(), orientation);
        }
        else if (re.type(localKey.subEntity(), localKey.codim()).isQuadrilateral() and (order()>2))
        {
          auto orientation = computeQuadrilateralOrientation(globalVertexIndicesOfFace(localKey.subEntity(), localKey.codim()));
          globalIndex = globalBaseIndex + flipQuadrilateralDOF(localKey.index(), orientation);
        }
      }
      *it = {{ (size_type)(globalIndex) }};
      ++it;
    }
    return it;
  }

  //! Polynomial order used in the local Lagrange finite-elements
  unsigned int order() const
  {
    return (useDynamicOrder) ? order_ : k;
  }

protected:
  GridView gridView_;

  // Run-time order, only valid if k<0
  unsigned int order_;

  Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper_;
  std::array<std::vector<std::size_t>, 8> triangleFlipTable_;
  std::array<std::vector<std::size_t>, 8> quadrilateralFlipTable_;
};



template<typename GV, int k, typename R>
class LagrangeNode :
  public LeafBasisNode
{
  // Stores LocalFiniteElement implementations with run-time order as a function of GeometryType
  template<typename Domain, typename Range, int dim>
  class LagrangeRunTimeLFECache
  {
  public:
    using FiniteElementType = LagrangeLocalFiniteElement<EquidistantPointSet,dim,Domain,Range>;

    const FiniteElementType& get(GeometryType type)
    {
      auto i = data_.find(type);
      if (i==data_.end())
        i = data_.emplace(type,FiniteElementType(type,order_)).first;
      return (*i).second;
    }

    std::map<GeometryType, FiniteElementType> data_;
    unsigned int order_;
  };

  static constexpr int dim = GV::dimension;
  static constexpr bool useDynamicOrder = (k<0);

  using FiniteElementCache = std::conditional_t<(useDynamicOrder),
                                                       LagrangeRunTimeLFECache<typename GV::ctype, R, dim>,
                                                       LagrangeLocalFiniteElementCache<typename GV::ctype, R, dim, std::max(k,0)>
                                                      >;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  //! Constructor without order (uses the compile-time value)
  LagrangeNode() :
    finiteElement_(nullptr),
    element_(nullptr)
  {}

  //! Constructor with a run-time order
  LagrangeNode(unsigned int order) :
    order_(order),
    finiteElement_(nullptr),
    element_(nullptr)
  {
    // Only the cache for the run-time-order case (i.e., k<0), has the 'order_' member
    if constexpr (useDynamicOrder)
      cache_.order_ = order;
  }

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const
  {
    return *finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(cache_.get(element_->type()));
    this->setSize(finiteElement_->size());
  }

protected:

  unsigned int order() const
  {
    return (useDynamicOrder) ? order_ : k;
  }

  // Run-time order, only valid if k<0
  unsigned int order_;

  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
};



namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a  Lagrange pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam k   The polynomial order of the ansatz functions; -1 means 'order determined at run-time'
 * \tparam R   The range type of the local basis
 */
template<std::size_t k, typename R=double>
auto lagrange()
{
  return [](const auto& gridView) {
    return LagrangePreBasis<std::decay_t<decltype(gridView)>, k, R>(gridView);
  };
}

/**
 * \brief Create a pre-basis factory that can create a  Lagrange pre-basis with a run-time order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R   The range type of the local basis
 */
template<typename R=double>
auto lagrange(int order)
{
  return [=](const auto& gridView) {
    return LagrangePreBasis<std::decay_t<decltype(gridView)>, -1, R>(gridView, order);
  };
}

} // end namespace BasisFactory



/** \brief Nodal basis of a scalar k-th-order Lagrangean finite element space
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * All arguments passed to the constructor will be forwarded to the constructor
 * of LagrangePreBasis.
 *
 * \warning The implementation of the basis with run-time order order uses the
 *   LagrangeFiniteElement implementation of dune-localfunctions, which is known
 *   to violate strict-aliasing rules
 *   (see https://gitlab.dune-project.org/core/dune-localfunctions/issues/14)
 *   Keep this in mind if ever you experience difficult-to-explain crashes
 *   or wrong results.
 *
 * \warning For pyramid elements in 3d, the shape functions are different for
 *    run-time and compile-time order. While the former are defined using the
 *    Duffy-transformation the latter are continuous and piecewise polynomial
 *    with discontinuous gradients.
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis; -1 means 'order determined at run-time'
 * \tparam R The range type of the local basis
 */
template<typename GV, int k=-1, typename R=double>
using LagrangeBasis = DefaultGlobalBasis<LagrangePreBasis<GV, k, R> >;





} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH
