// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH

#include <array>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/typelist.hh>
#include <dune/common/indices.hh>
#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/std/mdarray.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lagrange/lagrangecube.hh>
#include <dune/localfunctions/lagrange/lagrangeprism.hh>
#include <dune/localfunctions/lagrange/lagrangepyramid.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/functions/common/faceorientation.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/leafprebasismixin.hh>
#include <dune/functions/functionspacebases/transformed/bindcontext.hh>
#include <dune/functions/functionspacebases/transformed/geometryderivative.hh>
#include <dune/functions/functionspacebases/transformed/localfiniteelement.hh>
#include <dune/functions/functionspacebases/transformed/pipeline.hh>

#include <dune/grid/common/capabilities.hh>

namespace Dune {
namespace Functions {



namespace Experimental {


/**
 * \brief A class for permuting Lagrange DOFs of arbitrary order for dimensions 1,2,3
 *
 * \tparam IndexIdSet An IndexSet or IdSet used to define globally unique face orientations
 *
 * This class can be used to permute local Lagrange DOFs to get local per-face
 * DOF indices that are unique among adjacent elements.
 * While permutation of edge DOFs corresponds to a simple flip, permutation of
 * DOFs associated to triangular and quadrilateral faces is more involved and
 * thus looked up from pre-computed tables.
 */
template<class IndexIdSet>
class LagrangeFaceDOFPermutation
{
  static constexpr int dim = IndexIdSet::dimension;

  // The following methods pre-compute the tables for renumbering the face DOFs.
  // The rows of a table correspond to the local orientations which are indexed
  // using the face orientation index as provided by the FaceOrientations class.
  // Each row stores the local indices of the face DOF indices wrt the global orientation.
  // I.e. given the table T, a local face DOF index i with respect to the
  // local orientation j, the globally unique index of this DOF with respect to
  // the global face orientation is T[j][i].
  //
  // The computation of the table is based on the following observation:
  // Assume that F_j is the transformation mapping the global orientation
  // of the face to its local orientation with index j. Now let p_i an
  // interpolation point in the face with local index i. Then the
  // globally oriented index of the basis function associated to p_i
  // is obtained as the local index of the point F_j(p_i). Hence we can
  // compute the renumbering table T[j][i] as follows:
  // For given local orientation j and locally oriented DOF index i
  // we first compute the point p_i, then we apply the transformation
  // F_j (not its inverse as one might guess at first sight) which
  // consists of rotations and zero or one reflection and use the index
  // of the obtained point.

  // Compute table of all possible permuted triangle DOF indices
  static auto globallyOrientedTriangleDOFTable(int order)
  {
    auto dofPermutation = std::array<std::vector<std::size_t>, 8>();
    // Number of DOFs on the base line
    std::size_t m = order-2;
    // Total number of DOFs
    std::size_t n = m*(m+1)/2;
    if (order<3)
      return dofPermutation;
    // In order to reflect we need to apply the transformation
    // from global to local orientation as documented in the
    // FaceOrientations class for each face orientation index.
    // To do this we first map flat DOF indices to barycentric
    // multi-indices wrt. the vertices and summing up to m-1.
    // Then all transformations can be represented by flips
    // of the components for the barycentric multi-index.
    using BarycentricIndex = std::array<std::size_t, 3>;
    auto barycentricIndices = std::vector<BarycentricIndex>();
    for(auto i : Dune::range(m))
      for(auto j : Dune::range(m-i))
        barycentricIndices.push_back(BarycentricIndex{m-1-(i+j), j, i});
    auto flatToBarycentric = [&](auto i) { return barycentricIndices[i]; };
    auto barycentricToFlat = [&](auto b) { return b[1] + m*b[2] - b[2]*(b[2]-1)/2; };
    // Loop of all 8 orientations
    for(auto j : Dune::range(8))
    {
      dofPermutation[j].resize(n);
      for(auto i : Dune::range(n))
      {
        auto b = flatToBarycentric(i);
        // Bit 0 and 1: Number of counter-clockwise vertex-wise rotations
        if ((j & 0b011) == 1)
        {
          // The counter-clockwise rotation by one vertex can be decomposed
          // into two reflections.
          std::swap(b[0], b[1]);
          std::swap(b[0], b[2]);
        }
        if ((j & 0b011) == 2)
        {
          // The counter-clockwise rotation by two vertices is the inverse
          // of the rotation by one vertex and can thus be decomposed into
          // the same two reflections in opposite order.
          std::swap(b[0], b[2]);
          std::swap(b[0], b[1]);
        }
        // Bit 2: Reflect across xy-diagonal
        if (j & 0b100)
          std::swap(b[1], b[2]);
        dofPermutation[j][i] = barycentricToFlat(b);
      }
    }
    return dofPermutation;
  }

  // Compute table of all possible permuted quadrilateral DOF indices
  static auto globallyOrientedQuadrilateralDOFTable(int order)
  {
    auto dofPermutation = std::array<std::vector<std::size_t>, 8>();
    if (order<2)
      return dofPermutation;
    // Number of DOFs on the base line
    std::size_t m = order-1;
    // Total number of DOFs
    std::size_t n = m*m;
    // In order to reflect we map the flat DOF index into a cartesian
    // multi-index wrt. the x- and y-axis.
    using CartesianIndex = std::array<std::size_t, 2>;
    auto flatToCartesian = [&](auto i) { return CartesianIndex{i % m, i / m}; };
    auto cartesianToFlat = [&](auto c) { return c[0] + m*c[1]; };
    for(auto j : Dune::range(8))
    {
      dofPermutation[j].resize(n);
      for(auto i : Dune::range(n))
      {
        auto c = flatToCartesian(i);
        // Bit 0 and 1: Number of counter-clockwise vertex-wise rotations
        if ((j & 0b011) == 1)
        {
          std::swap(c[0], c[1]);
          c[0] = m-1-c[0];
        }
        if ((j & 0b011) == 2)
        {
          c[0] = m-1-c[0];
          c[1] = m-1-c[1];
        }
        if ((j & 0b011) == 3)
        {
          c[0] = m-1-c[0];
          std::swap(c[0], c[1]);
        }
        // Bit 2: Reflect across xy-diagonal
        if (j & 0b100)
          std::swap(c[0], c[1]);
        dofPermutation[j][i] = cartesianToFlat(c);
      }
    }
    return dofPermutation;
  }

  std::size_t order() const
  {
    return order_;
  }

  /*
   * Type used to store local per-face indices.
   *
   * Since the maximal number of face DOFs is n=(order-1)*(order-1),
   * this type must be able to represent the range 0,...,n-1.
   */
  using FaceIndexType = unsigned char;

public:

  /**
   * \brief Maximal supported polynomial order
   *
   * This is limited by the internally used type to represent
   * face-local indices of DOFs.
   */
  static constexpr unsigned int maxOrder3d = 17;

  /**
   * \param indexSet The gridView indexSet containing the elements and vertices
   * \param order The polynomial order of the shape functions
   */
  LagrangeFaceDOFPermutation(const IndexIdSet& indexSet, int order)
    : indexIdSet_(&indexSet)
    , order_(order)
    , facetFlipTable_((dim < 3) ? 0 : std::max(0, (order-1)*(order-1)))
  {
    if constexpr (dim == 3)
    {
      // Fill table with all possible permuted 2d DOF indices.
      // This is indexed with the faceOrientationIndex which first
      // enumerates the orientations of triangles and then the ones
      // of quadrilateral.
      auto triangleFlipTable = globallyOrientedTriangleDOFTable(order);
      auto quadrilateralFlipTable = globallyOrientedQuadrilateralDOFTable(order);
      for(std::size_t k : Dune::range(8))
        for(std::size_t i : Dune::range(triangleFlipTable[k].size()))
          facetFlipTable_(k,i) = triangleFlipTable[k][i];
      for(std::size_t k : Dune::range(8))
        for(std::size_t i : Dune::range(quadrilateralFlipTable[k].size()))
          facetFlipTable_(8+k,i) = quadrilateralFlipTable[k][i];
    }
  }

  /**
   * \brief Compute \ref FaceOrientations for grid elements
   *
   * The returned \ref FaceOrientations object encodes the orientations
   * of all faces of the given element wrt. to the globally unique
   * orientation provided by the index/id set and can be used to
   * compute permuted face DOF indices using \ref permuteFaceDOF
   */
  template <class Element>
  FaceOrientations<dim> computeFaceOrientations(const Element& element) const
  {
    const auto& re = Dune::referenceElement<double,dim>(element.type());
    auto vertexIds = Dune::transformedRangeView(re.subEntities(0,0,dim), [&](auto localVertexIndex) {
      if constexpr (requires { indexIdSet_->subIndex(element, 0, dim); })
        return indexIdSet_->subIndex(element, localVertexIndex, dim);
      else if constexpr (requires { indexIdSet_->subId(element, 0, dim); })
        return indexIdSet_->subId(element, localVertexIndex, dim);
    });
    using namespace Dune::Indices;
    auto k = order();
    if ((dim==1) or (k<3))
      return FaceOrientations<dim>(element.type(), vertexIds);
    else if ((k==3) and element.type().isTetrahedron())
      return FaceOrientations<dim>(element.type(), vertexIds, _2);
    else
      return FaceOrientations<dim>(element.type(), vertexIds, _2, _1);
  }

  /**
   * \brief Compute permuted index of a face DOF
   *
   * \param subEntity The index of the subentity the DOF is associated to
   * \param codim The codim of the subentity the DOF is associated to
   * \param index The local index of the DOF within the subentity
   * \param orientations A \ref FaceOrientations object obtained by computeFaceOrientations
   *
   * \returns The permuted local index of the DOF within the subentity
   *
   * The resulting local DOF index is permuted to match the global orientation
   * induced by the index/id set. I.e. it is guaranteed to be consistent among
   * all adjacent grid elements.
   */
  unsigned int permuteFaceDOF(unsigned int subEntity, unsigned int codim, unsigned int index, const FaceOrientations<dim>& orientations) const
  {
    auto k = order();
    if constexpr(dim == 1)
      return index;
    if constexpr(dim == 2)
    {
      // Flip edge DOFs if reorientation is needed
      if ((k>2) and (codim == 1) and orientations.faceOrientationIndex(subEntity, codim))
        return k-2-index;
      return index;
    }
    if constexpr(dim == 3)
    {
      if (k>2)
      {
        // Flip edge DOFs if reorientation is needed
        if ((codim == 2) and orientations.faceOrientationIndex(subEntity, codim))
          return k-2-index;
        // Permute 2d facet DOFs according to lookup table
        if (codim == 1)
          return facetFlipTable_(orientations.faceOrientationIndex(subEntity, codim), index);
      }
      return index;
    }
  }

  /**
   * \brief Compute permuted index of a face DOF
   *
   * \param localKey A \ref Dune::LocalKey identifying a local DOF
   * \param orientations A \ref FaceOrientations object obtained by computeFaceOrientations
   *
   * \returns The permuted local index of the DOF within the subentity
   *
   * The resulting local DOF index is permuted to match the global orientation
   * induced by the index/id set. I.e. it is guaranteed to be consistent among
   * all adjacent grid elements.
   */
  unsigned int permuteFaceDOF(const Dune::LocalKey& localKey, const FaceOrientations<dim>& orientations) const
  {
    return permuteFaceDOF(localKey.subEntity(), localKey.codim(), localKey.index(), orientations);
  }

private:
  const IndexIdSet* indexIdSet_;
  std::size_t order_;
  Std::mdarray<FaceIndexType, Std::extents<std::size_t, 16, std::dynamic_extent>> facetFlipTable_;
};

} // namespace Experimental



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
 * \note
 * For 3d grids the order is limited to 2 if the grid contains
 * prism and pyramid elements, otherwise to 17.
 * For all grids in 1d and 2d any order is supported.
 *
 * \warning
 * For pyramid elements in 3d, the shape functions are piecewise polynomials
 * with discontinuous gradients along the diagonal through the origin.
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
        return type.dim() == (unsigned int)(dim);
      if (type.isSimplex())
        return Dune::binomial(order-1, type.dim());
      if (type.isCube())
        return Dune::power(order-1, type.dim());
      if (type.isPyramid())
        return 0;
      if (type.isPrism() and (order>1))
        return (order-1)*(order-1)*(order-2)/2;
      return 0;
    };
  }

  using FaceDOFPermutation = Experimental::LagrangeFaceDOFPermutation<typename GV::Grid::GlobalIdSet>;

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
  LagrangePreBasis(const GridView& gv, unsigned int runTimeOrder)
    : gridView_(gv)
    , order_(runTimeOrder)
    , mapper_(gridView_, dofLayout(useDynamicOrder ? runTimeOrder : k))
    , faceDOFPermutation_(gridView_.grid().globalIdSet(), useDynamicOrder ? runTimeOrder : k)
  {
    if (!useDynamicOrder && runTimeOrder!=std::numeric_limits<unsigned int>::max())
      DUNE_THROW(RangeError, "Template argument k has to be -1 when supplying a run-time order!");
    if ((order() > 2) and (gridView_.indexSet().size(Dune::GeometryTypes::pyramid) > 0))
      DUNE_THROW(RangeError, "Polynomial order >2 is not supported for grids containing pyramid elements.");
    if ((order() > 2) and (gridView_.indexSet().size(Dune::GeometryTypes::prism) > 0))
      DUNE_THROW(RangeError, "Polynomial order >2 is not supported for grids containing prism elements.");
    if ((dim == 3) and (order() > FaceDOFPermutation::maxOrder3d))
      DUNE_THROW(RangeError, "Polynomial order >" << FaceDOFPermutation::maxOrder3d << " is not supported in 3d");
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
    if ((order() > 2) and (gv.indexSet().size(Dune::GeometryTypes::pyramid) > 0))
      DUNE_THROW(RangeError, "Polynomial order >2 is not supported for grids containing pyramid elements.");
    if ((order() > 2) and (gv.indexSet().size(Dune::GeometryTypes::prism) > 0))
      DUNE_THROW(RangeError, "Polynomial order >2 is not supported for grids containing prism elements.");
    gridView_ = gv;
    mapper_.update(gridView_);
    faceDOFPermutation_ = FaceDOFPermutation(gridView_.grid().globalIdSet(), order());
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{order()};
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

  template<class Node, class It>
  It indices(const Node& node, It it) const
  {
    const auto& element = node.element();
    const auto& localCoefficients = node.finiteElement().localCoefficients();

    // Short-cut for order 1 since directly using the IndexSet is cheaper than the Mapper.
    if (order() == 1)
    {
      for(auto localIndex : Dune::range(localCoefficients.size()))
      {
        auto localKey = localCoefficients.localKey(localIndex);
        auto globalIndex = mapper_.gridView().indexSet().subIndex(element, localKey.subEntity(), localKey.codim());
        *it = {{ (size_type)(globalIndex) }};
        ++it;
      }
      return it;
    }

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
  FaceDOFPermutation faceDOFPermutation_;
};



template<typename GV, int k, typename R>
class LagrangeNode :
  public LeafBasisNode
{
  static constexpr int dim = GV::dimension;

  // A simple cache handing storing exactly one LFE.  This can be
  // used for grids with only a single element type. In contrast to
  // StaticLagrangeLocalFiniteElementCache this also supports
  // Lagrange*LocalFiniteElement with run-time order.
  template <class LFE>
  class SingleLocalFiniteElementCache
  {
    LFE lfe_;
  public:
    using FiniteElementType = LFE;

    template<class... Args>
    SingleLocalFiniteElementCache(Args&&... args)
      : lfe_(std::forward<Args>(args)...)
    {}

    //! Obtain the cached local finite-element.
    const FiniteElementType& get ([[maybe_unused]] Dune::GeometryType type) const
    {
      return lfe_;
    }
  };

  // Utility function to construct the FiniteElementCache type.
  // Since the function is just a helper to generate a type,
  // it hands out a MetaType<T> instead of a raw T.
  static constexpr auto makeCacheType()
  {
    using D = typename GV::ctype;
    if constexpr (Dune::Capabilities::hasSingleGeometryType<typename GV::Grid>::v)
    {
      constexpr auto type = Dune::GeometryType(Dune::Capabilities::hasSingleGeometryType<typename GV::Grid>::topologyId, GV::dimension);
      if constexpr (type.isSimplex())
        return Dune::MetaType<SingleLocalFiniteElementCache<Dune::LagrangeSimplexLocalFiniteElement<D,R,dim,k>>>{};
      else if constexpr (type.isCube())
        return Dune::MetaType<SingleLocalFiniteElementCache<Dune::LagrangeCubeLocalFiniteElement<D,R,dim,k>>>{};
      else if constexpr (type.isPrism())
        return Dune::MetaType<SingleLocalFiniteElementCache<Dune::LagrangePrismLocalFiniteElement<D,R,k>>>{};
      else if constexpr (type.isPyramid())
        return Dune::MetaType<SingleLocalFiniteElementCache<Dune::LagrangePyramidLocalFiniteElement<D,R,k>>>{};
    }
    else
      return Dune::MetaType<Dune::LagrangeLocalFiniteElementCache<D,R,dim,k>>{};
  }

  using FiniteElementCache = typename decltype(makeCacheType())::type;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using ReferenceFiniteElement = typename FiniteElementCache::FiniteElementType;
  using Context = ElementBindContext<Element>;
  using Transformation = BasisEvaluationPipeline<Context,
    GeometryDerivativeStage<typename Element::Geometry>>;
  using FiniteElement = TransformedLocalFiniteElement<ReferenceFiniteElement, Context,
    Transformation, NoInterpolationTransformation, LocalBasisMode::reference>;

  //! Constructor without order (uses the compile-time value)
  LagrangeNode() :
    LagrangeNode(k)
  {}

  //! Constructor with a run-time order
  LagrangeNode(unsigned int order) :
    cache_(order)
  {}

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return context_.element();
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const
  {
    return finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    context_.bind(e);
    finiteElement_.bind(cache_.get(context_.type()), context_);
    this->setSize(finiteElement_.size());
  }

protected:

  FiniteElementCache cache_;
  Context context_;
  FiniteElement finiteElement_;
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
 * All arguments passed to the constructor will be forwarded to the constructor
 * of LagrangePreBasis.
 *
 * \note
 * For 3d grids containing prism and pyramid elements only order
 * 1 and 2 are supported. For all other grids in 1d, 2d, 3d any
 * order is supported.
 *
 * \warning
 * For pyramid elements in 3d, the shape functions are piecewise polynomials
 * with discontinuous gradients along the diagonal through the origin.
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
