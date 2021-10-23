// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH

#include <type_traits>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>


namespace Dune {
namespace Functions {

// *****************************************************************************
// This is the reusable part of the RefinedBasis. It contains
//
//   RefinedPreBasis
//   RefinedNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

template<class GV, class R=double>
class RefinedNode;

template<class GV, class MI, class R=double>
class RefinedPreBasis;



/**
 * \brief A pre-basis for a PQ-lagrange bases with given order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam MI  Type to be used for multi-indices
 * \tparam R   Range type used for shape function values
 */
template<class GV, class MI, class R>
class RefinedPreBasis
{
  static const int dim = GV::dimension;

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  using Node = RefinedNode<GV, R>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  //! Constructor for a given grid view object and run-time order
  RefinedPreBasis (const GV& gridView, Dune::RefinementIntervals intervals)
    : gridView_{gridView}
    , intervals_{intervals}
  {
    using Element = typename GV::template Codim<0>::Entity;
    using Geometry = typename Element::Geometry;
    for (GeometryTypes t : gridView.indexSet().types())
      refinements_[t] = buildRefinement<GV::dimension, typename Geometry::GlobalCoordinate>(t, t);

    for (int i=0; i<=dim; i++)
    {
      dofsPerCube_[i] = computeDofsPerCube(i);
      dofsPerSimplex_[i] = computeDofsPerSimplex(i);
    }
    dofsPerPrism_ = computeDofsPerPrism();
    dofsPerPyramid_ = computeDofsPerPyramid();
  }

  //! Initialize the global indices
  void initializeIndices ()
  {
    switch (dim)
    {
      case 1:
      {
        break;
      }
      case 2:
      {
        quadrilateralOffset_ = dofsPerSimplex_[2]
          * refinements_[Dune::GeometryTypes::triangle].nElements(intervals_)
          * gridView_.size(Dune::GeometryTypes::triangle);
        break;
      }
      case 3:
      {
        prismOffset_ = dofsPerSimplex_[3]
          * refinements_[Dune::GeometryTypes::tetrahedron].nElements(intervals_)
          * gridView_.size(Dune::GeometryTypes::tetrahedron);

        hexahedronOffset_ = prismOffset_ + dofsPerPrism_
          * refinements_[Dune::GeometryTypes::prism].nElements(intervals_)
          * gridView_.size(Dune::GeometryTypes::prism);

        pyramidOffset_ = hexahedronOffset_ + dofsPerCube_[3]
          * refinements_[Dune::GeometryTypes::hexahedron].nElements(intervals_)
          * gridView_.size(Dune::GeometryTypes::hexahedron);
        break;
      }
    }
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView () const
  {
    return gridView_;
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update (const GridView& gv)
  {
    gridView_ = gv;
  }

  /**
   * \brief Create tree node
   */
  Node makeNode () const
  {
    return Node{order_};
  }

  //! Same as size(prefix) with empty prefix
  size_type size () const
  {
    switch (dim)
    {
      case 1:
      {
        return dofsPerCube_[1]
          * refinements_[Dune::GeometryTypes::line].nElements(intervals_)
          * gridView_.size(0);
      }
      case 2:
      {
        return dofsPerSimplex_[2]
            * refinements_[Dune::GeometryTypes::triangle].nElements(intervals_)
            * gridView_.size(Dune::GeometryTypes::triangle)
          + dofsPerCube_[2]
            * refinements_[Dune::GeometryTypes::quadrilateral].nElements(intervals_)
            * gridView_.size(Dune::GeometryTypes::quadrilateral);
      }
      case 3:
      {
        return dofsPerSimplex_[3]
            * refinements_[Dune::GeometryTypes::tetrahedron].nElements(intervals_)
            * gridView_.size(Dune::GeometryTypes::tetrahedron)
          + dofsPerPyramid_
            * refinements_[Dune::GeometryTypes::pyramid].nElements(intervals_)
            * gridView_.size(Dune::GeometryTypes::pyramid)
          + dofsPerPrism_
            * refinements_[Dune::GeometryTypes::prism].nElements(intervals_)
            * gridView_.size(Dune::GeometryTypes::prism)
          + dofsPerCube_[3]
            * refinements_[Dune::GeometryTypes::hexahedron].nElements(intervals_)
            * gridView_.size(Dune::GeometryTypes::hexahedron);
      }
    }
    DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
  }

  //! Return number of possible values for next position in multi index
  size_type size (const SizePrefix prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? size() : 0;
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension () const
  {
    return size();
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize () const
  {
    // That cast to unsigned int is necessary because GV::dimension is an enum,
    // which is not recognized by the power method as an integer type...
    return power(order()+1, (unsigned int)GV::dimension);
  }

  template<typename It>
  It indices (const Node& node, It it) const
  {
    auto const& element = node.element();
    auto const& refinement = refinements_[element.type()];
    auto const& fe = node.finiteElement();

    auto const eIndex = gridIndexSet.index(element);

    // traverse all refined elements
    auto eIt = refinement.eBegin(intervals_), end_eIt != refinement.eEnd(intervals_);
    auto nSubElements = refinement.nElements(intervals_);
    for (; eIt != end_eIt ++eIt)
    {
      auto shift = (eIndex * nSubElements + eIt.index());

      // traverse all DOFs per element
      for (size_type i = 0, end = fe.size(); i < end ; ++it, ++i)
      {
        switch (dim)
        {
          case 1:
          {
            *it = {dofsPerSimplex_[1] * shift + i};
            continue;
          }
          case 2:
          {
            if (element.type().isTriangle())
            {
              *it = {dofsPerSimplex_[2] * shift + i};
              continue;
            }
            else if (element.type().isQuadrilateral())
            {
              *it = {quadrilateralOffset_ + dofsPerCube[2] * shift + i};
              continue;
            }
            break;
          }
          case 3:
          {
            if (element.type().isTetrahedron())
            {
              *it = {dofsPerSimplex_[3] * shift + i};
              continue;
            }
            else if (element.type().isPrism())
            {
              *it = {prismOffset_ + dofsPerPrism_ * shift + i};
              continue;
            }
            else if (element.type().isHexahedron())
            {
              *it = {hexahedronOffset_ + dofsPerCube_[3] * shift + i};
              continue;
            }
            else if (element.type().isPyramid())
            {
              *it = {pyramidOffset_ + dofsPerPyramid_ * shift + i};
              continue;
            }
            break;
          }
        }
        DUNE_THROW(Dune::NotImplemented,
          "No index method for " << element.type() << " available yet!");
      }
    }
    return it;
  }

  //! Polynomial order used in the local Refined finite-elements
  unsigned int order() const
  {
    return (useDynamicOrder) ? order_ : k;
  }

protected:
  GridView gridView_;

  // Run-time order, only valid if k<0
  const unsigned int order_;

  //! Number of degrees of freedom assigned to a simplex (without the ones assigned to its faces!)
  size_type computeDofsPerSimplex(std::size_t simplexDim) const
  {
    return order() == 0 ? (dim == simplexDim ? 1 : 0) : Dune::binomial(std::size_t(order()-1),simplexDim);
  }

  //! Number of degrees of freedom assigned to a cube (without the ones assigned to its faces!)
  size_type computeDofsPerCube(std::size_t cubeDim) const
  {
    return order() == 0 ? (dim == cubeDim ? 1 : 0) : Dune::power(order()-1, cubeDim);
  }

  size_type computeDofsPerPrism() const
  {
    return order() == 0 ? (dim == 3 ? 1 : 0) : (order()-1)*(order()-1)*(order()-2)/2;
  }

  size_type computeDofsPerPyramid() const
  {
    return order() == 0 ? (dim == 3 ? 1 : 0) : (order()-2)*(order()-1)*(2*order()-3)/6;
  }

  // When the order is given at run-time, the following numbers are pre-computed:
  std::array<size_type,dim+1> dofsPerSimplex_;
  std::array<size_type,dim+1> dofsPerCube_;
  size_type dofsPerPrism_;
  size_type dofsPerPyramid_;

  size_type quadrilateralOffset_;
  size_type prismOffset_;
  size_type hexahedronOffset_;
  size_type pyramidOffset_;

};


template<class D, class R>
class RefinedP1LocalBasis
{
public:
  //! \brief export type traits for function signature
  using Traits = LocalBasisTraits<D,1,FieldVector<D,1>,R,1,FieldVector<R,1>,
      FieldMatrix<R,1,1> >;

  //! \brief number of shape functions
  static constexpr unsigned int size ()
  {
    return 3;
  }

  //! \brief Evaluate all shape functions
  inline void evaluateFunction (const typename Traits::DomainType& in,
                                std::vector<typename Traits::RangeType>& out) const
  {
    out.resize(3);

    int subElement;
    typename Traits::DomainType local;
    this->getSubElement(in, subElement, local);

    switch (subElement) {
    case 0 :

      out[0] = 1 - local[0];
      out[1] = local[0];
      out[2] = 0;
      break;

    case 1 :

      out[0] = 0;
      out[1] = 1 - local[0];
      out[2] = local[0];
      break;

    }

  }

  //! \brief Evaluate Jacobian of all shape functions
  inline void
  evaluateJacobian (const typename Traits::DomainType& in,         // position
                    std::vector<typename Traits::JacobianType>& out) const      // return value
  {
    out.resize(3);

    int subElement;
    typename Traits::DomainType local;
    this->getSubElement(in, subElement, local);

    switch (subElement) {
    case 0 :

      out[0][0][0] = -2;
      out[1][0][0] =  2;
      out[2][0][0] =  0;
      break;

    case 1 :

      out[0][0][0] =  0;
      out[1][0][0] = -2;
      out[2][0][0] =  2;
      break;

    }
  }

  //! \brief Evaluate partial derivatives of all shape functions
  void partial (const std::array<unsigned int, 1>& order,
                const typename Traits::DomainType& in,         // position
                std::vector<typename Traits::RangeType>& out) const      // return value
  {
    auto totalOrder = order[0];
    if (totalOrder == 0) {
      evaluateFunction(in, out);
    } else if (totalOrder == 1)
    {
      out.resize(3);

      int subElement;
      typename Traits::DomainType local;
      this->getSubElement(in, subElement, local);

      switch (subElement) {
        case 0:
          out[0] = -2;
          out[1] =  2;
          out[2] =  0;
          break;
        case 1:
          out[0] =  0;
          out[1] = -2;
          out[2] =  2;
          break;
      }
    } else {
      out.resize(3);
      out[0] = out[1] = out[2] = 0;
    }
  }

  /** \brief Polynomial order of the shape functions
      Doesn't really apply: these shape functions are only piecewise linear
    */
  static constexpr unsigned int order ()
  {
    return 1;
  }
};



/** \brief Piecewise linear continuous Lagrange functions on a uniformly refined simplex element
 *
 * \tparam D Number type used for domain coordinates
 * \tparam R Number type used for shape function values
 * \tparam dim Dimension of the domain
 */
template<class D, class R, int dim>
class RefinedP1LocalFiniteElement
{
  using LB = RefinedP1LocalBasis<D,R,dim>;
  using LC = RefinedP1LocalCoefficients<dim>;
  using LI = RefinedP1LocalInterpolation<LB>;

public:
  /** \brief Export all types used by this implementation */
  using Traits = LocalFiniteElementTraits<LB,LC,LI>;

  /** \brief Default constructor */
  RefinedP1LocalFiniteElement (GeometryType type, RefinementIntervals intervals)
    : type_{type}
    , intervals_{intervals}
    , basis_{type, intervals}
    , coefficients_{type, intervals}
    , interpolation_{type, intervals}
  {}

  /** \brief The set of shape functions */
  const LB& localBasis () const
  {
    return basis_;
  }

  /** \brief Produces the assignments of the degrees of freedom to the element subentities */
  const LC& localCoefficients () const
  {
    return coefficients_;
  }

  /** \brief Evaluates all degrees of freedom for a given function */
  const LI& localInterpolation () const
  {
    return interpolation_;
  }

  /** \brief Number of shape functions of this finite element */
  unsigned int size () const
  {
    return basis_.size();
  }

  /** \brief The element type that this finite element is defined on */
  GeometryType type () const
  {
    return type_;
  }

private:
  GeometryType type_;
  RefinementIntervals intervals_;

  LB basis_;
  LC coefficients_;
  LI interpolation_;
};





template<typename GV, int k, typename R>
class RefinedNode :
  public LeafBasisNode
{
  // Stores LocalFiniteElement implementations with run-time order as a function of GeometryType
  template<typename Domain, typename Range, int dim>
  class RefinedRunTimeLFECache
  {
  public:
    using FiniteElementType = RefinedLocalFiniteElement<EquidistantPointSet,dim,Domain,Range>;

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

  using FiniteElementCache = typename std::conditional<(useDynamicOrder),
                                                       RefinedRunTimeLFECache<typename GV::ctype, R, dim>,
                                                       PQkLocalFiniteElementCache<typename GV::ctype, R, dim, k>
                                                      >::type;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  //! Constructor without order (uses the compile-time value)
  RefinedNode() :
    finiteElement_(nullptr),
    element_(nullptr)
  {}

  //! Constructor with a run-time order
  RefinedNode(unsigned int order) :
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

namespace Imp {

template<int k, typename R=double>
class RefinedPreBasisFactory
{
  static const bool useDynamicOrder = (k<0);
public:
  static const std::size_t requiredMultiIndexSize = 1;

  // \brief Constructor for factory with compile-time order
  RefinedPreBasisFactory() : order_(0)
  {}

  // \brief Constructor for factory with run-time order (template argument k is disregarded)
  RefinedPreBasisFactory(unsigned int order)
  : order_(order)
  {}

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    return (useDynamicOrder)
      ? RefinedPreBasis<GridView, k, MultiIndex, R>(gridView, order_)
      : RefinedPreBasis<GridView, k, MultiIndex, R>(gridView);
  }

private:
  unsigned int order_;
};

} // end namespace BasisFactory::Imp



/**
 * \brief Create a pre-basis factory that can create a  Refined pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam k   The polynomial order of the ansatz functions; -1 means 'order determined at run-time'
 * \tparam R   The range type of the local basis
 */
template<std::size_t k, typename R=double>
auto lagrange()
{
  return Imp::RefinedPreBasisFactory<k,R>();
}

/**
 * \brief Create a pre-basis factory that can create a  Refined pre-basis with a run-time order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R   The range type of the local basis
 */
template<typename R=double>
auto lagrange(int order)
{
  return Imp::RefinedPreBasisFactory<-1,R>(order);
}

} // end namespace BasisFactory



/** \brief Nodal basis of a scalar k-th-order Refinedan finite element space
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * All arguments passed to the constructor will be forwarded to the constructor
 * of RefinedPreBasis.
 *
 * \warning The implementation of the basis with run-time order order uses the
 *   RefinedFiniteElement implementation of dune-localfunctions, which is known
 *   to violate strict-aliasing rules
 *   (see https://gitlab.dune-project.org/core/dune-localfunctions/issues/14)
 *   Keep this in mind if ever you experience difficult-to-explain crashes
 *   or wrong results.
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis; -1 means 'order determined at run-time'
 * \tparam R The range type of the local basis
 */
template<typename GV, int k=-1, typename R=double>
using RefinedBasis = DefaultGlobalBasis<RefinedPreBasis<GV, k, FlatMultiIndex<std::size_t>, R> >;





} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH
