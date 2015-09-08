// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKNODALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKNODALBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/std/final.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/functions/functionspacebases/gridviewfunctionspacebasis.hh>


namespace Dune {
namespace Functions {

template<typename GV, int k>
class PQkNodalBasisLocalView;

template<typename GV, int k, class LocalIndexFunctor>
class PQkNodalBasisLeafNode;

template<typename GV, int k, class MI>
class PQkIndexSet;

template<typename GV, int k, class MI, class LocalIndexFunctor>
class PQkNodeIndexSet;

template<typename GV, int k, class MI>
class PQkNodeFactory
{
  static const int dim = GV::dimension;

public:



  // Precompute the number of dofs per entity type
  const static int dofsPerEdge        = k-1;
  const static int dofsPerTriangle    = (k-1)*(k-2)/2;
  const static int dofsPerQuad        = (k-1)*(k-1);
  const static int dofsPerTetrahedron = ((k-3)*(k-2)*(k-1)/6 > 0) ? (k-3)*(k-2)*(k-1)/6  : 0;
  const static int dofsPerPrism       = (k-1)*(k-1)*(k-2)/2;
  const static int dofsPerHexahedron  = (k-1)*(k-1)*(k-1);
  const static int dofsPerPyramid     = ((k-2)*(k-1)*(2*k-3))/6;



  template<class LocalIndexFunctor>
  using Node = PQkNodalBasisLeafNode<GV, k, LocalIndexFunctor>;

  template<class LocalIndexFunctor>
  using IndexSet = PQkNodeIndexSet<GV, k, MI, LocalIndexFunctor>;

  /** \brief The grid view that the FE space is defined on */
  typedef GV GridView;
  typedef std::size_t size_type;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  /** \brief Constructor for a given grid view object */
  PQkNodeFactory(const GridView& gv) :
    gridView_(gv)
  {
  }


  void initializeIndices()
  {
    vertexOffset_        = 0;
    edgeOffset_          = vertexOffset_          + gridView_.size(dim);
    triangleOffset_      = edgeOffset_            + dofsPerEdge * gridView_.size(dim-1);

    GeometryType triangle;
    triangle.makeTriangle();
    quadrilateralOffset_ = triangleOffset_        + dofsPerTriangle * gridView_.size(triangle);

    Dune::GeometryType quadrilateral;
    quadrilateral.makeQuadrilateral();
    if (dim==3) {
      tetrahedronOffset_   = quadrilateralOffset_ + dofsPerQuad * gridView_.size(quadrilateral);

      GeometryType tetrahedron;
      tetrahedron.makeSimplex(3);
      prismOffset_         = tetrahedronOffset_   +   dofsPerTetrahedron * gridView_.size(tetrahedron);

      GeometryType prism;
      prism.makePrism();
      hexahedronOffset_    = prismOffset_         +   dofsPerPrism * gridView_.size(prism);

      GeometryType hexahedron;
      hexahedron.makeCube(3);
      pyramidOffset_       = hexahedronOffset_    +   dofsPerHexahedron * gridView_.size(hexahedron);
    }
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gridView_;
  }

  template<class LocalIndexFunctor>
  Node<LocalIndexFunctor> node(const LocalIndexFunctor& localIndexFunctor) const
  {
    return Node<LocalIndexFunctor>{localIndexFunctor};
  }

  template<class LocalIndexFunctor>
  IndexSet<LocalIndexFunctor> indexSet() const
  {
    return IndexSet<LocalIndexFunctor>{*this};
  }

  std::size_t size() const
  {
    switch (dim)
    {
      case 1:
        return gridView_.size(1) + dofsPerEdge*gridView_.size(0);
      case 2:
      {
        GeometryType triangle, quad;
        triangle.makeTriangle();
        quad.makeQuadrilateral();
        return gridView_.size(dim) + dofsPerEdge*gridView_.size(1)
             + dofsPerTriangle*gridView_.size(triangle) + dofsPerQuad*gridView_.size(quad);
      }
      case 3:
      {
        GeometryType triangle, quad, tetrahedron, pyramid, prism, hexahedron;
        triangle.makeTriangle();
        quad.makeQuadrilateral();
        tetrahedron.makeTetrahedron();
        pyramid.makePyramid();
        prism.makePrism();
        hexahedron.makeCube(3);
        return gridView_.size(dim) + dofsPerEdge*gridView_.size(2)
             + dofsPerTriangle*gridView_.size(triangle) + dofsPerQuad*gridView_.size(quad)
             + dofsPerTetrahedron*gridView_.size(tetrahedron) + dofsPerPyramid*gridView_.size(pyramid)
             + dofsPerPrism*gridView_.size(prism) + dofsPerHexahedron*gridView_.size(hexahedron);
      }
    }
    DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
  }



//protected:
  const GridView gridView_;

  size_t vertexOffset_;
  size_t edgeOffset_;
  size_t triangleOffset_;
  size_t quadrilateralOffset_;
  size_t tetrahedronOffset_;
  size_t pyramidOffset_;
  size_t prismOffset_;
  size_t hexahedronOffset_;

};


template<typename GV, int k, class MI, class LocalIndexFunctor>
class PQkNodeIndexSet
{
  enum {dim = GV::dimension};

public:
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
//  typedef PQkNodalBasisLocalView<GV,k> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;

  using NodeFactory = PQkNodeFactory<GV, k, MI>;

  using Node = typename NodeFactory::template Node<LocalIndexFunctor>;

  PQkNodeIndexSet(const NodeFactory& nodeFactory) :
    nodeFactory_(&nodeFactory)
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Node& node)
  {
    node_ = &node;
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    node_ = nullptr;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    return node_->finiteElement().size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  const MultiIndex index(size_type i) const
  {
    Dune::LocalKey localKey = node_->finiteElement().localCoefficients().localKey(i);
    const auto& gridIndexSet = nodeFactory_->gridView().indexSet();
    const auto& element = node_->element();

    // The dimension of the entity that the current dof is related to
    size_t dofDim = dim - localKey.codim();

    if (dofDim==0) {  // vertex dof
      return {{ gridIndexSet.subIndex(element,localKey.subEntity(),dim) }};
    }

    if (dofDim==1)
    {  // edge dof
      if (dim==1)   // element dof -- any local numbering is fine
        return {{ nodeFactory_->edgeOffset_ + (k-1)*gridIndexSet.subIndex(element,0,0) + localKey.index() }};
      else
      {
        const Dune::ReferenceElement<double,dim>& refElement
            = Dune::ReferenceElements<double,dim>::general(element.type());

        // we have to reverse the numbering if the local triangle edge is
        // not aligned with the global edge
        size_t v0 = gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),0,dim),dim);
        size_t v1 = gridIndexSet.subIndex(element,refElement.subEntity(localKey.subEntity(),localKey.codim(),1,dim),dim);
        bool flip = (v0 > v1);
        return {{ (flip)
          ? nodeFactory_->edgeOffset_ + (k-1)*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()) + (k-2)-localKey.index()
              : nodeFactory_->edgeOffset_ + (k-1)*gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()) + localKey.index() }} ;
      }
    }

    if (dofDim==2)
    {
      if (dim==2)   // element dof -- any local numbering is fine
      {
        if (element.type().isTriangle())
        {
          const int interiorLagrangeNodesPerTriangle = (k-1)*(k-2)/2;
          return {{ nodeFactory_->triangleOffset_ + interiorLagrangeNodesPerTriangle*gridIndexSet.subIndex(element,0,0) + localKey.index() }};
        }
        else if (element.type().isQuadrilateral())
        {
          const int interiorLagrangeNodesPerQuadrilateral = (k-1)*(k-1);
          return {{ nodeFactory_->quadrilateralOffset_ + interiorLagrangeNodesPerQuadrilateral*gridIndexSet.subIndex(element,0,0) + localKey.index() }};
        }
        else
          DUNE_THROW(Dune::NotImplemented, "2d elements have to be triangles or quadrilaterals");
      } else
      {
        const Dune::ReferenceElement<double,dim>& refElement
            = Dune::ReferenceElements<double,dim>::general(element.type());

        if (k>3)
          DUNE_THROW(Dune::NotImplemented, "PQkNodalBasis for 3D grids is only implemented if k<=3");

        if (k==3 and !refElement.type(localKey.subEntity(), localKey.codim()).isTriangle())
          DUNE_THROW(Dune::NotImplemented, "PQkNodalBasis for 3D grids with k==3 is only implemented if the grid is a simplex grid");

        return {{ nodeFactory_->triangleOffset_ + gridIndexSet.subIndex(element,localKey.subEntity(),localKey.codim()) }};
      }
    }

    if (dofDim==3)
    {
      if (dim==3)   // element dof -- any local numbering is fine
      {
        if (element.type().isTetrahedron())
          return {{ nodeFactory_->tetrahedronOffset_ + NodeFactory::dofsPerTetrahedron*gridIndexSet.subIndex(element,0,0) + localKey.index() }};
        else if (element.type().isHexahedron())
          return {{ nodeFactory_->hexahedronOffset_ + NodeFactory::dofsPerHexahedron*gridIndexSet.subIndex(element,0,0) + localKey.index() }};
        else if (element.type().isPrism())
          return {{ nodeFactory_->prismOffset_ + NodeFactory::dofsPerPrism*gridIndexSet.subIndex(element,0,0) + localKey.index() }};
        else if (element.type().isPyramid())
          return {{ nodeFactory_->pyramidOffset_ + NodeFactory::dofsPerPyramid*gridIndexSet.subIndex(element,0,0) + localKey.index() }};
        else
          DUNE_THROW(Dune::NotImplemented, "3d elements have to be tetrahedra, hexahedra, prisms, or pyramids");
      } else
        DUNE_THROW(Dune::NotImplemented, "Grids of dimension larger than 3 are no supported");
    }
    DUNE_THROW(Dune::NotImplemented, "Grid contains elements not supported for the PQkNodalBasis");
  }

  const NodeFactory* nodeFactory_;

  const Node* node_;
};




template<typename GV, int k, class MI>
class PQkLocalIndexSet
{
  enum {dim = GV::dimension};

public:
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef PQkNodalBasisLocalView<GV,k> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex = MI;


  using NodeIndexSet = PQkNodeIndexSet<GV,k,MultiIndex,ShiftedIdentity<size_type>>;

  template<class NI>
  PQkLocalIndexSet(NI&& nodeIndexSet) :
    nodeIndexSet_(std::forward<NI>(nodeIndexSet))
  {}

//  PQkLocalIndexSet(const NodeIndexSet& nodeIndexSet) :
//    nodeIndexSet_(nodeIndexSet)
//  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const PQkNodalBasisLocalView<GV,k>& localView)
  {
    localView_ = &localView;
    nodeIndexSet_.bind(localView_->tree());
  }

  /** \brief Unbind the view
   */
  void unbind()
  {
    localView_ = nullptr;
    nodeIndexSet_.unbind();
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const
  {
    return nodeIndexSet_.size();
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  const MultiIndex index(size_type i) const
  {
    return nodeIndexSet_.index(i);
  }

  /** \brief Return the local view that we are attached to
   */
  const LocalView& localView() const
  {
    return *localView_;
  }

protected:

  const PQkNodalBasisLocalView<GV,k>* localView_;

  NodeIndexSet nodeIndexSet_;
};

template<typename GV, int k, class MI>
class PQkIndexSet
{
  static const int dim = GV::dimension;

  // Needs the mapper
  friend class PQkLocalIndexSet<GV,k,MI>;

  // Precompute the number of dofs per entity type
  const static int dofsPerEdge        = k-1;
  const static int dofsPerTriangle    = (k-1)*(k-2)/2;
  const static int dofsPerQuad        = (k-1)*(k-1);
  const static int dofsPerTetrahedron = ((k-3)*(k-2)*(k-1)/6 > 0) ? (k-3)*(k-2)*(k-1)/6  : 0;
  const static int dofsPerPrism       = (k-1)*(k-1)*(k-2)/2;
  const static int dofsPerHexahedron  = (k-1)*(k-1)*(k-1);
  const static int dofsPerPyramid     = ((k-2)*(k-1)*(2*k-3))/6;
public:

  typedef PQkLocalIndexSet<GV,k,MI> LocalIndexSet;

  using size_type = std::size_t;
  using LocalIndexFunctor = ShiftedIdentity<size_type>;

  using NodeFactory = PQkNodeFactory<GV,k,MI>;

  PQkIndexSet(const NodeFactory& nodeFactory) :
    nodeFactory_(&nodeFactory),
    gridView_(nodeFactory_->gridView())
  {
    vertexOffset_        = 0;
    edgeOffset_          = vertexOffset_          + gridView_.size(dim);
    triangleOffset_      = edgeOffset_            + dofsPerEdge * gridView_.size(dim-1);

    GeometryType triangle;
    triangle.makeTriangle();
    quadrilateralOffset_ = triangleOffset_        + dofsPerTriangle * gridView_.size(triangle);

    Dune::GeometryType quadrilateral;
    quadrilateral.makeQuadrilateral();
    if (dim==3) {
      tetrahedronOffset_   = quadrilateralOffset_ + dofsPerQuad * gridView_.size(quadrilateral);

      GeometryType tetrahedron;
      tetrahedron.makeSimplex(3);
      prismOffset_         = tetrahedronOffset_   +   dofsPerTetrahedron * gridView_.size(tetrahedron);

      GeometryType prism;
      prism.makePrism();
      hexahedronOffset_    = prismOffset_         +   dofsPerPrism * gridView_.size(prism);

      GeometryType hexahedron;
      hexahedron.makeCube(3);
      pyramidOffset_       = hexahedronOffset_    +   dofsPerHexahedron * gridView_.size(hexahedron);
    }

  }

  std::size_t size() const
  {
    switch (dim)
    {
      case 1:
        return gridView_.size(1) + dofsPerEdge*gridView_.size(0);
      case 2:
      {
        GeometryType triangle, quad;
        triangle.makeTriangle();
        quad.makeQuadrilateral();
        return gridView_.size(dim) + dofsPerEdge*gridView_.size(1)
             + dofsPerTriangle*gridView_.size(triangle) + dofsPerQuad*gridView_.size(quad);
      }
      case 3:
      {
        GeometryType triangle, quad, tetrahedron, pyramid, prism, hexahedron;
        triangle.makeTriangle();
        quad.makeQuadrilateral();
        tetrahedron.makeTetrahedron();
        pyramid.makePyramid();
        prism.makePrism();
        hexahedron.makeCube(3);
        return gridView_.size(dim) + dofsPerEdge*gridView_.size(2)
             + dofsPerTriangle*gridView_.size(triangle) + dofsPerQuad*gridView_.size(quad)
             + dofsPerTetrahedron*gridView_.size(tetrahedron) + dofsPerPyramid*gridView_.size(pyramid)
             + dofsPerPrism*gridView_.size(prism) + dofsPerHexahedron*gridView_.size(hexahedron);
      }

    }

    DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
  }

  LocalIndexSet localIndexSet() const
  {
    return LocalIndexSet(nodeFactory_->template indexSet<LocalIndexFunctor>());
  }

private:

  size_t vertexOffset_;
  size_t edgeOffset_;
  size_t triangleOffset_;
  size_t quadrilateralOffset_;
  size_t tetrahedronOffset_;
  size_t pyramidOffset_;
  size_t prismOffset_;
  size_t hexahedronOffset_;

  const NodeFactory* nodeFactory_;
  const GV gridView_;
};

/** \brief Nodal basis of a scalar third-order Lagrangean finite element space
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - Grids must be 1d, 2d, or 3d
 * - 3d grids must be simplex grids
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k>
class PQkNodalBasis
: public GridViewFunctionSpaceBasis<GV,
                                    PQkNodalBasisLocalView<GV,k>,
                                    PQkIndexSet<GV,k,std::array<std::size_t, 1> >,
                                    std::array<std::size_t, 1> >
{
  static const int dim = GV::dimension;

public:

  /** \brief The grid view that the FE space is defined on */
  typedef GV GridView;
  typedef std::size_t size_type;

  /** \brief Type of the local view on the restriction of the basis to a single element */
  typedef PQkNodalBasisLocalView<GV,k> LocalView;

  /** \brief Type used for global numbering of the basis vectors */
  typedef std::array<size_type, 1> MultiIndex;


  using NodeFactory = PQkNodeFactory<GV, k, MultiIndex >;

  /** \brief Constructor for a given grid view object */
  PQkNodalBasis(const GridView& gv) :
    gridView_(gv),
    nodeFactory_(gridView_),
    indexSet_(nodeFactory_)
  {
    nodeFactory_.initializeIndices();
  }

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const DUNE_FINAL
  {
    return gridView_;
  }

  PQkIndexSet<GV,k, MultiIndex> indexSet() const
  {
    return indexSet_;
  }

  std::size_t size() const
  {
    return nodeFactory_.size();
  }

  /** \brief Return local view for basis
   *
   */
  LocalView localView() const
  {
    return LocalView(this);
  }

  const NodeFactory& nodeFactory() const
  {
    return nodeFactory_;
  }

protected:
  const GridView gridView_;

  NodeFactory nodeFactory_;
  PQkIndexSet<GV,k, MultiIndex> indexSet_;
};


/** \brief The restriction of a finite element basis to a single element */
template<typename GV, int k>
class PQkNodalBasisLocalView
{
public:
  /** \brief The global FE basis that this is a view on */
  typedef PQkNodalBasis<GV,k> GlobalBasis;

  typedef typename GlobalBasis::GridView GridView;

  /** \brief The type used for sizes */
  typedef typename GlobalBasis::size_type size_type;

  /** \brief Type used to number the degrees of freedom
   *
   * In the case of mixed finite elements this really can be a multi-index, but for a standard
   * P3 space this is only a single-digit multi-index, i.e., it is an integer.
   */
  typedef typename GlobalBasis::MultiIndex MultiIndex;

  /** \brief Type of the grid element we are bound to */
  typedef typename GridView::template Codim<0>::Entity Element;

  /** \brief Tree of local finite elements / local shape function sets
   *
   * In the case of a P3 space this tree consists of a single leaf only,
   * i.e., Tree is basically the type of the LocalFiniteElement
   */
  typedef PQkNodalBasisLeafNode<GV,k,ShiftedIdentity<size_type> > Tree;

  /** \brief Construct local view for a given global finite element basis */
  PQkNodalBasisLocalView(const GlobalBasis* globalBasis) :
    globalBasis_(globalBasis),
    tree_(globalBasis->nodeFactory().node(ShiftedIdentity<size_type>(0)))
  {}

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Element& e)
  {
    element_ = &e;
    tree_.bind(e);
  }

  /** \brief Return the grid element that the view is bound to
   *
   * \throws Dune::Exception if the view is not bound to anything
   */
  const Element& element() const
  {
    if (element_)
      return *element_;
    else
      DUNE_THROW(Dune::Exception, "Can't query element of unbound local view");
  }

  /** \brief Unbind from the current element
   *
   * Calling this method should only be a hint that the view can be unbound.
   * And indeed, in the PQkNodalBasisView implementation this method does nothing.
   */
  void unbind()
  {}

  /** \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns Tree // This is tree
   */
  const Tree& tree() const
  {
  return tree_;
  }

  /** \brief Number of degrees of freedom on this element
   */
  size_type size() const
  {
  // We have subTreeSize==lfe.size() because we're in a leaf node.
  return tree_.finiteElement_->size();
  }

  /**
   * \brief Maximum local size for any element on the GridView
   *
   * This is the maximal size needed for local matrices
   * and local vectors, i.e., the result is
   *
   * The method returns k^dim, which is the number of degrees of freedom you get for cubes.
   */
  size_type maxSize() const
  {
  return StaticPower<(k+1),GV::dimension>::power;
  }

  /** \brief Return the global basis that we are a view on
   */
  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }

protected:
  const GlobalBasis* globalBasis_;
  const Element* element_;
  Tree tree_;
};


template<typename GV, int k, class LocalIndexFunctor>
class PQkNodalBasisLeafNode :
  public GridFunctionSpaceBasisLeafNodeInterface<
    typename GV::template Codim<0>::Entity,
    typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, GV::dimension, k>::FiniteElementType,
    typename std::result_of<LocalIndexFunctor(int)>::type>
{
  typedef PQkNodalBasis<GV,k> GlobalBasis;
  static const int dim = GV::dimension;
  static const int maxSize = StaticPower<(k+1),GV::dimension>::power;

  typedef typename GV::template Codim<0>::Entity E;
  typedef typename Dune::PQkLocalFiniteElementCache<typename GV::ctype, double, dim, k> FiniteElementCache;
  typedef typename FiniteElementCache::FiniteElementType FE;
  typedef typename std::result_of<LocalIndexFunctor(int)>::type ST;
  typedef typename GlobalBasis::MultiIndex MI;

  typedef typename GlobalBasis::LocalView LocalView;

  friend LocalView;
  friend class PQkLocalIndexSet<GV,k,MI>;

public:
  typedef GridFunctionSpaceBasisLeafNodeInterface<E,FE,ST> Interface;
  typedef typename Interface::size_type size_type;
  typedef typename Interface::Element Element;
  typedef typename Interface::FiniteElement FiniteElement;

  PQkNodalBasisLeafNode(const LocalIndexFunctor& indexFunctor) :
    finiteElement_(nullptr),
    element_(nullptr),
    indexFunctor_(indexFunctor)
  {}

  //! Return current element, throw if unbound
  const Element& element() const DUNE_FINAL
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const DUNE_FINAL
  {
    return *finiteElement_;
  }

  /** \brief Size of subtree rooted in this node (element-local)
   */
  size_type size() const DUNE_FINAL
  {
    // We have subTreeSize==lfe.size() because we're in a leaf node.
    return finiteElement_->size();
  }

  //! Maps from subtree index set [0..subTreeSize-1] into root index set (element-local) [0..localSize-1]
  size_type localIndex(size_type i) const DUNE_FINAL
  {
    return indexFunctor_(i);
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(cache_.get(element_->type()));
  }

  const LocalIndexFunctor& localIndexFunctor() const
  {
    return indexFunctor_;
  }

  LocalIndexFunctor& localIndexFunctor()
  {
    return indexFunctor_;
  }

protected:

  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
  LocalIndexFunctor indexFunctor_;
};

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PQKNODALBASIS_HH
