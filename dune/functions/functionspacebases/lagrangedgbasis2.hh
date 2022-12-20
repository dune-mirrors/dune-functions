// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dune/common/math.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>




namespace Dune {
namespace Functions {

template<typename GV, int k, typename R>
class LagrangeSimplexNode;

// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   LagrangeSimplexDGPreBasis
//   LagrangeSimplexDGNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename R>
class LagrangeSimplexDGPreBasis
{
  static const int dim = GV::dimension;

  // Requires simplex geometries
  using SGT = Dune::Capabilities::hasSingleGeometryType<typename GV::Grid>;
  static_assert(SGT::v && GeometryType(SGT::topologyId,dim).isSimplex());

public:

  /** \brief The grid view that the FE space is defined on */
  using GridView = GV;
  using size_type = std::size_t;
  using Node = LagrangeSimplexNode<GV, k, R>;

  static constexpr size_type maxMultiIndexSize = 2;
  static constexpr size_type minMultiIndexSize = 2;
  static constexpr size_type multiIndexBufferSize = 2;

  /** \brief Constructor for a given grid view object */
  LagrangeSimplexDGPreBasis(const GridView& gv) :
    gridView_(gv)
  {}


  void initializeIndices() {}

  /** \brief Obtain the grid view that the basis is defined on
   */
  const GridView& gridView() const
  {
    return gridView_;
  }

  void update(const GridView& gv)
  {
    gridView_ = gv;
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{};
  }

  size_type size() const
  {
    return gridView_.size(0);
  }

  //! Return number possible values for next position in multi index
  template<class SizePrefix>
  size_type size(const SizePrefix& prefix) const
  {
    return prefix.size() == 0 ? gridView_.size(0) :
           prefix.size() == 1 ? binomial(k+dim,dim) : 0;
  }

  auto sizeTree() const
  {
    return UniformSizeTree{gridView_.size(), StaticFlatSizeTree<binomial(k+dim,dim)>{}};
  }

  /** \todo This method has been added to the interface without prior discussion. */
  size_type dimension() const
  {
    return size();
  }

  size_type maxNodeSize() const
  {
    return binomial(k+dim,dim);
  }

  template<typename It>
  It indices(const Node& node, It it) const
  {
    auto elementIndex = gridView().indexSet().index(node.element());
    for (size_type i = 0, end = node.finiteElement().size() ; i < end ; ++it, ++i)
      *it = {elementIndex, i};
    return it;
  }

  //! Polynomial order used in the local Lagrange finite-elements
  unsigned int order() const
  {
    return k;
  }

protected:
  GridView gridView_;
};

template<typename GV, int k, typename R>
class LagrangeSimplexNode :
  public LeafBasisNode
{
  static constexpr int dim = GV::dimension;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = Dune::LagrangeSimplexLocalFiniteElement<typename GV::ctype,R,dim,k>;

  //! Constructor without order (uses the compile-time value)
  LagrangeSimplexNode() :
    finiteElement_(),
    element_(nullptr)
  {
    this->setSize(finiteElement_.size());
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
    return finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
  }

protected:
  FiniteElement finiteElement_;
  const Element* element_;
};


namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a LagrangeDG pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam k   The polynomial order of the ansatz functions
 */
template<std::size_t k>
auto lagrangeSimplexDG()
{
  return [](const auto& gridView) {
    return LagrangeSimplexDGPreBasis<std::decay_t<decltype(gridView)>, k>(gridView);
  };
}

} // end namespace BasisFactory



// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/** \brief Basis of a scalar k-th-order Lagrangean-DG finite element space
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis
 */
template<typename GV, int k>
using LagrangeSimplexDGBasis = DefaultGlobalBasis<LagrangeSimplexDGPreBasis<GV, k> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH
