// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TAYLORHOODBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TAYLORHOODBASIS_HH

#include <dune/common/exceptions.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/indices.hh>

#include <dune/typetree/powernode.hh>
#include <dune/typetree/compositenode.hh>

#include <dune/functions/functionspacebases/nodes.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

namespace Dune {
namespace Functions {


// *****************************************************************************
// This is the reusable part of the basis. It contains
//
//   TaylorHoodPreBasis
//   TaylorHoodBasisTree
//   TaylorHoodVelocityTree
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

template<typename GV>
class TaylorHoodVelocityTree;

template<typename GV>
class TaylorHoodBasisTree;

/**
 * \brief Pre-basis for lowest order Taylor-Hood basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV The grid view that the FE basis is defined on
 * \tparam MI Type to be used for multi-indices
 * \tparam HI Flag to select hybrid indices
 *
 * \note This mainly serves as an example, since you can construct a pre-basis with
 * the same functionality manually using
 * \code
 * static const int k = 1;
 * using VelocityPreBasis = PowerPreBasis<MI,IMS,LagrangePreBasis<GV,k+1,MI>,dim>;
 * using PressurePreBasis = LagrangePreBasis<GV,k,MI>;
 * using TaylorHoodKPreBasis = CompositePreBasis<MI, BlockedLexicographic, VelocityPreBasis, PressurePreBasis>;
 * \endcode
 * Where IMS is BlockedInterleaved if HI is set and
 * FlatInterleaved otherwise.
 */
template<typename GV, class MI, bool HI=false>
class TaylorHoodPreBasis
{
  static const bool useHybridIndices = HI;

  static const int dim = GV::dimension;

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  using Node = TaylorHoodBasisTree<GV>;

  //! Type of created tree node index set. \deprecated
  using IndexSet = Impl::DefaultNodeIndexSet<TaylorHoodPreBasis>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, 2>;

private:

  using PQMultiIndex = std::array<size_type, 1>;

public:

  //! PreBasis type of the pressure in the taylorhood space
  using PQ1PreBasis = LagrangePreBasis<GV,1,PQMultiIndex>;

  //! PreBasis type of the velocity components in the taylorhood space
  using PQ2PreBasis = LagrangePreBasis<GV,2,PQMultiIndex>;

public:

  //! Constructor for a given grid view object
  TaylorHoodPreBasis(const GridView& gv) :
    gridView_(gv),
    pq1PreBasis_(gv),
    pq2PreBasis_(gv)
  {}

  //! Initialize the global indices
  void initializeIndices()
  {
    pq1PreBasis_.initializeIndices();
    pq2PreBasis_.initializeIndices();
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return gridView_;
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update (const GridView& gv)
  {
    pq1PreBasis_.update(gv);
    pq2PreBasis_.update(gv);
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{};
  }

  /**
   * \brief Create tree node index set
   *
   * Create an index set suitable for the tree node obtained
   * by makeNode().
   * \deprecated
   */
  [[deprecated("Warning: The IndexSet typedef and the makeIndexSet method are deprecated. "\
               "As a replacement use the indices() method of the PreBasis directly.")]]
  IndexSet makeIndexSet() const
  {
    return IndexSet{*this};
  }

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    return 2;
  }

  //! Return number of possible values for next position in multi index
  size_type size(const SizePrefix prefix) const
  {
    return sizeImp<useHybridIndices>(prefix);
  }

private:

  template<bool hi,
    typename std::enable_if<not hi,int>::type = 0>
  size_type sizeImp(const SizePrefix prefix) const
  {
    if (prefix.size() == 0)
      return 2;
    if (prefix.size() == 1)
    {
      if (prefix[0] == 0)
        return dim * pq2PreBasis_.size();
      if (prefix[0] == 1)
        return pq1PreBasis_.size();
    }
    assert(prefix.size() == 2);
    return 0;
  }

  template<bool hi,
    typename std::enable_if<hi,int>::type = 0>
  size_type sizeImp(const SizePrefix prefix) const
  {
    if (prefix.size() == 0)
      return 2;
    if (prefix.size() == 1)
    {
      if (prefix[0] == 0)
        return pq2PreBasis_.size();
      if (prefix[0] == 1)
        return pq1PreBasis_.size();
    }
    if (prefix.size() == 2)
    {
      if (prefix[0] == 0)
        return dim;
      if (prefix[0] == 1)
        return 0;
    }
    assert(prefix.size() == 3);
    return 0;
  }

public:

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return dim * pq2PreBasis_.size() + pq1PreBasis_.size();
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    return dim * pq2PreBasis_.maxNodeSize() + pq1PreBasis_.maxNodeSize();
  }

  template<typename It>
  It indices(const Node& node, It it) const
  {
    return indicesImp<useHybridIndices>(node, it);
  }

protected:

  static const void multiIndexPushFront(MultiIndex& M, size_type M0)
  {
    M.resize(M.size()+1);
    for(std::size_t i=M.size()-1; i>0; --i)
      M[i] = M[i-1];
    M[0] = M0;
  }

  template<bool hi, class It,
    typename std::enable_if<not hi,int>::type = 0>
  It indicesImp(const Node& node, It multiIndices) const
  {
    using namespace Dune::Indices;
    for(std::size_t child=0; child<dim; ++child)
    {
      size_type subTreeSize = node.child(_0, 0).size();
      pq2PreBasis_.indices(node.child(_0, 0), multiIndices);
      for (std::size_t i = 0; i<subTreeSize; ++i)
      {
        multiIndexPushFront(multiIndices[i], 0);
        multiIndices[i][1] = multiIndices[i][1]*dim + child;
      }
      multiIndices += subTreeSize;
    }
    size_type subTreeSize = node.child(_1).size();
    pq1PreBasis_.indices(node.child(_1), multiIndices);
    for (std::size_t i = 0; i<subTreeSize; ++i)
      multiIndexPushFront(multiIndices[i], 1);
    multiIndices += subTreeSize;
    return multiIndices;
  }

  template<bool hi, class It,
    typename std::enable_if<hi,int>::type = 0>
  It indicesImp(const Node& node, It multiIndices) const
  {
    using namespace Dune::Indices;
    for(std::size_t child=0; child<dim; ++child)
    {
      size_type subTreeSize = node.child(_0, 0).size();
      pq2PreBasis_.indices(node.child(_0, 0), multiIndices);
      for (std::size_t i = 0; i<subTreeSize; ++i)
      {
        multiIndexPushFront(multiIndices[i], 0);
        multiIndices[i].push_back(i);
      }
      multiIndices += subTreeSize;
    }
    size_type subTreeSize = node.child(_1).size();
    pq1PreBasis_.indices(node.child(_1), multiIndices);
    for (std::size_t i = 0; i<subTreeSize; ++i)
      multiIndexPushFront(multiIndices[i], 1);
    multiIndices += subTreeSize;
    return multiIndices;
  }

  GridView gridView_;

  PQ1PreBasis pq1PreBasis_;
  PQ2PreBasis pq2PreBasis_;
};



template<typename GV>
class TaylorHoodVelocityTree :
    public PowerBasisNode<LagrangeNode<GV,2>, GV::dimension>
{
  using PQ2Node = LagrangeNode<GV,2>;
  using Base = PowerBasisNode<PQ2Node, GV::dimension>;

public:
  TaylorHoodVelocityTree()
  {
    for(int i=0; i<GV::dimension; ++i)
      this->setChild(i, std::make_shared<PQ2Node>());
  }
};

template<typename GV>
class TaylorHoodBasisTree :
    public CompositeBasisNode<
      TaylorHoodVelocityTree<GV>,
      LagrangeNode<GV,1>
    >
{
  using VelocityNode=TaylorHoodVelocityTree<GV>;
  using PressureNode=LagrangeNode<GV,1>;

  using Base=CompositeBasisNode<VelocityNode, PressureNode>;

public:
  TaylorHoodBasisTree()
  {
    this->template setChild<0>(std::make_shared<VelocityNode>());
    this->template setChild<1>(std::make_shared<PressureNode>());
  }
};



namespace BasisFactory {

namespace Imp {

class TaylorHoodPreBasisFactory
{
public:
  static const std::size_t requiredMultiIndexSize=2;

  template<class MultiIndex, class GridView>
  auto makePreBasis(const GridView& gridView) const
  {
    return TaylorHoodPreBasis<GridView, MultiIndex>(gridView);
  }

};

} // end namespace BasisFactory::Imp

/**
 * \brief Create a pre-basis factory that can create a Taylor-Hood pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 */
auto taylorHood()
{
  return Imp::TaylorHoodPreBasisFactory();
}

} // end namespace BasisFactory

// *****************************************************************************
// This is the actual global basis implementation based on the reusable parts.
// *****************************************************************************

/**
 * \brief Nodal basis for a lowest order Taylor-Hood Lagrangean finite element space
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV The GridView that the space is defined on.
 *
 * \note This mainly serves as an example, since you can construct a basis with
 * the same functionality manually using
 * \code
 * static const int k = 1;
 * auto taylorHoodBasis = makeBasis(
 *   gridView,
 *   composite(
 *     power<dim>(
 *       lagrange<k+1>(),
 *       flatInterleaved()),
 *     lagrange<k>()
 *   ));
 * \endcode
 */
template<typename GV>
using TaylorHoodBasis = DefaultGlobalBasis<TaylorHoodPreBasis<GV, Dune::ReservedVector<std::size_t, 2>> >;



} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TAYLORHOODBASIS_HH
