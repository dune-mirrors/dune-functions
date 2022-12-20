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
namespace BasisFactory {

struct ElementBlocked : public IndexMergingStrategy {};


/**
 * \brief Create a blocking based on element DOFs.
 *
 * \ingroup FunctionSpaceBasesUtilities
 */
constexpr ElementBlocked elementBlocked()
{
  return {};
}

} // end namespace BasisFactory


template<class IMS, class SPB>
class DGPreBasis
{
  static const int dim = GV::dimension;

  // Requires single geometry type
  using SGT = Dune::Capabilities::hasSingleGeometryType<typename SPB::GridView::Grid>;
  static_assert(SGT::v);

  static constexpr bool isBlocked = std::is_same_v<IMS,BasisFactory::ElementBlocked>;

public:
  //! Strategy used to merge the global indices of the child factories
  using IndexMergingStrategy = IMS;

  //! The child pre-basis
  using SubPreBasis = SPB;

  //! The grid view that the FE basis is defined on
  using GridView = typename SubPreBasis::GridView;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  using Node = PowerBasisNode<typename SubPreBasis::Node,1>;

  static constexpr size_type maxMultiIndexSize = SubPreBasis::maxMultiIndexSize + isBlocked;
  static constexpr size_type minMultiIndexSize = SubPreBasis::minMultiIndexSize + isBlocked;
  static constexpr size_type multiIndexBufferSize = SubPreBasis::multiIndexBufferSize + isBlocked;

  /** \brief Constructor for a given grid view object */
  DGPreBasis(const SPB& spb) :
    subPreBasis_(spb)
  {}

  //! Initialize the global indices
  void initializeIndices()
  {
    subPreBasis_.initializeIndices();
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return subPreBasis_.gridView();
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update(const GridView& gv)
  {
    subPreBasis_.update(gv);
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{subPreBasis_.makeNode()};
  }

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    return size(Dune::ReservedVector<size_type, multiIndexBufferSize>{});
  }

  //! Return number of possible values for next position in multi index

  template<class SizePrefix>
  size_type size(const SizePrefix& prefix) const
  {
    return size(prefix, IndexMergingStrategy{});
  }

  auto sizeTree() const
  {
    return sizeTree(IndexMergingStrategy{});
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
    return indices(node, it, IndexMergingStrategy{});
  }

protected:

  template<class SizePrefix>
  size_type size(const SizePrefix& prefix, BasisFactory::ElementBlocked) const
  {
    if (prefix.size() == 0)
      return gridView().size();

    return subPreBasis_.maxNodeSize();
  }

  template<class SizePrefix, class FlatIMS>
  size_type size(const SizePrefix& prefix, FlatIMS) const
  {
    if (prefix.size() == 0)
      return gridView().size() * subPreBasis_.maxNodeSize();
    else
      return 0;
  }

  auto size(BasisFactory::ElementBlocked) const
  {
    if constexpr(hasStaticNodeSize)
      return UniformSizeTree{gridView().size(), StaticFlatSizeTree<SubPreBasis::maxNodeSize()>{}};
    else
      return UniformSizeTree{gridView().size(), FlatSizeTree{subPreBasis_.maxNodeSize()}};
  }

  template<class FlatIMS>
  auto sizeTree(const SizePrefix& prefix, FlatIMS) const
  {
    return FlatSizeTree{gridView().size() * subPreBasis_.maxNodeSize()};
  }

  template<class It>
  It indices(const Node& node, It it, BasisFactory::ElementBlocked) const
  {
    auto elementIndex = gridView().indexSet().index(node.element());
    for (size_type i = 0, end = node.size() ; i < end ; ++it, ++i)
      *it = {elementIndex, i};
    return it;
  }

  template<class It, class FlatIMS>
  It indices(const Node& node, It it, FlatIMS) const
  {
    auto elementIndex = gridView().indexSet().index(node.element());
    auto blockSize = subPreBasis_.maxNodeSize();
    for (size_type i = 0, end = node.size() ; i < end ; ++it, ++i)
      *it = {elementIndex * blockSize + i};
    return it;
  }

protected:
  SubPreBasis subPreBasis_;
};

namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can build a DGPreBasis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam ChildPreBasisFactory Types of child pre-basis factory
 * \tparam IndexMergingStrategy An IndexMergingStrategy type
 * \param childPreBasisFactory Child pre-basis factory
 * \param ims IndexMergingStrategy to be used
 *
 * This overload can be used to explicitly supply an IndexMergingStrategy.
 */
template<class ChildPreBasisFactory, class IndexMergingStrategy>
auto dg(ChildPreBasisFactory&& childPreBasisFactory, const IndexMergingStrategy&)
{
  return [childPreBasisFactory](const auto& gridView) {
    auto childPreBasis = childPreBasisFactory(gridView);
    return DGPreBasis<IndexMergingStrategy, decltype(childPreBasis)>(std::move(childPreBasis));
  };
}

/**
 * \brief Create a factory builder that can build a DGPreBasis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam ChildPreBasisFactory Types of child pre-basis factory
 * \param childPreBasisFactory Child pre-basis factory
 *
 * This overload will select the BasisFactory::ElementBlocked strategy.
 */
template<class ChildPreBasisFactory>
auto dg(ChildPreBasisFactory&& childPreBasisFactory)
{
  return [childPreBasisFactory](const auto& gridView) {
    auto childPreBasis = childPreBasisFactory(gridView);
    return PowerPreBasis<ElementBlocked, decltype(childPreBasis)>(std::move(childPreBasis));
  };
}

} // end namespace BasisFactory
} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEDGBASIS_HH
