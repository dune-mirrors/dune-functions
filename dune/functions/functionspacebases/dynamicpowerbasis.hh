// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DYNAMICPOWERBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DYNAMICPOWERBASIS_HH

#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/indices.hh>

#include <dune/functions/common/utility.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>



namespace Dune {
namespace Functions {

/**
 * \brief A pre-basis for dynamic power bases
 *
 * This pre-basis represents a power of a given pre-basis.
 * Its node type is a DynamicPowerBasisNodes for the given subnode.
 *
 * \tparam MI  Type to be used for multi-indices
 * \tparam IMS An IndexMergingStrategy used to merge the global indices of the child factories
 * \tparam SPB  The child pre-basis
 */
template<class MI, class IMS, class SPB>
class DynamicPowerPreBasis
  : public PowerPreBasis<MI,IMS,SPB,-1>
{
public:
  //! The child pre-basis
  using SubPreBasis = SPB;

  //! Node type of the children
  using SubNode = typename SubPreBasis::Node;

  //! Template mapping root tree path to type of created tree node
  using Node =  DynamicPowerBasisNode<SubNode>;

public:
  /**
   * \brief Constructor for given child pre-basis objects
   *
   * The child factories will be stored as copies
   */
  template<class... SFArgs,
    enableIfConstructible<SubPreBasis, SFArgs...> = 0>
  explicit DynamicPowerPreBasis (std::size_t children, SFArgs&&... sfArgs)
    : PowerPreBasis(children, std::forward<SFArgs>(sfArgs)...)
  {}

  /**
   * \brief Create tree node
   */
  Node makeNode () const
  {
    auto node = Node{children_};
    for (std::size_t i=0; i<children_; ++i)
      node.setChild(i, subPreBasis_.makeNode());
    return node;
  }
};


namespace BasisFactory {
namespace Imp {

template<class IMS, class ChildFactory>
struct DynamicPowerPreBasisFactory
{
  static const bool isBlocked = std::is_same_v<IMS,BlockedLexicographic>
                             || std::is_same_v<IMS,BlockedInterleaved>;

  static const std::size_t maxChildIndexSize = Factory::requiredMultiIndexSize;
  static const std::size_t requiredMultiIndexSize = isBlocked ? maxChildIndexSize+1 : maxChildIndexSize;

  template<class MultiIndex, class GridView>
  auto makePreBasis (const GridView& gridView) const
  {
    auto childPreBasis = childPreBasisFactory_.template makePreBasis<MultiIndex>(gridView);
    using ChildPreBasis = decltype(childPreBasis);

    return DynamicPowerPreBasis<MultiIndex, IMS, ChildPreBasis>{children_, std::move(childPreBasis)};
  }

  std::size_t children_;
  ChildFactory childPreBasisFactory_;
};

} // end namespace BasisFactory::Imp



/**
 * \brief Create a pre-basis factory that can build a DynamicPowerPreBasis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam ChildPreBasisFactory Types of child pre-basis factory
 * \tparam IndexMergingStrategy An IndexMergingStrategy type
 *
 * \param cpb Child pre-basis factory
 * \param k   Number of children in this power basis
 * \param ims IndexMergingStrategy to be used [default: BasisFactory::blockedInterleaved]
 */
template<class ChildPreBasisFactory, class IndexMergingStrategy = BasisFactory::BlockedInterleaved>
auto dynpower (ChildPreBasisFactory&& cpb, std::size_t k, IndexMergingStrategy ims = {})
{
  using Factory = Imp::DynamicPowerPreBasisFactory<IndexMergingStrategy, ChildPreBasisFactory>;
  return Factory{k, std::forward<ChildPreBasisFactory>(cpb)};
}

} // end namespace BasisFactory

// Backward compatibility
namespace BasisBuilder {

  using namespace BasisFactory;

}


} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DYNAMICPOWERBASIS_HH
