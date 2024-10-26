// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_MULTIDOMAIN_RESTRICTED_HH
#define DUNE_FUNCTIONS_MULTIDOMAIN_RESTRICTED_HH

#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/indices.hh>

#include <dune/functions/common/utility.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/dynamicpowerbasis.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/functionspacebases/containerdescriptors.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

#include "domaininfo.hh"


namespace Dune {
namespace Functions {


// *****************************************************************************
// This is the reusable part of the restricted bases. It contains
//
//   RestrictedPreBasis
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

#define RESTRICTED_PARANOIA_CLEAR 0

template<typename Node, typename GridView>
class RestrictedNode :
    public Node
{

public:

  using Element = typename Node::Element;

  RestrictedNode() = default;

  // we should actually also get some kind of
  // subDomainInfo
  RestrictedNode(Node&& subNode, const GridView & gridView) :
    Node(std::forward<Node>(subNode)),
    _gridView(gridView)
  {}

  void bind(const Element& entity, std::size_t& offset)
  {
    // cast to actual implementation
    Node& node = *this;

    // TODO
    // where do we get access to _subDomainInfo?
    if (_gridView.contains(entity))
    {
      // forward to sub node and do the full bind there
      Impl::callNodeBind(node, entity, offset);
    }
    else
    {
#if RESTRICTED_PARANOIA_CLEAR
      // actually we should not need to reset the size, as this was already done in the first visitor of bind
      node.setOffset(offset);
      node.setSize(0);
#endif
    }
  }

private:
  const GridView & _gridView;
};


/**
 * \brief A pre-basis for power bases
 *
 * This pre-basis represents a power of a given pre-basis.
 * Its node type is a PowerBasisNodes for the given subnode.
 *
 * \tparam SPB  The child pre-basis
 */
template<class SubPreBasis>
class RestrictedPreBasis :
    public SubPreBasis
{
  // static const bool isBlocked = std::is_same_v<IMS,BasisFactory::BlockedLexicographic> or std::is_same_v<IMS,BasisFactory::BlockedInterleaved>;

public:

  //! The grid view that the FE basis is defined on
  using GridView = typename SubPreBasis::GridView;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Strategy used to merge the global indices of the child factories
  // using IndexMergingStrategy = IMS;

  //! Template mapping root tree path to type of created tree node
  using Node = RestrictedNode<typename SubPreBasis::Node, GridView>;

  static constexpr size_type maxMultiIndexSize = SubPreBasis::maxMultiIndexSize;
  static constexpr size_type minMultiIndexSize = SubPreBasis::minMultiIndexSize;
  static constexpr size_type multiIndexBufferSize = SubPreBasis::multiIndexBufferSize;

  /**
   * \brief Constructor for given child pre-basis objects
   *
   * The child factories will be stored as copies
   */
  // template<class ,
  //   disableCopyMove<RestrictedPreBasis, SFArgs...> = 0,
  //   enableIfConstructible<SubPreBasis, SFArgs...> = 0>
  // explicit RestrictedPreBasis(SubPreBasis&& subPreBasis) :
  //   SubPreBasis(std::forward<SubPreBasis>(subPreBasis))
  explicit RestrictedPreBasis(const SubPreBasis& subPreBasis) :
    SubPreBasis(subPreBasis)
  {
    static_assert(models<Concept::PreBasis<GridView>, SubPreBasis>(), "Subprebasis passed to RestrictedPreBasis does not model the PreBasis concept.");
  }

  // RestrictedPreBasis(RestrictedPreBasis&&) = default;
  // RestrictedPreBasis(const RestrictedPreBasis&) = default;

  //! Initialize the global indices
  void initializeIndices()
  {
    // is there something we need to update first?
    SubPreBasis::initializeIndices();
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update(const GridView& gv)
  {
#warning (LATER) we get a multidomain gridview and then have to re-extract our subdomain gridview
    assert(false);
    // this gridview comes from the CompositeMultiDomainPreBasis
    // or the PowerMultiDomainPreBasis, this we don't need to update
    // any indexSets, but the pre basis must be made aware of
    // a possible renumbering
    SubPreBasis::update(gv);
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    auto subNode = SubPreBasis::makeNode();
    return Node{std::move(subNode), this->gridView()};
  }

  template<typename It>
  It indices(const Node& node, It it) const
  {
    if (node.size() == 0)
      return it;
    return SubPreBasis::indices(node,it);
  }

  // TODO

  // Do the inherited size methods work?

  // //! Same as size(prefix) with empty prefix
  // size_type size() const
  // {
  //   return size(Dune::ReservedVector<size_type, multiIndexBufferSize>{});
  // }

  // //! Return number of possible values for next position in multi index
  // template<class SizePrefix>
  // size_type size(const SizePrefix& prefix) const
  // {
  //   return sizeImpl(prefix, children_, IndexMergingStrategy{});
  // }

};

/**
 * \brief A pre-basis for power bases
 *
 * This pre-basis represents a power of a given pre-basis.
 * Its node type is a PowerBasisNodes for the given subnode.
 *
 * \tparam SPB  The child pre-basis
 */
template<typename CompositePowerPreBasis, typename MultiDomainGridView>
class MultiDomainPreBasis :
    public CompositePowerPreBasis
{
  std::shared_ptr<MultiDomainGridView> _multiDomainGridView;
public:
  using HostGridView = typename MultiDomainGridView::HostGridView;

  //! The grid view that the FE basis is defined on
  using GridView = MultiDomainGridView;

  explicit MultiDomainPreBasis(
    // CompositePowerPreBasis&& subPreBasis,
    const CompositePowerPreBasis& subPreBasis,
    std::shared_ptr<MultiDomainGridView>& mdgv
    ) :
    //CompositePowerPreBasis(std::forward<CompositePowerPreBasis>(subPreBasis)),
    CompositePowerPreBasis(subPreBasis),
    _multiDomainGridView(mdgv)
  {
  }

  // MultiDomainPreBasis(MultiDomainPreBasis&&) = default;
  // MultiDomainPreBasis(const MultiDomainPreBasis&) = default;

  //! Initialize the global indices
  void initializeIndices()
  {
    // is there something we need to update first?
    CompositePowerPreBasis::initializeIndices();
  }

  //! Obtain the grid view that the basis is defined on
  const MultiDomainGridView& gridView() const
  {
    return * _multiDomainGridView;
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update(const HostGridView& gv)
  {
    _multiDomainGridView->update(gv);
    // this gridview comes from the CompositeMultiDomainPreBasis
    // or the PowerMultiDomainPreBasis, this we don't need to update
    // any indexSets, but the pre basis must be made aware of
    // a possible renumbering
    CompositePowerPreBasis::update(*_multiDomainGridView);
  }

};

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_MULTIDOMAIN_POWERBASIS_HH
