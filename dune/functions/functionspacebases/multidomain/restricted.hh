// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_MULTIDOMAIN_RESTRICTED_HH
#define DUNE_FUNCTIONS_MULTIDOMAIN_RESTRICTED_HH

#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/indices.hh>
#include <dune/common/tupleutility.hh>

#include <dune/localfunctions/restricted.hh>

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

// DEBUG Helper
template<typename Node>
std::string nodeName()
{
  auto substitute = [](std::string& str, const std::string& from, std::string to)
  {
    size_t start_pos = str.find(from);
    while(start_pos != std::string::npos)
    {
      str.replace(start_pos, from.length(), to);
      start_pos = str.find(from);
    }
  };
  std::map<std::string,std::string> translation;
  translation["Dune::Functions::Impl::RestrictedNodeBase<"] = "Restrict< ";
  translation["Dune::Functions::Impl::RestrictedLeafBase<"] = "RestrictLeaf< ";
  translation["Dune::Functions::RestrictedNode<"] = "Restrict< ";
  translation["Dune::Functions::PowerBasisNode<"] = "Power< ";
  translation["Dune::Functions::CompositeBasisNode<"] = "Composite< ";
  translation["Dune::Functions::MultiDomain::RestrictedGridView<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> > >, "] = "";
  translation["Dune::Functions::LagrangeNode<1, double>"] = "Lagrange<P1>";
  translation["Dune::Functions::LagrangeNode<2, double>"] = "Lagrange<P2>";
  translation[", Dune::Functions::MultiDomain::RestrictedGridView<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::YaspGrid<2, Dune::EquidistantCoordinates<double, 2> > const> > >"] = "";
  translation[", 2ul>"] = " >";
  std::string name = className<Node>();

  for (const auto & [from,to] : translation)
  {
    substitute(name,from,to);
  }
  return name;
}

template<typename Node>
std::string nodeName(const Node& node)
{
  return nodeName<Node>();
}

// forward declaration
template<typename Node, typename GridView>
class RestrictedNode;

namespace Impl {

  template<typename Node, typename GridView>
  class RestrictedNodeBase :
    public Node
  {
  public:

    static constexpr bool isRestricted = true;

    using Element = typename Node::Element;

    RestrictedNodeBase() = default;

    // we should actually also get some kind of
    // subDomainInfo
    RestrictedNodeBase(Node&& subNode, const GridView & gridView) :
      Node(std::forward<Node>(subNode)),
      _gridView(gridView)
    {}

    void bind(const Element& entity, std::size_t& offset)
    {
      // cast to actual implementation
      Node& node = *this;

      // std::cout << "RESTRICTED contains? "
      //           << _gridView.contains(entity) << "\t"
      //           << nodeName(*this) << std::endl;

      // forward to sub node and do the full
      // or restricted bind there
      Impl::callNodeBind(node, entity, offset);
    }

  protected:
    const GridView & _gridView;
  };

  template<typename Node, typename GridView>
  class RestrictedLeafNode :
    public RestrictedNodeBase<Node, GridView>
  {
  public:

    using Base = RestrictedNodeBase<Node, GridView>;
    using FiniteElement = RestrictedLocalFiniteElementView<typename Node::FiniteElement>;

    // subDomainInfo
    RestrictedLeafNode(Node&& subNode, const GridView & gridView) :
      Base(std::forward<Node>(subNode), gridView),
      finiteElement_(subNode.finiteElement())
    {}

    bool active() const {
      return this->size() > 0;
    };

    using Element = typename Base::Element;

    void bind(const Element& entity, std::size_t& offset)
    {
      std::cout << "BIND LEAF " << nodeName(*this) << "\n";
      if (this->_gridView.contains(entity))
      {
        std::cout << "ACTIVE " << nodeName(*this) << "\n";
        Base::bind(entity, offset);
        finiteElement_ = Base::finiteElement();
        finiteElement_.activate();
      }
      else
      {
        std::cout << "DEACTIVATE " << nodeName(*this) << "\n";
        finiteElement_.setDefaultGeometry(entity.type());
        finiteElement_.deactivate();
        // set sizes to 0 in this sub tree
        clearSize(*this, offset);
      }
    }

    /** \brief Return the LocalFiniteElement for the element we are bound to
     *
     * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
     *
     * If the current element is not
     */
    const FiniteElement& finiteElement() const
    {
      std::cout << "RESTRICTED finiteElement() ... " << nodeName(*this) << std::endl;
      finiteElement_ = ((const Node&)(*this)).finiteElement();
      std::cout << "NAME " << className<Node>() << std::endl;
      std::cout << "ADDRESS " << &(((const Node&)(*this)).finiteElement()) << std::endl;
      return finiteElement_;
    }

  private:
    mutable FiniteElement finiteElement_;
  };

  template<typename Node, typename GridView>
  using RestrictedNodeBaseClass =
    std::conditional_t<
      std::is_same_v<typename Node::NodeTag, TypeTree::LeafNodeTag>,
      Impl::RestrictedLeafNode<Node,GridView>,
      Impl::RestrictedNodeBase<Node,GridView>
    >;

  template<typename Tree, typename GridView>
  auto restrictTree(Tree && node, const GridView & gridView);

  struct CreateOnStack {};
  struct CreateSharedPtr {};

  template<typename T, typename... Args>
  auto createObject(CreateOnStack, Args&&... args)
  {
    return T{std::forward<Args>(args)...};
  }

  template<typename T, typename... Args>
  auto createObject(CreateSharedPtr, Args&&... args)
  {
    return std::make_shared<T>(std::forward<Args>(args)...);
  }

  template<typename Tree, typename GridView>
  auto transformSubTree(Tree && node, const GridView & gridView, TypeTree::LeafNodeTag)
  {
    return std::forward<Tree>(node);
  }

  template<typename GridView, typename T, std::size_t n>
  auto transformSubTree(PowerBasisNode<T,n> && node, const GridView & gridView, TypeTree::PowerNodeTag)
  {
    // ALTERNATIVE: node.setChild(restrictTree(...), index_constant<k> = {})

    using Node = PowerBasisNode<T,n>;
    // unpack the power node, restrict all children and
    // pack everything again
    const auto & children = node.nodeStorage();
    using ChildType = decltype(restrictTree(std::declval<T>(), gridView, CreateOnStack()));
    using TransformedPowerNode = PowerBasisNode<ChildType,Node::degree()>;
    typename TransformedPowerNode::NodeStorage restrictedChildren;
    for (std::size_t i=0; i<children.size(); i++)
      restrictedChildren[i] = restrictTree(*children[i], gridView, CreateSharedPtr());
    return TransformedPowerNode{restrictedChildren};
  }

  template<typename... Children>
  auto createCompositeBasisNode(std::tuple<std::shared_ptr<Children>...> && storage)
  {
    return CompositeBasisNode<Children...>(
      std::forward<std::tuple<std::shared_ptr<Children>...>>(storage));
  }

  template<typename GridView, typename... T>
  auto transformSubTree(CompositeBasisNode<T...> && node, const GridView & gridView, TypeTree::CompositeNodeTag)
  {
    const auto & children = node.nodeStorage();
    // unpack the composite node, restrict all children and
    // pack everything again
    auto restrictedChildren = transformTuple(
      [&](auto && childNode){
        return restrictTree(*childNode, gridView, CreateSharedPtr());
      }, children);
    return createCompositeBasisNode(std::move(restrictedChildren));
  }

#warning TODO avoid wrapping nodes twice!
  template<typename Tree, typename GridView, typename CreatePolicy>
  auto restrictTree(Tree && node, const GridView & gridView, CreatePolicy policy)
  {
    using Node = std::decay_t<Tree>;
    typename Node::NodeTag tag;
    auto transformedNode = transformSubTree(std::forward<Node>(node), gridView, tag);
    using TransformedNode = decltype(transformedNode);
    using Restricted = RestrictedNode<TransformedNode, GridView>;
    std::cout << "Create " << nodeName<Restricted>() << std::endl;
    return createObject<Restricted>(CreatePolicy(), std::move(transformedNode), gridView);
  }

} // end namespace Impl

template<typename Node, typename GridView>
class RestrictedNode :
    public Impl::RestrictedNodeBaseClass<Node,GridView>
{
  using Base = Impl::RestrictedNodeBaseClass<Node,GridView>;
public:
  using Base::Base;
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
  auto makeNode() const
  {
    auto subNode = SubPreBasis::makeNode();
    return Impl::restrictTree(subNode, this->gridView(), Impl::CreateOnStack());
    // // transform node!
    // // ...
    // return Node{std::move(subNode), this->gridView()};
  }

  using SubNode = typename SubPreBasis::Node;
  //! Template mapping root tree path to type of created tree node
  using Node = decltype(Impl::restrictTree(std::declval<SubNode>(), std::declval<const GridView&>(), std::declval<Impl::CreateOnStack>()));
  // RestrictedNode<typename SubPreBasis::Node, GridView>;

  template<typename It>
  It indices(const Node& node, It it) const
  {
    if (node.size() == 0)
      return it;
    return SubPreBasis::indices((const SubNode&)node,it);
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
private:

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
