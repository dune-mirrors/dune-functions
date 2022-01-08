// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_TYPETREE2_TREECONTAINER_HH
#define DUNE_TYPETREE2_TREECONTAINER_HH

#include <type_traits>
#include <utility>
#include <functional>

#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/typetree/treepath.hh>
#include <dune/typetree2/transformtree.hh>
#include <dune/typetree2/typetree.hh>

namespace Dune {
namespace TypeTree2 {
namespace Detail {

/**
  * \brief Wrap a tree and transform it into a vector backend
  */
template<class TreeContainer>
class TreeContainerVectorBackend
{
  template<class Tree, class SizeTree>
  static void resizeImpl(Tree& tree, const SizeTree& sizeTree)
  {
    if constexpr(not SizeTree::isLeaf) {
      if constexpr(not SizeTree::isStatic)
        tree.resize(sizeTree.degree());

      Dune::Hybrid::forEach(Dune::range(sizeTree.degree()), [&](auto i) {
        resizeImpl(tree.child(i), sizeTree.child(i));
      });
    }
  }

  template<class T>
  using SizeTreeConcept = decltype((
    std::declval<T>().degree(),
    T::isLeaf,
    T::isStatic,
  true));

public:
  //! Move the passed tree into the internal storage
  explicit TreeContainerVectorBackend(TreeContainer&& treeContainer) :
    treeContainer_{std::move(treeContainer)}
  {}

  //! Move the passed tree into the internal storage
  template<class SizeTree,
    Dune::disableCopyMove<TreeContainerVectorBackend, SizeTree> = 0,
    SizeTreeConcept<SizeTree> = true>
  explicit TreeContainerVectorBackend(const SizeTree& sizeTree) :
    TreeContainerVectorBackend{}
  {
    resize(sizeTree);
  }

  //! Default constructor. The stored container might need to be resized before usage.
  template <class T = TreeContainer,
    std::enable_if_t<std::is_default_constructible_v<T>, bool> = true>
  TreeContainerVectorBackend() :
    treeContainer_{}
  {}

  //! Access the tree childs using a tree-path (const access)
  template<class... T>
  const auto& operator[](const TypeTree::HybridTreePath<T...>& path) const
  {
    return TypeTree::child(treeContainer_, path);
  }

  //! Access the tree childs using a tree-path (mutable access)
  template<class... T>
  auto& operator[](const TypeTree::HybridTreePath<T...>& path)
  {
    return TypeTree::child(treeContainer_, path);
  }

  //! Resize the (nested) container depending on the degree of the tree nodes
  template<class SizeTree, SizeTreeConcept<SizeTree> = true>
  void resize(const SizeTree& sizeTree)
  {
    resizeImpl(treeContainer_, sizeTree);
  }

  //! Get the raw container
  const TreeContainer& data() const
  {
    return treeContainer_;
  }

  //! Get the raw container
  TreeContainer& data()
  {
    return treeContainer_;
  }

private:
  TreeContainer treeContainer_;
};


/**
  * \brief A simple lambda for creating default constructible values from a node
  *
  * This simply returns LeafToValue<Node>{} for a given Node. It's needed
  * because using a lambda expression in a using declaration is not allowed
  * because it's an unevaluated context.
  */
template<template<class Node> class LeafToValue>
struct LeafToDefaultConstructibleValue
{
  template<class Node>
  auto operator()(const Node& node) const
  {
    return LeafToValue<Node>{};
  }
};

} // namespace Detail


/** \addtogroup TypeTree
 *  \{
 */

/**
 * \brief Create container having the same structure as the given tree
 *
 * This class allows to create a nested hybrid container having the same structure
 * as a given type tree. Power nodes are represented as std::array's while composite
 * nodes are represented as Dune::TupleVector's. The stored values for the leaf nodes
 * are creating using a given predicate. For convenience the created container is
 * not returned directly. Instead, the returned object stores the container and
 * provides operator[] access using a HybridTreePath.
 *
 * \param tree The tree which should be mapper to a container
 * \param leafToValue A mapping `Node -> Value` used to generate the stored values
 *                    for the leaves
 *
 * \returns A container matching the tree structure
 */
template <class Tree, class LeafToValue>
auto makeTreeContainer(Tree const& tree, LeafToValue leafToValue)
{
  auto mapNodes = [leafToValue](auto&& node) {
    using Node = std::decay_t<decltype(node)>;
    if constexpr(Node::isLeaf)
      return leafToValue(std::forward<decltype(node)>(node));
    else
      return node;
  };

  using TreeContainer = std::decay_t<decltype(transformTree(tree, mapNodes))>;
  return Detail::TreeContainerVectorBackend<TreeContainer>{transformTree(tree, mapNodes)};
}

/**
 * \brief Create container havin the same structure as the given tree
 *
 * This class allows to create a nested hybrid container having the same structure
 * as a given type tree. Power nodes are represented as std::array's while composite
 * nodes are represented as Dune::TupleVector's. The stored values for the leaf nodes
 * are of the given type Value. For convenience the created container is
 * not returned directly. Instead, the returned object stores the container and
 * provides operator[] access using a HybridTreePath.
 *
 * \tparam Value Type of the values to be stored for the leafs. Should be default constructible.
 * \param leafToValue A predicate used to generate the stored values for the leaves
 *
 * \returns A container matching the tree structure
 */
template<class Value, class Tree>
auto makeTreeContainer(const Tree& tree)
{
  return makeTreeContainer(tree, [](const auto&) {return Value{};});
}

/**
 * \brief Alias to container type generated by makeTreeContainer for given tree type and uniform value type
 */
template<class Value, class Tree>
using UniformTreeContainer = decltype(makeTreeContainer<Value>(std::declval<const Tree&>()));

/**
 * \brief Alias to container type generated by makeTreeContainer for give tree type and when using LeafToValue to create values
 */
template<template<class Node> class LeafToValue, class Tree>
using TreeContainer = decltype(makeTreeContainer(std::declval<const Tree&>(), std::declval<Detail::LeafToDefaultConstructibleValue<LeafToValue>>()));

//! \} group TypeTree

} // namespace TypeTree
} //namespace Dune

#endif // DUNE_TYPETREE2_TREECONTAINER_HH
