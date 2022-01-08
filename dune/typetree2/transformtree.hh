// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_TYPETREE2_TRANSFORMTREE_HH
#define DUNE_TYPETREE2_TRANSFORMTREE_HH

#include <type_traits>
#include <utility>

#include <dune/common/indices.hh>
#include <dune/typetree2/typetree.hh>

namespace Dune {
namespace TypeTree2 {

//! Transformation of the nodes in a type-tree
template<class Tree, class MapNode>
auto transformTree(Tree const& tree, MapNode mapNode)
{
  if constexpr(Tree::isLeaf)
    return mapNode(tree);
  else {
    // not leaf
    auto transformedTree = [&]{
      auto subTree = [&](auto i) { return transformTree(tree[i], mapNode); };
      if constexpr(Tree::isTypeUniform) {
        using SubTree = decltype(subTree(Dune::index_constant<0>{}));
        if constexpr(Tree::isUniform) {
          if constexpr(Tree::isStatic)
            return StaticUniformTypeTree{Tree::degree(), subTree(0)};
          else
            return DynamicUniformTypeTree{Tree::degree(), subTree(0)};
        }
        else {
          // not uniform
          if constexpr(Tree::isStatic) {
            return Dune::unpackIntegerSequence([subTree](auto... ii) {
              return StaticNonUniformTypeTree<SubTree,Tree::degree()>{subTree(ii)...}; },
              std::make_index_sequence<Tree::degree()>{});
          }
          else {
            DynamicNonUniformTypeTree<SubTree> container;
            container.reserve(tree.degree());
            for (std::size_t i = 0; i < tree.degree(); ++i)
              container.emplace_back(subTree(i));
            return container;
          }
        }
      }
      else {
        // not type-uniform
        return Dune::unpackIntegerSequence(
          [subTree](auto... ii) {
            return VariadicNonUniformTypeTree<decltype(subTree(ii))...>{subTree(ii)...};
          }, std::make_index_sequence<Tree::degree()>());
      }
    };

    return mapNode(transformedTree());
  }
}

//! Class that adds a data member and corresponding access methods to a Node base class
template<class Node, class Data>
struct DataNodeMixin
    : public Node
{
  template<class N, class D>
  DataNodeMixin(N&& node, D&& data)
    : Node{std::forward<N>(node)}
    , data_{std::forward<D>(data)}
  {}

  const Data& data() const { return data_; }
  Data& data() { return data_; }

private:
  Data data_;
};

//! Add a `data` member to all tree nodes
template<class Data, class Tree,
  std::enable_if_t<std::is_default_constructible_v<Data>, int> = 0>
auto attachDataToTree(const Tree& tree)
{
  return transformTree(tree, [](auto&& node)
  {
    using DataNode = DataNodeMixin<std::decay_t<decltype(node)>,Data>;
    return DataNode{std::forward<decltype(node)>(node), Data{}};
  });
}

//! Add a `data` member to all tree nodes, constructed using the function `DataFactory(Node)`
template<class Tree, class DataFactory>
auto attachDataToTree(const Tree& tree, const DataFactory& dataFactory)
{
  return transformTree(tree, [dataFactory](auto const& node)
  {
    auto data = dataFactory(node);
    using DataNode = DataNodeMixin<std::decay_t<decltype(node)>,std::decay_t<decltype(data)> >;
    return DataNode{node, std::move(data)};
  });
}

} // namespace TypeTree2
} //namespace Dune

#endif // DUNE_TYPETREE2_TRANSFORMTREE_HH
