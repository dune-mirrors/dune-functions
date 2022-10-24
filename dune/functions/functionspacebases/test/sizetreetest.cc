// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/sizeinfo.hh>
#include <dune/typetree/treepath.hh>

#include "basisfactories.hh"

using namespace Dune;
using namespace Dune::Functions::BasisFactory;

template <class Node>
std::integral_constant<std::size_t, Node::size()> hybridSizeImpl (Node const& node, Dune::PriorityTag<2>)
{
  return {};
}

template <class Node>
std::size_t hybridSizeImpl (Node const& node, Dune::PriorityTag<1>)
{
  return node.size();
}

template <class Node>
auto hybridSize (Node const& node)
{
  return hybridSizeImpl(node, Dune::PriorityTag<3>{});
}

template <class T>
void printType ()
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
}

template <class T>
void printType(const T&)
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
}

// check at run time whether index is a valid child index
template <class Node, class Index>
std::true_type checkAccessIndex (Node const& node, Index i)
{
  assert(std::size_t(i) < node.size() && "Child index out of range");
  return {};
}

// check at compile time whether index is a valid index
template <class Node, std::size_t i>
std::bool_constant<(i < Node::size())> checkAccessIndex (Node const& node, index_constant<i>)
{
  static_assert(i < Node::size(), "Child index out of range");
  return {};
}

// finally return the node itself if no further indices are provided. Break condition
// for the recursion over the node childs.
template<class Node>
decltype(auto) accessImpl (Node&& node)
{
  return std::forward<Node>(node);
}

// recursively call `node[i]` with the given indices
template<class Node, class I0, class... I>
decltype(auto) accessImpl (Node&& node, I0 i0, [[maybe_unused]] I... i)
{
  auto valid = checkAccessIndex(node,i0);
  if constexpr (valid)
    return accessImpl(node[i0],i...);
  else
    return;
}

// forward to the impl methods by extracting the indices from the treepath
template<class Tree, class... Indices, std::size_t... i>
decltype(auto) access (Tree&& tree, [[maybe_unused]] TypeTree::HybridTreePath<Indices...> tp, std::index_sequence<i...>)
{
  return accessImpl(std::forward<Tree>(tree),TypeTree::treePathEntry<i>(tp)...);
}

// access a tree using a hybridTreePath
template<class Tree, class... Indices>
decltype(auto) access (Tree&& tree, TypeTree::HybridTreePath<Indices...> tp)
{
  return access(std::forward<Tree>(tree),tp,std::index_sequence_for<Indices...>{});
}

// convert a HybridTreePath into a ReservedVector
template <class SizePrefix, class... Indices>
auto sizePrefix (TypeTree::HybridTreePath<Indices...> tp)
{
  using value_type = typename SizePrefix::value_type;
  return Dune::unpackIntegerSequence([&](auto... i) {
    return SizePrefix{value_type(tp[i])...};
  }, std::index_sequence_for<Indices...>{});
}


template<class SizeTree, class SizeProvider, class... Indices>
void checkSize (TestSuite& test, const SizeTree& sizeTree, const SizeProvider& sizeProvider,
                TypeTree::HybridTreePath<Indices...> prefix)
{
  using SizePrefix = typename SizeProvider::SizePrefix;
  auto size1 = sizeProvider.size(sizePrefix<SizePrefix>(prefix));
  auto size2 = hybridSize(access(sizeTree, prefix));
  test.require(std::size_t(size1) == std::size_t(size2));

  if constexpr(sizeof...(Indices) < SizePrefix::max_size()) {
    Hybrid::forEach(Dune::range(size2), [&](auto i) {
      checkSize(test, sizeTree, sizeProvider, push_back(prefix,i));
    });
  }
}


int main (int argc, char *argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test;

  using Factories = Functions::BasisFactories<2>;
  Factories factories;
  Hybrid::forEach(Dune::StaticIntegralRange<std::size_t,Factories::num_bases>{},
  [&](auto i) {
    std::cout << std::size_t(i) << ") ";
    auto basis = factories.basis(i);
    printType(basis.preBasis());
    checkSize(test, basis.preBasis().sizeTree(), sizeInfo(basis), TypeTree::HybridTreePath<>{});
  });

  return test.exit();
}
