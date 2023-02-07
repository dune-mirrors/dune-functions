// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/classname.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/typetree/treepath.hh>

#include "basisfactories.hh"

using namespace Dune;

namespace TestImpl {

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
// for the recursion over the node children.
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

// forward to the impl methods by extracting the indices from the tree path
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


template<class T, class Value, class Branches>
constexpr void hybridSwitchCases(IntegralRange<T> cases,
                                 const Value& value, Branches&& branches)
{
  branches(value);
}

template<class T, T to, T from, class Value, class Branches>
constexpr void hybridSwitchCases(StaticIntegralRange<T,to,from> cases,
                                 const Value& value, Branches&& branches)
{
  using seq = typename decltype(cases)::integer_sequence;
  Dune::Hybrid::switchCases(seq{}, value, std::forward<Branches>(branches), []() {});
}

} // end namespace TestImpl

// ---------------------------------------------------

// check that the sizes of an index-tree correspond to the sizes provided by the
// basis (size-provider) directly.
template<class IndexTree, class SizeProvider, class... Indices>
void checkSize (TestSuite& test, const IndexTree& indexTree, const SizeProvider& sizeProvider,
                TypeTree::HybridTreePath<Indices...> prefix)
{
  using SizePrefix = typename SizeProvider::SizePrefix;
  auto size1 = sizeProvider.size(TestImpl::sizePrefix<SizePrefix>(prefix));
  auto size2 = Dune::Hybrid::size(TestImpl::access(indexTree, prefix));
  test.require(std::size_t(size1) == std::size_t(size2));

  if constexpr(sizeof...(Indices) < SizePrefix::max_size()) {
    Hybrid::forEach(Dune::range(size2), [&](auto i) {
      checkSize(test, indexTree, sizeProvider, push_back(prefix,i));
    });
  }
}

// check a specific multi-index by traversing all its components and the index-tree simultaneously
template <class IndexTree, class MultiIndex>
void checkMultiIndex (TestSuite& test, const IndexTree& indexTree,
                      const MultiIndex& mi, std::size_t j = 0)
{
  if (j < mi.size()) {
    test.require(mi[j] < indexTree.size());
    auto size = Dune::Hybrid::size(indexTree);
    TestImpl::hybridSwitchCases(Dune::range(size), mi[j],
      [&](auto jj) { checkMultiIndex(test,indexTree[jj],mi,j+1); });
  }
}

// check that all multi-indices of a global basis are within the range of the index-tree
template <class IndexTree, class Basis>
void checkMultiIndices (TestSuite& test, const IndexTree& indexTree, const Basis& basis)
{
  auto localView = basis.localView();
  for (auto const& e : elements(basis.gridView()))
  {
    localView.bind(e);
    for (std::size_t i = 0; i < localView.size(); ++i) {
      auto mi = localView.index(i);
      checkMultiIndex(test,indexTree,mi);
    }
  }
}

int main (int argc, char *argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test;

  using Grid = Dune::YaspGrid<2>;
  Grid grid({1.0, 1.0}, {2, 2});

  using namespace Dune::Functions::BasisFactory;
  Hybrid::forEach(Dune::StaticIntegralRange<std::size_t,BasisFactories<2>::num_bases>{},
  [&](auto i) {
    std::cout << std::size_t(i) << ") ";
    auto basis = makeBasis(grid.leafGridView(), BasisFactories<2>::basis(i));
    std::cout << Dune::className(basis.preBasis()) << std::endl;
    checkSize(test, basis.preBasis().indexTree(), basis, TypeTree::HybridTreePath<>{});
    checkMultiIndices(test, basis.preBasis().indexTree(), basis);
  });

  return test.exit();
}
