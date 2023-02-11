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

template <class SizeProvider, class... Indices>
void checkSize (TestSuite& test, const Dune::Functions::UnknownIndexTree& indexTree,
                const SizeProvider& sizeProvider, TypeTree::HybridTreePath<Indices...> prefix)
{
  test.require(false) << "Got an UnknownIndexTree!";
}

// check that the sizes of an index-tree correspond to the sizes provided by the
// basis (size-provider) directly.
template<class IndexTree, class SizeProvider, class... Indices>
void checkSize (TestSuite& test, const IndexTree& indexTree, const SizeProvider& sizeProvider,
                TypeTree::HybridTreePath<Indices...> prefix)
{
  using SizePrefix = typename SizeProvider::SizePrefix;
  auto size1 = sizeProvider.size(TestImpl::sizePrefix<SizePrefix>(prefix));
  auto size2 = Dune::Hybrid::size(TestImpl::access(indexTree, prefix));
  test.require(std::size_t(size1) == std::size_t(size2), "size1 == size2");

  if constexpr(sizeof...(Indices) < SizePrefix::max_size()) {
    Hybrid::forEach(Dune::range(size2), [&](auto i) {
      checkSize(test, indexTree, sizeProvider, push_back(prefix,i));
    });
  }
}

template <class MultiIndex>
void checkMultiIndex (TestSuite& test, const Dune::Functions::UnknownIndexTree& indexTree,
                      const MultiIndex& mi, std::size_t j = 0)
{
  test.require(false) << "Got an UnknownIndexTree!";
}

// check a specific multi-index by traversing all its components and the index-tree simultaneously
template <class IndexTree, class MultiIndex>
void checkMultiIndex (TestSuite& test, const IndexTree& indexTree,
                      const MultiIndex& mi, std::size_t j = 0)
{
  if (j < mi.size()) {
    test.check(mi[j] < indexTree.size(), "mi[j] < indexTree.size");
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


void testMergeIndexTrees (TestSuite& test)
{
  using namespace Dune::Functions;
  using namespace Dune::Functions::BasisFactory;
  using FI = FlatInterleaved;
  using FL = FlatLexicographic;

  using Dune::Functions::Impl::mergeIndexTrees;

  { // merge no trees
    auto it1 = mergeIndexTrees<FI>();
    static_assert(std::is_same_v<decltype(it1), UnknownIndexTree>);
  }

  { // merge single tree
    auto it1 = mergeIndexTrees<FI>(StaticFlatIndexTree<4>{});
    static_assert(std::is_same_v<decltype(it1), StaticFlatIndexTree<4>>);

    auto it2 = mergeIndexTrees<FI>(FlatIndexTree{4});
    static_assert(std::is_same_v<decltype(it2), FlatIndexTree>);
  }

  { // merge flat index-trees
    auto it1a = mergeIndexTrees<FI>(StaticFlatIndexTree<3>{},StaticFlatIndexTree<3>{});
    static_assert(std::is_same_v<decltype(it1a), StaticFlatIndexTree<6>>);

    auto it1b = mergeIndexTrees<FL>(StaticFlatIndexTree<3>{},StaticFlatIndexTree<4>{});
    static_assert(std::is_same_v<decltype(it1b), StaticFlatIndexTree<7>>);

    auto it2a = mergeIndexTrees<FI>(FlatIndexTree{2},FlatIndexTree{2});
    static_assert(std::is_same_v<decltype(it2a), FlatIndexTree>);
    test.check(it2a.size() == 4, "merge dynamic flat index-trees");

    auto it2b = mergeIndexTrees<FL>(FlatIndexTree{2},FlatIndexTree{3});
    static_assert(std::is_same_v<decltype(it2b), FlatIndexTree>);
    test.check(it2b.size() == 5, "merge dynamic flat index-trees");

    auto it3a = mergeIndexTrees<FI>(FlatIndexTree{1},StaticFlatIndexTree<1>{});
    static_assert(std::is_same_v<decltype(it3a), FlatIndexTree>);
    test.check(it3a.size() == 2, "merge mixed flat index-trees");

    auto it3b = mergeIndexTrees<FL>(FlatIndexTree{1},StaticFlatIndexTree<2>{});
    static_assert(std::is_same_v<decltype(it3b), FlatIndexTree>);
    test.check(it3b.size() == 3, "merge mixed flat index-trees");

    auto it4a = mergeIndexTrees<4,FI>(FlatIndexTree{2});
    static_assert(std::is_same_v<decltype(it4a), FlatIndexTree>);
    test.check(it4a.size() == 8, "merge 4 dynamic flat index-trees");

    auto it4b = mergeIndexTrees<4,FL>(FlatIndexTree{2});
    static_assert(std::is_same_v<decltype(it4b), FlatIndexTree>);
    test.check(it4b.size() == 8, "merge 4 dynamic flat index-trees");
  }

  {
    auto it1a = mergeIndexTrees<FL>(UniformIndexTree{2,FlatIndexTree{3}},
                                    UniformIndexTree{3,FlatIndexTree{2}});
    static_assert(std::is_same_v<decltype(it1a), TypeUniformIndexTree<FlatIndexTree>>);
    test.check(it1a.size() == 5, "merge uniform index-trees FlatLexicographic");
    test.check(it1a[0].size() == 3, "merge uniform index-trees, tree[0].size");
    test.check(it1a[1].size() == 3, "merge uniform index-trees, tree[1].size");
    test.check(it1a[2].size() == 2, "merge uniform index-trees, tree[2].size");
    test.check(it1a[3].size() == 2, "merge uniform index-trees, tree[3].size");
    test.check(it1a[4].size() == 2, "merge uniform index-trees, tree[4].size");

    auto it1b = mergeIndexTrees<FI>(UniformIndexTree{2,FlatIndexTree{3}},
                                    UniformIndexTree{2,FlatIndexTree{2}});
    static_assert(std::is_same_v<decltype(it1b), TypeUniformIndexTree<FlatIndexTree>>);
    test.check(it1b.size() == 4, "merge uniform index-trees FlatInterleaved");
    test.check(it1b[0].size() == 3, "merge uniform index-trees, tree[0].size");
    test.check(it1b[1].size() == 2, "merge uniform index-trees, tree[1].size");
    test.check(it1b[2].size() == 3, "merge uniform index-trees, tree[2].size");
    test.check(it1b[3].size() == 2, "merge uniform index-trees, tree[3].size");
  }
}

int main (int argc, char *argv[])
{
  MPIHelper::instance(argc, argv);

  TestSuite test;

  using Grid = Dune::YaspGrid<2>;
  Grid grid({1.0, 1.0}, {2, 2});

  using namespace Dune::Functions::BasisFactory;
  Hybrid::forEach(Dune::StaticIntegralRange<std::size_t,BasisFactories::size>{},
  [&](auto i) {
    std::cout << std::size_t(i) << ") ";
    auto basis = makeBasis(grid.leafGridView(), BasisFactories::basis(i));
    std::cout << Dune::className(basis.preBasis()) << std::endl;
    checkSize(test, basis.preBasis().indexTree(), basis, TypeTree::HybridTreePath<>{});
    checkMultiIndices(test, basis.preBasis().indexTree(), basis);
  });

  testMergeIndexTrees(test);

  return test.exit();
}
