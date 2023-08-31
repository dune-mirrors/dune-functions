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

// template <class SizeProvider, class... Indices>
// void checkSize (TestSuite& test, const Dune::Functions::UnknownIndexTree& containerDescriptor,
//                 const SizeProvider& sizeProvider, TypeTree::HybridTreePath<Indices...> prefix)
// {
//   test.require(false) << "Got an UnknownIndexTree!";
// }

// check that the sizes of an container descriptor correspond to the sizes provided by the
// basis (size-provider) directly.
template<class ContainerDescriptor, class SizeProvider, class... Indices>
void checkSize (TestSuite& test, const ContainerDescriptor& containerDescriptor, const SizeProvider& sizeProvider, TypeTree::HybridTreePath<Indices...> prefix)
{
  using SizePrefix = typename SizeProvider::SizePrefix;
  auto size1 = sizeProvider.size(TestImpl::sizePrefix<SizePrefix>(prefix));
  auto size2 = Dune::Hybrid::size(TestImpl::access(containerDescriptor, prefix));
  test.require(std::size_t(size1) == std::size_t(size2), "size1 == size2");

  if constexpr(sizeof...(Indices) < SizePrefix::max_size()) {
    Hybrid::forEach(Dune::range(size2), [&](auto i) {
      checkSize(test, containerDescriptor, sizeProvider, push_back(prefix,i));
    });
  }
}

// template <class MultiIndex>
// void checkMultiIndex (TestSuite& test, const Dune::Functions::UnknownIndexTree& containerDescriptor,
//                       const MultiIndex& mi, std::size_t j = 0)
// {
//   test.require(false) << "Got an UnknownIndexTree!";
// }

// check a specific multi-index by traversing all its components and the container descriptor
//  simultaneously
template <class ContainerDescriptor, class MultiIndex>
void checkMultiIndex (TestSuite& test, const ContainerDescriptor& containerDescriptor,
                      const MultiIndex& mi, std::size_t j = 0)
{
  if (j < mi.size()) {
    test.check(mi[j] < containerDescriptor.size(), "mi[j] < containerDescriptor.size");
    auto size = Dune::Hybrid::size(containerDescriptor);
    TestImpl::hybridSwitchCases(Dune::range(size), mi[j],
      [&](auto jj) { checkMultiIndex(test,containerDescriptor[jj],mi,j+1); });
  }
}

// check that all multi-indices of a global basis are within the range of the container descriptor
template <class ContainerDescriptor, class Basis>
void checkMultiIndices (TestSuite& test, const ContainerDescriptor& containerDescriptor, const Basis& basis)
{
  auto localView = basis.localView();
  for (auto const& e : elements(basis.gridView()))
  {
    localView.bind(e);
    for (std::size_t i = 0; i < localView.size(); ++i) {
      auto mi = localView.index(i);
      checkMultiIndex(test,containerDescriptor,mi);
    }
  }
}


void testMergeTrees (TestSuite& test)
{
  using namespace Dune::Functions;
  using namespace Dune::Functions::BasisFactory;
  using FI = FlatInterleaved;
  using FL = FlatLexicographic;

  using Dune::Functions::ContainerDescriptors::Impl::mergeTrees;
  using Dune::Functions::ContainerDescriptors::Impl::mergeIdenticalTrees;

  { // merge no trees
    using CD1 = decltype(mergeTrees<FI>());
    static_assert(std::is_void_v<CD1>);
  }

  { // merge single tree
    using CD1 = decltype(mergeTrees<FI>(ContainerDescriptors::FlatArray<4>{}));
    static_assert(std::is_same_v<CD1, ContainerDescriptors::FlatArray<4>>);

    using CD2 = decltype(mergeTrees<FI>(ContainerDescriptors::FlatVector{4}));
    static_assert(std::is_same_v<CD2, ContainerDescriptors::FlatVector>);
  }

  { // merge flat trees
    using CD1a = decltype(mergeTrees<FI>(ContainerDescriptors::FlatArray<3>{},ContainerDescriptors::FlatArray<3>{}));
    static_assert(std::is_same_v<CD1a, ContainerDescriptors::FlatArray<6>>);

    using CD1b = decltype(mergeTrees<FL>(ContainerDescriptors::FlatArray<3>{},ContainerDescriptors::FlatArray<4>{}));
    static_assert(std::is_same_v<CD1b, ContainerDescriptors::FlatArray<7>>);

    auto cd2a = mergeTrees<FI>(ContainerDescriptors::FlatVector{2},ContainerDescriptors::FlatVector{2});
    static_assert(std::is_same_v<decltype(cd2a), ContainerDescriptors::FlatVector>);
    test.check(cd2a.size() == 4, "merge dynamic flat container descriptor");

    auto cd2b = mergeTrees<FL>(ContainerDescriptors::FlatVector{2},ContainerDescriptors::FlatVector{3});
    static_assert(std::is_same_v<decltype(cd2b), ContainerDescriptors::FlatVector>);
    test.check(cd2b.size() == 5, "merge dynamic flat container descriptor");

    auto cd3a = mergeTrees<FI>(ContainerDescriptors::FlatVector{1},ContainerDescriptors::FlatArray<1>{});
    static_assert(std::is_same_v<decltype(cd3a), ContainerDescriptors::FlatVector>);
    test.check(cd3a.size() == 2, "merge mixed flat container descriptor");

    auto cd3b = mergeTrees<FL>(ContainerDescriptors::FlatVector{1},ContainerDescriptors::FlatArray<2>{});
    static_assert(std::is_same_v<decltype(cd3b), ContainerDescriptors::FlatVector>);
    test.check(cd3b.size() == 3, "merge mixed flat container descriptor");

    auto cd4a = mergeIdenticalTrees<4,FI>(ContainerDescriptors::FlatVector{2});
    static_assert(std::is_same_v<decltype(cd4a), ContainerDescriptors::FlatVector>);
    test.check(cd4a.size() == 8, "merge 4 dynamic flat container descriptor");

    auto cd4b = mergeIdenticalTrees<4,FL>(ContainerDescriptors::FlatVector{2});
    static_assert(std::is_same_v<decltype(cd4b), ContainerDescriptors::FlatVector>);
    test.check(cd4b.size() == 8, "merge 4 dynamic flat container descriptor");
  }

  {
    auto cd1a = mergeTrees<FL>(
      ContainerDescriptors::UniformVector{2,ContainerDescriptors::FlatVector{3}},
      ContainerDescriptors::UniformVector{3,ContainerDescriptors::FlatVector{2}});
    static_assert(std::is_same_v<decltype(cd1a), ContainerDescriptors::Vector<ContainerDescriptors::FlatVector>>);
    test.check(cd1a.size() == 5, "merge uniform container descriptor FlatLexicographic");
    test.check(cd1a[0].size() == 3, "merge uniform container descriptor, tree[0].size");
    test.check(cd1a[1].size() == 3, "merge uniform container descriptor, tree[1].size");
    test.check(cd1a[2].size() == 2, "merge uniform container descriptor, tree[2].size");
    test.check(cd1a[3].size() == 2, "merge uniform container descriptor, tree[3].size");
    test.check(cd1a[4].size() == 2, "merge uniform container descriptor, tree[4].size");

    auto cd1b = mergeTrees<FI>(
      ContainerDescriptors::UniformVector{2,ContainerDescriptors::FlatVector{3}},
      ContainerDescriptors::UniformVector{2,ContainerDescriptors::FlatVector{2}});
    static_assert(std::is_same_v<decltype(cd1b), ContainerDescriptors::Vector<ContainerDescriptors::FlatVector>>);
    test.check(cd1b.size() == 4, "merge uniform container descriptor FlatInterleaved");
    test.check(cd1b[0].size() == 3, "merge uniform container descriptor, tree[0].size");
    test.check(cd1b[1].size() == 2, "merge uniform container descriptor, tree[1].size");
    test.check(cd1b[2].size() == 3, "merge uniform container descriptor, tree[2].size");
    test.check(cd1b[3].size() == 2, "merge uniform container descriptor, tree[3].size");
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
    checkSize(test, basis.preBasis().containerDescriptor(), basis, TypeTree::HybridTreePath<>{});
    checkMultiIndices(test, basis.preBasis().containerDescriptor(), basis);
  });

  testMergeTrees(test);

  return test.exit();
}
