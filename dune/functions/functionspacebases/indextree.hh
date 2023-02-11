// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INDEXTREE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INDEXTREE_HH

#include <array>
#include <cassert>
#include <functional>
#include <vector>

#include <dune/common/filledarray.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/tuplevector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/basistags.hh>


/**
 * \file indextree.hh
 * \brief Lightweight representation of (hierarchic) size and blocking structure
 * The index-tree can be used to define types for data-structures, like vectors
 * or matrices, that can be accessed by the multi-indices provided by a basis.
 *
 * An index-tree encodes the dimensions of the index-space in a hierarchic basis.
 * This means especially, for each component of a multi-index it encodes the
 * possible range of indices.
 *
 * The structure of an index-tree is as follows
 * \code
  struct IndexTree
  {
    static constexpr bool isUniform = [true|false];     // Whether all children are identical
    static constexpr bool isTypeUniform = [true|false]; // Whether all children have the same type

    template<class Index>
    SubTree operator[](Index i) const;  // return the i-th sub-tree

    [static constexpr] std::size_t size() [const];  // return the number of sub-nodes
  };
 * \endcode
 *
 * The property `isUniform` specifies that all nodes of that tree are identical,
 * especially have the same size, the same properties and the same children.
 *
 * The property `isTypeUniform` is a bet less, it specifies that only the type
 * of all nodes of that tree are identical. This means that all static information
 * of the children are identical, e.g., the static size.
 *
 * With the `operator[]` the children can be accessed. Thereby, the `Index` type
 * is either an integral value or an `integral_constant` for non-uniform nodes.
 *
 * Size is either a static property, or a runtime value.
 **/

namespace Dune {
namespace Functions {

  //! Nested size-info that cannot be assigned one of the other index-tree types.
  struct UnknownIndexTree {};

  //! An index-tree of size zero
  struct EmptyIndexTree
  {
    static constexpr bool isUniform = true;
    static constexpr bool isTypeUniform = true;

    template <class Index>
    EmptyIndexTree operator[] (const Index&) const { return {}; }

    static constexpr std::size_t size () { return 0; }
  };


  //! Non-uniform index-tree with all sub-trees of different type
  template<class... SubTrees>
  struct StaticNonUniformIndexTree
      : private Dune::TupleVector<SubTrees...>
  {
    using Super = Dune::TupleVector<SubTrees...>;
    static constexpr bool isUniform = false;
    static constexpr bool isTypeUniform = false;

    template<class ST = std::tuple<SubTrees...>,
      std::enable_if_t<std::is_default_constructible_v<ST>, int> = 0>
    StaticNonUniformIndexTree ()
      : StaticNonUniformIndexTree{SubTrees{}...}
    {}

    explicit StaticNonUniformIndexTree (const SubTrees&... subTrees)
      : Super{subTrees...}
    {}

    using Super::operator[];
    using Super::size;
  };

  template<class SubTree0, class... SubTrees,
    std::enable_if_t<(sizeof...(SubTrees) > 0), int> = 0,
    std::enable_if_t<(...|| (not std::is_same_v<SubTree0, SubTrees>)), int> = 0>
  auto makeNonUniformIndexTree (SubTree0 subTree0, SubTrees... subTrees)
  {
    using IndexTree = StaticNonUniformIndexTree<SubTree0,SubTrees...>;
    return IndexTree{std::move(subTree0),std::move(subTrees)...};
  }


  //! Non-uniform index-tree with all sub-trees of different type but dynamic size
  //! NOTE, this cannot easily be implemented, maybe using type-erasure
  struct NonUniformIndexTree {};


  //! Non-uniform index-tree with all sub-trees of the same type and static size.
  template<class SubTree, std::size_t n>
  struct StaticTypeUniformIndexTree
      : private std::array<SubTree, n>
  {
    using Super = std::array<SubTree, n>;
    static constexpr bool isUniform = false;
    static constexpr bool isTypeUniform = true;

    template<class ST = SubTree,
      std::enable_if_t<std::is_default_constructible_v<ST>, int> = 0>
    StaticTypeUniformIndexTree ()
      : StaticTypeUniformIndexTree{SubTree{}}
    {}

    explicit StaticTypeUniformIndexTree (SubTree subTree)
      : Super{Dune::filledArray<n>(std::move(subTree))}
    {}

    StaticTypeUniformIndexTree (std::integral_constant<std::size_t,n>, SubTree subTree)
      : Super{Dune::filledArray<n>(std::move(subTree))}
    {}

    template <class... SubTrees,
      std::enable_if_t<(std::is_same_v<SubTrees,SubTree> &&...), int> = 0>
    explicit StaticTypeUniformIndexTree (SubTrees... subTrees)
      : Super{std::move(subTrees)...}
    {}

    using Super::operator[];

    static constexpr std::size_t size () { return n; }
  };

  template<class SubTree0, class... SubTrees,
    std::enable_if_t<(...&& std::is_same_v<SubTree0, SubTrees>), int> = 0>
  auto makeNonUniformIndexTree (SubTree0 subTree0, SubTrees... subTrees)
  {
    using IndexTree = StaticTypeUniformIndexTree<SubTree0,1+sizeof...(SubTrees)>;
    return IndexTree{std::move(subTree0),std::move(subTrees)...};
  }


  //! Non-uniform index-tree with all sub-trees of the same type and dynamic size.
  template<class SubTree>
  struct TypeUniformIndexTree
      : private std::vector<SubTree>
  {
    using Super = std::vector<SubTree>;
    static constexpr bool isUniform = false;
    static constexpr bool isTypeUniform = true;

    using Super::Super;
    using Super::operator[];
    using Super::size;
  };


  //! Uniform index-tree with static size.
  template<class SubTree, std::size_t n>
  struct StaticUniformIndexTree
  {
    static constexpr bool isUniform = true;
    static constexpr bool isTypeUniform = true;

    template<class ST = SubTree,
      std::enable_if_t<std::is_default_constructible_v<ST>, int> = 0>
    StaticUniformIndexTree ()
      : subTree_{}
    {}

    explicit StaticUniformIndexTree (SubTree subTree)
      : subTree_{std::move(subTree)}
    {}

    StaticUniformIndexTree (std::integral_constant<std::size_t,n>, SubTree subTree)
      : subTree_{std::move(subTree)}
    {}

    template <class Index>
    const SubTree& operator[] (const Index&) const { return subTree_; }

    static constexpr std::size_t size () { return n; }

  private:
    SubTree subTree_;
  };

  template <std::size_t n>
  using StaticFlatIndexTree = StaticUniformIndexTree<EmptyIndexTree,n>;

  template<class SubTree, std::size_t n>
  auto makeUniformIndexTree (std::integral_constant<std::size_t,n>, SubTree subTree)
  {
    return StaticUniformIndexTree<SubTree,n>{std::move(subTree)};
  }


  //! Uniform index-tree with dynamic size.
  template<class SubTree>
  struct UniformIndexTree
  {
    static constexpr bool isUniform = true;
    static constexpr bool isTypeUniform = true;

    template<class ST = SubTree,
      std::enable_if_t<std::is_default_constructible_v<ST>, int> = 0>
    UniformIndexTree ()
      : size_{0}
      , subTree_{}
    {}

    template<class ST = SubTree,
      std::enable_if_t<std::is_default_constructible_v<ST>, int> = 0>
    explicit UniformIndexTree (std::size_t size)
      : size_{size}
      , subTree_{}
    {}

    UniformIndexTree (std::size_t size, SubTree subTree)
      : size_{size}
      , subTree_{std::move(subTree)}
    {}

    template <class Index>
    const SubTree& operator[] (const Index&) const { return subTree_; }

    std::size_t size () const { return size_; }

  private:
    std::size_t size_;
    SubTree subTree_;
  };

  using FlatIndexTree = UniformIndexTree<EmptyIndexTree>;

  template<class SubTree>
  auto makeUniformIndexTree (std::size_t n, SubTree subTree)
  {
    return UniformIndexTree<SubTree>{n,std::move(subTree)};
  }


  namespace Impl {

  // -----------------------------------------------------------------------------------------------
  // Some utilities for generating index-trees
  // -----------------------------------------------------------------------------------------------

  class FlatIndexAccess
  {
    struct Increment
    {
      template<class T>
      constexpr auto operator() (T i) const { return i+1; }
    };

    static constexpr auto lt_ = Hybrid::hybridFunctor(std::less<>{});
    static constexpr auto incr_ = Hybrid::hybridFunctor(Increment{});
    static constexpr auto minus_ = Hybrid::hybridFunctor(std::minus<>{});

    template<class FlatIndex, class OuterOffsetIndex, class IMS>
    static auto getEntry (UnknownIndexTree, FlatIndex, OuterOffsetIndex, IMS)
    {
      return UnknownIndexTree{};
    }

    template<class FlatIndex, class OuterOffsetIndex, class IMS>
    static auto getEntry (EmptyIndexTree, FlatIndex, OuterOffsetIndex, IMS)
    {
      return UnknownIndexTree{};
    }

    // overload for flat-lexicographic index-merging strategy
    template<class IndexTree, class FlatIndex, class OuterOffsetIndex>
    static auto getEntry (const IndexTree& tree, FlatIndex i, OuterOffsetIndex o,
                          Dune::Functions::BasisFactory::FlatLexicographic ims)
    {
      const auto treeSize = Hybrid::size(tree);
      const auto outerOffsetSize = Hybrid::size(tree[o]);
      return Hybrid::ifElse(lt_(i,outerOffsetSize), // i < outerOffsetSize
        [&](auto id) {
          return id(tree)[o][i];
        },
        [&](auto id) {
          assert((o+1) < treeSize && i >= outerOffsetSize);
          return getEntry(id(tree), minus_(i,outerOffsetSize), incr_(o), ims);
        });
    }

    // overload for flat-interleaved index-merging strategy
    template<class IndexTree, class FlatIndex, class InnerIndex>
    static auto getEntry (const IndexTree& tree, FlatIndex i, InnerIndex o,
                          Dune::Functions::BasisFactory::FlatInterleaved ims)
    {
      const auto rootSize = Hybrid::size(tree);
      return Hybrid::ifElse(lt_(i,rootSize), // i < rootSize
        [&](auto id) {
          return id(tree)[i][o];
        },
        [&](auto id) {
          assert(i >= rootSize);
          return getEntry(id(tree), minus_(i,rootSize), incr_(o), ims);
        });
    }

  public:
    // Return the sub-node in a given index-tree associated to a flat index
    // after merging the nodes with a given index-merging strategy
    template<class IMS, class IndexTree, class FlatIndex>
    static auto getEntry (const IndexTree& tree, FlatIndex i)
    {
      return getEntry(tree,i,Indices::_0,IMS{});
    }
  };

  template<class IMS, class Size, class IndexTree, bool allUniform, bool allTypeUniform>
  auto mergeIndexTreesImpl (Size size, const IndexTree& tree,
                            std::bool_constant<allUniform>,
                            std::bool_constant<allTypeUniform>)
  {
    auto child = [&](auto ii) { return FlatIndexAccess::getEntry<IMS>(tree,ii); };

    if constexpr(allUniform)
      return makeUniformIndexTree(size, tree[Indices::_0][Indices::_0]);
    else if constexpr(IsIntegralConstant<Size>::value)
      return unpackIntegerSequence([&](auto... ii) {
        return makeNonUniformIndexTree(child(ii)...);
      }, std::make_index_sequence<Size::value>{});
    else if constexpr(allTypeUniform) {
      TypeUniformIndexTree<decltype(child(0))> result(size);
      for (std::size_t i = 0; i < size; ++i)
        result[i] = std::move(child(i));
      return result;
    }
    else {
      DUNE_THROW(Dune::NotImplemented,
        "Merging of dynamic non-uniform index-trees not implemented");
      return UnknownIndexTree{};
    }
  }


  // Overload for zero index-trees, return an unknown tree.
  inline auto mergeIndexTrees ()
  {
    return UnknownIndexTree{};
  }

  // Overload for one index-tree, return the tree itself.
  template<class IT>
  const auto& mergeIndexTrees (const IT& indexTree)
  {
    return indexTree;
  }

  // Generate a sum of a different index-trees
  template<class IMS, class... IT,
    std::enable_if_t<(sizeof...(IT) > 1), int> = 0>
  auto mergeIndexTrees (const IT&... indexTrees)
  {
    auto sumSizes = Hybrid::plus(Hybrid::size(indexTrees)...);
    return mergeIndexTreesImpl<IMS>(sumSizes, makeNonUniformIndexTree(indexTrees...),
      std::bool_constant<(... && IT::isUniform)>{},
      std::bool_constant<(... && IT::isTypeUniform)>{}
    );
  }

  // Generate a sum of an index-tree consisting of `n` identical `indexTree`s
  template<class IMS, class Size, class IT,
    std::enable_if_t<std::is_convertible_v<Size,std::size_t>, int> = 0>
  auto mergeIndexTrees (Size n, const IT& indexTree)
  {
    auto multiplies = Hybrid::hybridFunctor(std::multiplies<>{});
    auto sumSizes = multiplies(Hybrid::size(indexTree), n);
    return mergeIndexTreesImpl<IMS>(sumSizes, makeUniformIndexTree(n,indexTree),
      std::bool_constant<IT::isUniform>{},
      std::bool_constant<IT::isTypeUniform>{}
    );
  }

  // Generate a sum of an index-tree consisting of `n` identical `indexTree`s
  template<std::size_t n, class IMS, class IT>
  auto mergeIndexTrees (const IT& indexTree)
  {
    return mergeIndexTrees<IMS>(index_constant<n>{}, indexTree);
  }


  /*
   * Append a size to all children of an index-tree
   *
   * Transforming index-trees by appending a size, e.g.,
   * a blocked-interleaved index-merging strategy in a power-basis.
   *
   * append( FlatIndexTree it, size ) -> UniformIndexTree( it.size(), FlatIndexTree(size) )
   * append( IndexTree(child...), size ) -> IndexTree( append(child, size)... )
   */

  // Append the size `s` at the inner-most node of the tree
  template<class Size>
  auto appendToIndexTree (EmptyIndexTree, Size s)
  {
    return makeUniformIndexTree(s, EmptyIndexTree{});
  }

  // Append the size `s` at the inner-most node of the tree
  template<class IT, class Size>
  auto appendToIndexTree (const IT& indexTree, Size s)
  {
    auto child = [&](auto ii) { return appendToIndexTree(indexTree[ii], s); };

    if constexpr(IT::isUniform)
      return makeUniformIndexTree(Hybrid::size(indexTree), child(Indices::_0));
    else if constexpr(HasStaticSize_v<IT>)
      return unpackIntegerSequence([&](auto... ii) {
        return makeNonUniformIndexTree(child(ii)...);
      }, std::make_index_sequence<std::size_t(IT::size())>());
    else if constexpr(IT::isTypeUniform) {
      TypeUniformIndexTree<decltype(child(0))> result(indexTree.size());
      for (std::size_t i = 0; i < indexTree.size(); ++i)
        result[i] = std::move(child(i));
      return result;
    }
    else {
      DUNE_THROW(Dune::NotImplemented,
        "Merging of dynamic non-uniform index-trees not implemented");
      return UnknownIndexTree{};
    }
  }

  // Append the size `s` at the inner-most node of the tree
  template<std::size_t s, class IT>
  auto appendToIndexTree (const IT& indexTree)
  {
    return appendToIndexTree(indexTree, index_constant<s>{});
  }

  } // end namespace Impl

}} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INDEXTREE_HH
