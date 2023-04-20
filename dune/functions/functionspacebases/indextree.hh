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
 * \brief Lightweight representation of (hierarchical) size and block structure.
 *
 * The index-tree can be used to define types for data-structures, such as vectors
 * or matrices, that can be accessed by the multi-indices provided by a basis.
 *
 * An index-tree encodes the dimensions of the index-space in a hierarchical basis.
 * In particular, for each component of a multi-index, it encodes the
 * possible range of indices.
 *
 * The structure of an index-tree is as follows:
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
 * The `isUniform` property specifies that all nodes of this tree are identical,
 * specifically that they have the same size, the same properties, and the same children.
 *
 * The `isTypeUniform` property is a bit less restrictive. It specifies that only the type
 * of all nodes of this tree are identical. This means that all static information
 * of the children is identical, e.g., the static size.
 *
 * With the `operator[]` you can access the children. The `Index` type is either
 * an integral value or an `integral_constant` for non-uniform nodes.
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

  //! Generate a StaticNonUniformIndexTree in case the sub-trees are not all of the same type.
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

  //! Generate a StaticTypeUniformIndexTree in case the sub-trees are all of the same type.
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

  //! Generate a StaticUniformIndexTree in case the size is a static constant
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

  //! Generate a UniformIndexTree in case the size is a dynamic value
  template<class SubTree>
  auto makeUniformIndexTree (std::size_t n, SubTree subTree)
  {
    return UniformIndexTree<SubTree>{n,std::move(subTree)};
  }


  namespace Impl {

  // ---------------------------------------------------------------------------
  // Some utilities for generating index-trees
  // ---------------------------------------------------------------------------

  // Hybrid assert for static error messages.
  // Overload that checks the static constant bool.
  template <bool b>
  void hybridAssert(std::bool_constant<b>)
  {
    static_assert(b);
  }

  // Hybrid assert for dynamic error messages
  // Overload that checks the boolean value.
  inline void hybridAssert(bool b)
  {
    assert(b);
  }

  // Constexpr functor that is used to define a Hybrid::HybridFunctor
  struct LogicalAnd
  {
    template<class... T>
    constexpr auto operator() (T... b) const { return (...&& b); }
  };

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

    // overload for UnknownIndexTrees. Cannot access any component.
    template<class FlatIndex, class OuterOffsetIndex, class IMS>
    static auto getEntry (UnknownIndexTree, FlatIndex, OuterOffsetIndex, IMS)
    {
      return UnknownIndexTree{};
    }

    // overload for UnknownIndexTrees. Cannot access any component.
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
          return Hybrid::ifElse(lt_(incr_(o),treeSize), // (o+1) < treeSize
            [&](auto id_) {
              return getEntry(id_(id(tree)), minus_(i,outerOffsetSize), incr_(o), ims);
            },
            [&](auto id_) {
              // fallback condition to fix a return type
              return id_(id(tree))[Indices::_0][Indices::_0];
            });
        });
    }

    // overload for flat-interleaved index-merging strategy
    template<class IndexTree, class FlatIndex, class InnerIndex>
    static auto getEntry (const IndexTree& tree, FlatIndex i, InnerIndex o,
                          Dune::Functions::BasisFactory::FlatInterleaved ims)
    {
      const auto treeSize = Hybrid::size(tree);
      return Hybrid::ifElse(lt_(i,treeSize), // i < treeSize
        [&](auto id) {
          return id(tree)[i][o];
        },
        [&](auto id) {
          assert(i >= treeSize && o < treeSize);
          return Hybrid::ifElse(lt_(o,treeSize), // o < treeSize
            [&](auto id_) {
              return getEntry(id_(id(tree)), minus_(i,treeSize), incr_(o), ims);
            },
            [&](auto id_) {
              // fallback condition to fix a return type
              return id_(id(tree))[Indices::_0][Indices::_0];
            });
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


  /*
   * Generic implementation of the index-tree merging. The trees to merge are
   * collected into an artificial super index-tree `tree`. The number of nodes
   * after merging the sub-trees of `tree` is given by the parameter `size`. And
   * the tree properties as `isUniform` and `isTypeUniform`.
   *
   * \tparam IMS  The index-merging strategy used as strategy for merging the tree
   *    nodes, either BasisFactory::FlatLexicographic or BasisFactory::FlatInterleaved.
   *
   * \param tree  Collection of nodes to be merged
   * \param size  The size of the resulting index-tree
   * \param isUniform  The resulting index-tree has the property `isUniform`
   * \param isTypeUniform  The resulting index-tree has the property `isTypeUniform`
   */
  template<class IMS, class Tree, class Size, bool isUniform, bool isTypeUniform>
  auto mergeIndexTreesImpl (const Tree& tree, Size size,
                            std::bool_constant<isUniform>,
                            std::bool_constant<isTypeUniform>)
  {
    // access one of the sub-nodes of the nodes in `tree` by a flat index `ii`
    auto child = [&](auto ii) { return FlatIndexAccess::getEntry<IMS>(tree,ii); };
    using Child00 = std::decay_t<decltype(tree[Indices::_0][Indices::_0])>;

    if constexpr(isUniform)
      return makeUniformIndexTree(size, tree[Indices::_0][Indices::_0]);
    else if constexpr(IsIntegralConstant<Size>::value)
      return unpackIntegerSequence([&](auto... ii) {
        return makeNonUniformIndexTree(child(ii)...);
      }, std::make_index_sequence<Size::value>{});
    else if constexpr(isTypeUniform) {
      TypeUniformIndexTree<Child00> result(size);
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
  template<class IMS>
  inline auto mergeIndexTrees ()
  {
    return UnknownIndexTree{};
  }

  // Overload for one index-tree, return the tree itself.
  template<class IMS,class IT>
  const auto& mergeIndexTrees (const IT& indexTree)
  {
    return indexTree;
  }

  // Merge a variadic list of index-trees
  template<class IMS, class IT0, class... IT,
    std::enable_if_t<(sizeof...(IT) > 0), int> = 0>
  auto mergeIndexTrees (const IT0& indexTree0, const IT&... indexTrees)
  {
    auto and_ = Hybrid::hybridFunctor(LogicalAnd{});
    auto or_ = Hybrid::hybridFunctor(std::logical_or<>{});
    auto equal_ = Hybrid::hybridFunctor(std::equal_to<>{});

    using IT00 = std::decay_t<decltype(indexTree0[Indices::_0])>;

    // collect some properties of the passed index-trees
    const bool allUniform = (IT0::isUniform &&...&& IT::isUniform);
    const bool allTypeUniform = (IT0::isTypeUniform &&...&& IT::isTypeUniform);
    const bool allSubTypeUniform = (...&&
      std::is_same_v<IT00, std::decay_t<decltype(indexTrees[Indices::_0])>>);
    const auto allSameSize = and_(equal_(Hybrid::size(indexTree0), Hybrid::size(indexTrees))...);

    // The resulting index-tree can only be uniform if we can deduce at compile-time that
    // all nodes of all index-trees are identical. This is only possible in some cases,
    // e.g. if all nodes are `EmptyIndexTree`s.
    const bool allSubTypeEmpty = (std::is_same_v<IT00, EmptyIndexTree> &&...&&
      std::is_same_v<std::decay_t<decltype(indexTrees[Indices::_0])>, EmptyIndexTree>);

    // The merge of variadic index-trees currently is only implemented with lexicographic merge
    // or interleaved merge if all index-trees have the same size
    const bool isFlatLexicographic = std::is_same_v<IMS,BasisFactory::FlatLexicographic>;
    hybridAssert(or_(std::bool_constant<isFlatLexicographic>{}, allSameSize));

    auto sumSizes = Hybrid::plus(Hybrid::size(indexTree0),Hybrid::size(indexTrees)...);
    return mergeIndexTreesImpl<IMS>(makeNonUniformIndexTree(indexTree0,indexTrees...),
      sumSizes,
      std::bool_constant<allUniform && allSubTypeEmpty>{},
      std::bool_constant<allTypeUniform && allSubTypeUniform>{}
    );
  }

  // Merge `n` identical `indexTree`s
  template<class IMS, class Size, class IT,
    std::enable_if_t<std::is_convertible_v<Size,std::size_t>, int> = 0>
  auto mergeIndexTrees (Size n, const IT& indexTree)
  {
    auto multiplies = Hybrid::hybridFunctor(std::multiplies<>{});
    auto sumSizes = multiplies(Hybrid::size(indexTree), n);
    return mergeIndexTreesImpl<IMS>(makeUniformIndexTree(n,indexTree),
      sumSizes,
      std::bool_constant<IT::isUniform>{},
      std::bool_constant<IT::isTypeUniform>{}
    );
  }

  // Merge `n` identical `indexTree`s, redirects to the function with 2 arguments
  template<std::size_t n, class IMS, class IT>
  auto mergeIndexTrees (const IT& indexTree)
  {
    return mergeIndexTrees<IMS>(index_constant<n>{}, indexTree);
  }



  // Append the size `s` at the inner-most node of the tree
  template<class Size>
  auto appendToIndexTree (EmptyIndexTree, Size s)
  {
    return makeUniformIndexTree(s, EmptyIndexTree{});
  }

  /*
   * Append a size to the inner-most node of the index-tree
   *
   * This transforming of the given index-tree is used to implement
   * a blocked-interleaved index-merging strategy in a power-basis.
   *
   * Examples:
   * append( FlatIndexTree it, size ) -> UniformIndexTree( it.size(), FlatIndexTree(size) )
   * append( IndexTree(child...), size ) -> IndexTree( append(child, size)... )
   */
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

  // Append the size `s` at the inner-most node of the tree. Redirects to the
  // two-argument functions
  template<std::size_t s, class IT>
  auto appendToIndexTree (const IT& indexTree)
  {
    return appendToIndexTree(indexTree, index_constant<s>{});
  }

  } // end namespace Impl

}} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INDEXTREE_HH
