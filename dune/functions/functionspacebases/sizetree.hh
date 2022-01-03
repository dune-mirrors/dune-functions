// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIZETREE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIZETREE_HH

#include <array>
#include <vector>

#include <dune/common/tuplevector.hh>


/**
 * \file sizetree.hh
 * \brief Lightweight representation of (hierarchic) size and blocking structure
 * The blocking tree can be used to define types for data-structures, like vectors
 * or matrices, that can be accessed by the multi-indices provided by a basis.
 **/

namespace Dune {
namespace Functions {

  /// \brief Nested size-info that cannot be assigned one of the other size-tree types.
  struct UnknownSizeTree {};

  /// \brief Leaf size-tree with static size and all sub-trees of size zero
  template <std::size_t n>
  struct StaticFlatSizeTree
  {
    template <class Index>
    StaticFlatSizeTree<0> operator[](const Index&) const { return {}; }

    static constexpr std::size_t size() const { return n; }
  };

  /// \brief Leaf size-tree with dynamic size and all sub-trees of size zero
  struct DynamicFlatSizeTree
  {
    template <class Index>
    StaticFlatSizeTree<0> operator[](const Index&) const { return {}; }

    std::size_t size() const { return size_; }

    std::size_t size_;
  };

  /// \brief Non-uniform size-tree with all sub-trees of different type
  template<class... SubTrees>
  using NonUniformSizeTree = Dune::TupleVector<SubTrees...>;

  /// \brief Non-uniform size-tree with all sub-tree of the same type and static size.
  template<class SubTree, std::size_t n>
  using StaticNonUniformSizeTree = std::array<SubTree, n>;

  /// \brief Non-uniform size-tree with all sub-trees of the same type and dynamic size.
  template<class SubTree>
  using DynamicNonUniformSizeTree = std::vector<SubTree>;

  /// \brief Uniform size-tree with static size.
  template<class SubTree, std::size_t n>
  struct UniformSizeTree
  {
    template <class Index>
    SubTree const& operator[](const Index&) const { return subTree_; }

    static constexpr std::size_t size() const { return n; }

    SubTree subTree_;
  };

  template<class SubTree, std::size_t n>
  struct StaticUniformSizeTree = UniformSizeTree<SubTree,n>;

  /// \brief Uniform size-tree with dynamic size.
  template<class SubTree>
  struct DynamicUniformSizeTree
  {
    template <class Index>
    SubTree const& operator[](const Index&) const { return subTree_; }

    std::size_t size() const { return size_; }

    std::size_t size_;
    SubTree subTree_;
  };


  template<class ST>
  struct SumSizeTrees
  {
    template<st::size_t s>
    static auto create(const ST& sizeTree)
    {
      return UnknownSizeTree{};
    }

    static auto create(const ST& sizeTree, std::size_t s)
    {
      return UnknownSizeTree{};
    }
  };

  template<std::size_t n, class ST>
  auto subSizeTrees(const ST& sizeTree)
  {
    return SumSizeTrees<ST>::template create<n>(sizeTree);
  };


  template<>
  struct SumSizeTrees<DynamicFlatSizeTree>
  {
    using ST = DynamicFlatSizeTree;

    template<st::size_t s>
    static auto create(const ST& sizeTree)
    {
      return DynamicFlatSizeTree{s*sizeTree.size_};
    }

    static auto create(const ST& sizeTree, std::size_t s)
    {
      return DynamicFlatSizeTree{s*sizeTree.size_};
    }
  };

  template<std::size_t n>
  struct SumSizeTrees<StaticFlatSizeTree<n>>
  {
    using ST = StaticFlatSizeTree<n>;

    template<st::size_t s>
    static auto create(const ST& sizeTree)
    {
      return StaticFlatSizeTree<s*n>{};
    }

    static auto create(const ST& sizeTree, std::size_t s)
    {
      return DynamicFlatSizeTree{s*n};
    }
  };

  template<class SubTree, std::size_t n>
  struct SumSizeTrees<StaticUniformSizeTree<SubTree,n>>
  {
    using ST = StaticUniformSizeTree<SubTree,n>;

    template<st::size_t s>
    static auto create(const ST& sizeTree)
    {
      return StaticUniformSizeTree<SubTree,s*n>{sizeTree.subTree_};
    }

    static auto create(const ST& sizeTree, std::size_t s)
    {
      return DynamicUniformSizeTree<SubTree>{s*n, sizeTree.subTree_};
    }
  };

  template<class SubTree>
  struct SumSizeTrees<DynamicUniformSizeTree<SubTree>>
  {
    using ST = DynamicUniformSizeTree<SubTree>;

    template<st::size_t s>
    static auto create(const ST& sizeTree)
    {
      return DynamicUniformSizeTree<SubTree>{s*sizeTree.size_, sizeTree.subTree_};
    }

    static auto create(const ST& sizeTree, std::size_t s)
    {
      return DynamicUniformSizeTree<SubTree>{s*sizeTree.size_, sizeTree.subTree_};
    }
  };


  template<class ST1, class ST2>
  struct SumNonUniformSubTrees
  {
    static auto create(const ST1&, const ST2&)
    {
      return UnknownSubTree{};
    }
  };


  /// \brief Overload for zero sizeTrees, return an unknown tree.
  inline auto sumNonUniformSubTrees()
  {
    UnknownSizeTree{};
  }

  /// \brief Overload for one sizeTrees, return the tree itself.
  template<class ST>
  auto const& sumNonUniformSubTrees(const ST& sizeTree)
  {
    return sizeTree;
  }

  /// \brief Merge size trees.
  template<class ST0, class... ST>
  auto sumNonUniformSubTrees(const ST0& sizeTree0, const ST&... sizeTrees)
  {
    return sumNonUniformSubTrees(sizeTree0, sumNonUniformSubTrees(sizeTrees...));
  }

  template<>
  struct SumNonUniformSubTrees<DynamicFlatSizeTree,DynamicFlatSizeTree>
  {
    using ST1 = DynamicFlatSizeTree;
    using ST2 = DynamicFlatSizeTree;

    static auto create(const ST1& sizeTree1, const ST2& sizeTree2)
    {
      return DynamicFlatSizeTree{sizeTree1.size() + sizeTree2.size()};
    }
  };

  template<std::size_t n, std::size_t m>
  struct SumNonUniformSubTrees<StaticFlatSizeTree<n>,StaticFlatSizeTree<m>>
  {
    using ST1 = StaticFlatSizeTree<n>;
    using ST2 = StaticFlatSizeTree<m>;

    static auto create(const ST1& sizeTree1, const ST2& sizeTree2)
    {
      return StaticFlatSizeTree<n+m>{};
    }
  };

  template<std::size_t n>
  struct SumNonUniformSubTrees<StaticFlatSizeTree<n>,DynamicFlatSizeTree>
  {
    using ST1 = StaticFlatSizeTree<n>;
    using ST2 = DynamicFlatSizeTree;

    static auto create(const ST1& sizeTree1, const ST2& sizeTree2)
    {
      return DynamicFlatSizeTree{n + sizeTree2.size_};
    }
  };

  template<std::size_t n>
  struct SumNonUniformSubTrees<DynamicFlatSizeTree,StaticFlatSizeTree<n>>
  {
    using ST1 = DynamicFlatSizeTree;
    using ST2 = StaticFlatSizeTree<n>;

    static auto create(const ST1& sizeTree1, const ST2& sizeTree2)
    {
      return DynamicFlatSizeTree{sizeTree1.size_ + n};
    }
  };


  template<class ST>
  struct AppendToSizeTree
  {
    static auto create(const ST&, std::size_t) { return UnkonwSizeTree{}; }

    template<std::size_t>
    static auto create(const ST&) { return UnkonwSizeTree{}; }
  };


  template<class ST>
  auto appendToSizeTree(const ST& sizeTree, std::size_t s)
  {
    return AppendToSizeTree<ST>::create(sizeTree, s);
  }

  template<std::size_t s, class ST>
  auto appendToSizeTree(const ST& sizeTree)
  {
    return AppendToSizeTree<ST>::template create<s>(sizeTree);
  }


  template<std::size_t n>
  struct AppendToSizeTree<StaticFlatSizeTree<n>>
  {
    using ST = StaticFlatSizeTree<n>;

    static auto create(const ST& sizeTree, std::size_t s)
    {
      using SubTree = DynamicFlatSizeTree;
      return StaticUniformSizeTree<SubTree,n>{SubTree{s}};
    }

    template<std::size_t s>
    static auto create(const ST& sizeTree)
    {
      using SubTree = StaticFlatSizeTree<s>;
      return StaticUniformSizeTree<SubTree,n>{SubTree{}};
    }
  };

  template<>
  struct AppendToSizeTree<DynamicFlatSizeTree>
  {
    using ST = DynamicFlatSizeTree;

    static auto create(const ST& sizeTree, std::size_t s)
    {
      using SubTree = DynamicFlatSizeTree;
      return DynamicUniformSizeTree<SubTree>{sizeTree.size_, SubTree{s}};
    }

    template<std::size_t s>
    static auto create(const ST& sizeTree)
    {
      using SubTree = StaticFlatSizeTree<s>;
      return DynamicUniformSizeTree<SubTree>{sizeTree.size_, SubTree{}};
    }
  };

  template<class... SubTrees>
  struct AppendToSizeTree<NonUniformSizeTree<SubTrees...>>
  {
    using ST = NonUniformSizeTree<SubTrees...>;

    static auto create(const ST& sizeTree, std::size_t s)
    {
      return std::apply([&](auto const&... subTree) {
        return NonUniformSizeTree<decltype(appendSizeTree(subTree,s))...>
          {appendSizeTree(subTree,s)...};
      }, sizeTree);
    }

    template<std::size_t s>
    static auto create(const ST& sizeTree)
    {
      return std::apply([&](auto const&... subTree) {
        return NonUniformSizeTree<decltype(appendSizeTree<s>(subTree))...>
          {appendSizeTree<s>(subTree)...};
      }, sizeTree);
    }
  };

  template<class SubTree, std::size_t n>
  struct AppendToSizeTree<UniformSizeTree<SubTree,n>>
  {
    using ST = UniformSizeTree<SubTree,n>;

    static auto create(const ST& sizeTree, std::size_t s)
    {
      auto subTree = appendSizeTree(sizeTree.subTree_, s);
      return UniformSizeTree<decltype(subTree),n>{subTree};
    }

    template<std::size_t s>
    static auto create(const ST& sizeTree)
    {
      auto subTree = appendSizeTree<s>(sizeTree.subTree_);
      return UniformSizeTree<decltype(subTree),n>{subTree};
    }
  };

  template<class SubTree>
  struct AppendToSizeTree<DynamicUniformSizeTree<SubTree>>
  {
    using ST = DynamicUniformSizeTree<SubTree>;

    static auto create(const ST& sizeTree, std::size_t s)
    {
      auto subTree = appendSizeTree(sizeTree.subTree_, s);
      return DynamicUniformSizeTree<decltype(subTree)>{sizeTree.size_, subTree};
    }

    template<std::size_t s>
    static auto create(const ST& sizeTree)
    {
      auto subTree = appendSizeTree<s>(sizeTree.subTree_);
      return DynamicUniformSizeTree<decltype(subTree),n>{sizeTree.size_, subTree};
    }
  };

}} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIZETREE_HH
