// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIZETREE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIZETREE_HH

#include <array>
#include <vector>

#include <dune/common/tuplevector.hh>
#include <dune/common/typeutilities.hh>


/**
 * \file sizetree.hh
 * \brief Lightweight representation of (hierarchic) size and blocking structure
 * The blocking tree can be used to define types for data-structures, like vectors
 * or matrices, that can be accessed by the multi-indices provided by a basis.
 **/

namespace Dune {
namespace Functions {

  //! Nested size-info that cannot be assigned one of the other size-tree types.
  struct UnknownSizeTree {};

  //! Leaf size-tree with static size and all sub-trees of size zero
  template <std::size_t n>
  struct StaticFlatSizeTree
  {
    template <class Index>
    StaticFlatSizeTree<0> operator[](const Index&) const { return {}; }

    static constexpr std::size_t size() { return n; }
  };

  //! Leaf size-tree with dynamic size and all sub-trees of size zero
  struct DynamicFlatSizeTree
  {
    template <class Index>
    StaticFlatSizeTree<0> operator[](const Index&) const { return {}; }

    std::size_t size() const { return size_; }

    std::size_t size_;
  };

  //! Non-uniform size-tree with all sub-trees of different type
  template<class... SubTrees>
  using NonUniformSizeTree = Dune::TupleVector<SubTrees...>;

  //! Non-uniform size-tree with all sub-tree of the same type and static size.
  template<class SubTree, std::size_t n>
  using StaticNonUniformSizeTree = std::array<SubTree, n>;

  //! Non-uniform size-tree with all sub-trees of the same type and dynamic size.
  template<class SubTree>
  using DynamicNonUniformSizeTree = std::vector<SubTree>;

  //! Uniform size-tree with static size.
  template<class SubTree, std::size_t n>
  struct UniformSizeTree
  {
    template <class Index>
    SubTree const& operator[](const Index&) const { return subTree_; }

    static constexpr std::size_t size() { return n; }

    SubTree subTree_;

    UniformSizeTree(std::integral_constant<std::size_t,n>, SubTree subTree)
      : subTree_{std::move(subTree)}
    {}
  };

  template<class SubTree, std::size_t n>
  using StaticUniformSizeTree = UniformSizeTree<SubTree,n>;

  //! Uniform size-tree with dynamic size.
  template<class SubTree>
  struct DynamicUniformSizeTree
  {
    template <class Index>
    SubTree const& operator[](const Index&) const { return subTree_; }

    std::size_t size() const { return size_; }

    std::size_t size_;
    SubTree subTree_;

    DynamicUniformSizeTree(std::size_t size, SubTree subTree)
      : size_{size}
      , subTree_{std::move(subTree)}
    {}
  };


  namespace Impl
  {
    template<class PreBasis>
    auto sizeTreeImpl(const PreBasis& preBasis, PriorityTag<2>)
      -> decltype(preBasis.sizeTree())
    {
      return preBasis.sizeTree();
    }

    template<class PreBasis>
    auto sizeTreeImpl(const PreBasis& preBasis, PriorityTag<1>)
      -> decltype(preBasis.dimension())
    {
      return DynamicFlatSizeTree{preBasis.dimension()};
    }
  }

  //! Generate a SizeTree associated to a PreBasis
  template<class PreBasis>
  auto sizeTree(const PreBasis& preBasis)
  {
    return Impl::sizeTreeImpl(preBasis, PriorityTag<5>{});
  }


  // -----------------------------------------------------------------------------------------------
  // Some utilities for to generate SizeTrees
  // -----------------------------------------------------------------------------------------------

  namespace Impl
  {
    template<class ST>
    struct SumSizeTrees
    {
      template<std::size_t s>
      static auto create(const ST& sizeTree)
      {
        return UnknownSizeTree{};
      }

      static auto create(const ST& sizeTree, std::size_t s)
      {
        return UnknownSizeTree{};
      }
    };

  } // end namespace Impl

  //! Generate a sum of a SizeTree consisting of n identical `sizeTree`s
  template<std::size_t n, class ST>
  auto sumSizeTrees(const ST& sizeTree)
  {
    return Impl::SumSizeTrees<ST>::template create<n>(sizeTree);
  };

  //! Generate a sum of a SizeTree consisting of n identical `sizeTree`s
  template<class ST>
  auto sumSizeTrees(const ST& sizeTree, std::size_t n)
  {
    return Impl::SumSizeTrees<ST>::create(sizeTree, n);
  };


  namespace Impl
  {
    template<>
    struct SumSizeTrees<DynamicFlatSizeTree>
    {
      using ST = DynamicFlatSizeTree;

      template<std::size_t s>
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

      template<std::size_t s>
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

      template<std::size_t s>
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

      template<std::size_t s>
      static auto create(const ST& sizeTree)
      {
        return DynamicUniformSizeTree<SubTree>{s*sizeTree.size_, sizeTree.subTree_};
      }

      static auto create(const ST& sizeTree, std::size_t s)
      {
        return DynamicUniformSizeTree<SubTree>{s*sizeTree.size_, sizeTree.subTree_};
      }
    };

  } // end namespace Impl



  namespace Impl
  {
    template<class ST1, class ST2>
    struct SumNonUniformSubTrees
    {
      static auto create(const ST1&, const ST2&)
      {
        return UnknownSizeTree{};
      }
    };

  } // end namespace Impl

  //! Overload for zero sizeTrees, return an unknown tree.
  inline auto sumNonUniformSubTrees()
  {
    UnknownSizeTree{};
  }

  //! Overload for one sizeTrees, return the tree itself.
  template<class ST>
  auto const& sumNonUniformSubTrees(const ST& sizeTree)
  {
    return sizeTree;
  }

  //! Merge size trees.
  template<class ST0, class... ST>
  auto sumNonUniformSubTrees(const ST0& sizeTree0, const ST&... sizeTrees)
  {
    return sumNonUniformSubTrees(sizeTree0, sumNonUniformSubTrees(sizeTrees...));
  }


  namespace Impl
  {
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

  } // end namespace Impl


  namespace Impl
  {
    template<class ST>
    struct AppendToSizeTree
    {
      static auto create(const ST&, std::size_t) { return UnknownSizeTree{}; }

      template<std::size_t>
      static auto create(const ST&) { return UnknownSizeTree{}; }
    };

  } // end namespace Impl

  //! Append the size `s` at the inner-most node of the tree
  template<class ST>
  auto appendToSizeTree(const ST& sizeTree, std::size_t s)
  {
    return Impl::AppendToSizeTree<ST>::create(sizeTree, s);
  }

  //! Append the size `s` at the inner-most node of the tree
  template<std::size_t s, class ST>
  auto appendToSizeTree(const ST& sizeTree)
  {
    return Impl::AppendToSizeTree<ST>::template create<s>(sizeTree);
  }


  namespace Impl
  {
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
          return NonUniformSizeTree<decltype(appendToSizeTree(subTree,s))...>
            {appendToSizeTree(subTree,s)...};
        }, sizeTree);
      }

      template<std::size_t s>
      static auto create(const ST& sizeTree)
      {
        return std::apply([&](auto const&... subTree) {
          return NonUniformSizeTree<decltype(appendToSizeTree<s>(subTree))...>
            {appendToSizeTree<s>(subTree)...};
        }, sizeTree);
      }
    };

    template<class SubTree, std::size_t n>
    struct AppendToSizeTree<UniformSizeTree<SubTree,n>>
    {
      using ST = UniformSizeTree<SubTree,n>;

      static auto create(const ST& sizeTree, std::size_t s)
      {
        auto subTree = appendToSizeTree(sizeTree.subTree_, s);
        return UniformSizeTree<decltype(subTree),n>{subTree};
      }

      template<std::size_t s>
      static auto create(const ST& sizeTree)
      {
        auto subTree = appendToSizeTree<s>(sizeTree.subTree_);
        return UniformSizeTree<decltype(subTree),n>{subTree};
      }
    };

    template<class SubTree>
    struct AppendToSizeTree<DynamicUniformSizeTree<SubTree>>
    {
      using ST = DynamicUniformSizeTree<SubTree>;

      static auto create(const ST& sizeTree, std::size_t s)
      {
        auto subTree = appendToSizeTree(sizeTree.subTree_, s);
        return DynamicUniformSizeTree<decltype(subTree)>{sizeTree.size_, subTree};
      }

      template<std::size_t s>
      static auto create(const ST& sizeTree)
      {
        auto subTree = appendToSizeTree<s>(sizeTree.subTree_);
        return DynamicUniformSizeTree<decltype(subTree)>{sizeTree.size_, subTree};
      }
    };

  } // end namespace Impl

}} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIZETREE_HH
