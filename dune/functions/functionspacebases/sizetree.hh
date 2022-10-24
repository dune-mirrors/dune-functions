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
 * The size tree can be used to define types for data-structures, like vectors
 * or matrices, that can be accessed by the multi-indices provided by a basis.
 **/

namespace Dune {
namespace Functions {

  //! Define basic properties of a tree
  template<bool leaf, bool uniform, bool typeUniform, bool staticSize>
  struct SizeTreePropertiesBase
  {
    //! The tree has no children
    inline static constexpr bool isLeaf = leaf;

    //! All nodes of the tree can be represented by one object
    inline static constexpr bool isUniform = uniform;

    //! All node types are identical
    inline static constexpr bool isTypeUniform = typeUniform;

    //! The tree has static degree
    inline static constexpr bool isStatic = staticSize;

    // for backwards compatibility
    inline static constexpr bool isPower = !isLeaf && isTypeUniform && !isUniform;
    inline static constexpr bool isComposite = !isLeaf && !isTypeUniform && !isUniform && isStatic;
  };

  //! The properties class has to be specialized for each SizeTree
  template<class SizeTree>
  struct SizeTreeProperties {};


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

  template <std::size_t n>
  struct SizeTreeProperties<StaticFlatSizeTree<n>>
      : public SizeTreePropertiesBase<true,true,true,true>
  {};


  //! Leaf size-tree with dynamic size and all sub-trees of size zero
  struct DynamicFlatSizeTree
  {
    template <class Index>
    StaticFlatSizeTree<0> operator[](const Index&) const { return {}; }

    std::size_t size() const { return size_; }

  private:
    std::size_t size_;
  };

  template <>
  struct SizeTreeProperties<DynamicFlatSizeTree>
      : public SizeTreePropertiesBase<true,true,true,false>
  {};


  //! Non-uniform size-tree with all sub-trees of different type
  template<class... SubTrees>
  using NonUniformSizeTree = Dune::TupleVector<SubTrees...>;

  template<class... SubTrees>
  struct SizeTreeProperties<NonUniformSizeTree<SubTrees...>>
      : public SizeTreePropertiesBase<false,false,false,true>
  {};


  //! Non-uniform size-tree with all sub-tree of the same type and static size.
  template<class SubTree, std::size_t n>
  using StaticNonUniformSizeTree = std::array<SubTree, n>;

  template<class SubTree, std::size_t n>
  struct SizeTreeProperties<StaticNonUniformSizeTree<SubTree,n>>
      : public SizeTreePropertiesBase<false,false,true,true>
  {};


  //! Non-uniform size-tree with all sub-trees of the same type and dynamic size.
  template<class SubTree>
  using DynamicNonUniformSizeTree = std::vector<SubTree>;

  template<class SubTree>
  struct SizeTreeProperties<DynamicNonUniformSizeTree<SubTree>>
      : public SizeTreePropertiesBase<false,false,true,false>
  {};


  //! Uniform size-tree with static size.
  template<class SubTree, std::size_t n>
  struct UniformSizeTree
  {
    template <class Index>
    SubTree const& operator[](const Index&) const { return subTree_; }

    static constexpr std::size_t size() { return n; }

    UniformSizeTree(SubTree subTree)
      : subTree_{std::move(subTree)}
    {}

    UniformSizeTree(std::integral_constant<std::size_t,n>, SubTree subTree)
      : subTree_{std::move(subTree)}
    {}

  private:
    SubTree subTree_;
  };

  template<class SubTree, std::size_t n>
  using StaticUniformSizeTree = UniformSizeTree<SubTree,n>;

  template<class SubTree, std::size_t n>
  struct SizeTreeProperties<StaticUniformSizeTree<SubTree,n>>
      : public SizeTreePropertiesBase<false,true,true,true>
  {};


  //! Uniform size-tree with dynamic size.
  template<class SubTree>
  struct DynamicUniformSizeTree
  {
    template <class Index>
    SubTree const& operator[](const Index&) const { return subTree_; }

    std::size_t size() const { return size_; }

    DynamicUniformSizeTree(std::size_t size, SubTree subTree)
      : size_{size}
      , subTree_{std::move(subTree)}
    {}

  private:
    std::size_t size_;
    SubTree subTree_;
  };

  template<class SubTree>
  struct SizeTreeProperties<DynamicUniformSizeTree<SubTree>>
      : public SizeTreePropertiesBase<false,true,true,false>
  {};


  namespace Impl
  {
    template<class PreBasis>
    auto sizeTreeImpl(const PreBasis& preBasis, PriorityTag<2>)
      -> decltype(preBasis.sizeTree())
    {
      return preBasis.sizeTree();
    }

    template<class PreBasis,
      class = decltype(std::declval<PreBasis>().dimension())>
    auto sizeTreeImpl(const PreBasis& preBasis, PriorityTag<1>)
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
        return DynamicFlatSizeTree{s*sizeTree.size()};
      }

      static auto create(const ST& sizeTree, std::size_t s)
      {
        return DynamicFlatSizeTree{s*sizeTree.size()};
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
        return StaticUniformSizeTree<SubTree,s*n>{sizeTree[0]};
      }

      static auto create(const ST& sizeTree, std::size_t s)
      {
        return DynamicUniformSizeTree<SubTree>{s*n, sizeTree[0]};
      }
    };

    template<class SubTree>
    struct SumSizeTrees<DynamicUniformSizeTree<SubTree>>
    {
      using ST = DynamicUniformSizeTree<SubTree>;

      template<std::size_t s>
      static auto create(const ST& sizeTree)
      {
        return DynamicUniformSizeTree<SubTree>{s*sizeTree.size(), sizeTree[0]};
      }

      static auto create(const ST& sizeTree, std::size_t s)
      {
        return DynamicUniformSizeTree<SubTree>{s*sizeTree.size(), sizeTree[0]};
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

  template<class ST0, class ST1>
  auto sumNonUniformSubTrees(const ST0& sizeTree0, const ST1& sizeTree1)
  {
    return Impl::SumNonUniformSubTrees<ST0,ST1>::create(sizeTree0, sizeTree1);
  }

  //! Merge size trees.
  template<class ST0, class... ST,
    std::enable_if_t<(sizeof...(ST) > 1), int> = 0>
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
        return DynamicFlatSizeTree{n + sizeTree2.size()};
      }
    };

    template<std::size_t n>
    struct SumNonUniformSubTrees<DynamicFlatSizeTree,StaticFlatSizeTree<n>>
    {
      using ST1 = DynamicFlatSizeTree;
      using ST2 = StaticFlatSizeTree<n>;

      static auto create(const ST1& sizeTree1, const ST2& sizeTree2)
      {
        return DynamicFlatSizeTree{sizeTree1.size() + n};
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
        return DynamicUniformSizeTree<SubTree>{sizeTree.size(), SubTree{s}};
      }

      template<std::size_t s>
      static auto create(const ST& sizeTree)
      {
        using SubTree = StaticFlatSizeTree<s>;
        return DynamicUniformSizeTree<SubTree>{sizeTree.size(), SubTree{}};
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
        auto subTree = appendToSizeTree(sizeTree[0], s);
        return UniformSizeTree<decltype(subTree),n>{std::move(subTree)};
      }

      template<std::size_t s>
      static auto create(const ST& sizeTree)
      {
        auto subTree = appendToSizeTree<s>(sizeTree[0]);
        return UniformSizeTree<decltype(subTree),n>{std::move(subTree)};
      }
    };

    template<class SubTree>
    struct AppendToSizeTree<DynamicUniformSizeTree<SubTree>>
    {
      using ST = DynamicUniformSizeTree<SubTree>;

      static auto create(const ST& sizeTree, std::size_t s)
      {
        auto subTree = appendToSizeTree(sizeTree[0], s);
        return DynamicUniformSizeTree<decltype(subTree)>{sizeTree.size(), std::move(subTree)};
      }

      template<std::size_t s>
      static auto create(const ST& sizeTree)
      {
        auto subTree = appendToSizeTree<s>(sizeTree[0]);
        return DynamicUniformSizeTree<decltype(subTree)>{sizeTree.size(), std::move(subTree)};
      }
    };

  } // end namespace Impl

}} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIZETREE_HH
