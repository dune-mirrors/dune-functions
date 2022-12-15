// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIZETREE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIZETREE_HH

#include <any>
#include <array>
#include <vector>

#include <dune/common/filledarray.hh>
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

  // Some size-tree properties
  struct nonUniform_t {};
  struct typeUniform_t {};
  struct uniform_t : typeUniform_t {};
  struct flat_t : uniform_t {};

  template <class ST>
  inline constexpr bool isNonUniform = std::is_base_of_v<nonUniform_t, typename ST::State>;

  template <class ST>
  inline constexpr bool isTypeUniform = std::is_base_of_v<typeUniform_t, typename ST::State>;

  template <class ST>
  inline constexpr bool isUniform = std::is_base_of_v<uniform_t, typename ST::State>;

  template <class ST>
  inline constexpr bool isFlat = std::is_base_of_v<flat_t, typename ST::State>;

  template <class ST>
  inline constexpr bool isStatic = ST::isStatic;

  struct UnknownSizeTree {};
  struct EmptySizeTree
  {
    template <class Index>
    EmptySizeTree operator[] (const Index&) const { return {}; }

    static constexpr std::size_t size () { return 0; }
  };

  //! Leaf size-tree with static size and all sub-trees of size zero
  template <std::size_t n>
  struct StaticFlatSizeTree
  {
    using State = flat_t;
    static constexpr bool isStatic = true;

    template <class Index>
    EmptySizeTree operator[] (const Index&) const { return {}; }

    static constexpr std::size_t size () { return n; }
  };


  //! Leaf size-tree with dynamic size and all sub-trees of size zero
  struct FlatSizeTree
  {
    using State = flat_t;
    static constexpr bool isStatic = false;

    explicit FlatSizeTree (std::size_t size)
      : size_(size)
    {}

    template <class Index>
    EmptySizeTree operator[] (const Index&) const { return {}; }

    std::size_t size () const { return size_; }

  private:
    std::size_t size_;
  };


  //! Non-uniform size-tree with all sub-trees of different type
  template<class... SubTrees>
  struct StaticNonUniformSizeTree
      : Dune::TupleVector<SubTrees...>
  {
    using Super = Dune::TupleVector<SubTrees...>;
    using State = nonUniform_t;
    static constexpr bool isStatic = true;

    using Super::Super;
  };


  //! Non-uniform size-tree with all sub-trees of different type but dynamic size
  //! NOTE, this cannot easily be implemented, maybe using type-erasure
  struct NonUniformSizeTree
      : std::vector<std::any>
  {
    using Super = std::vector<std::any>;
    using State = nonUniform_t;
    static constexpr bool isStatic = false;

    template <class... SubTrees>
    NonUniformSizeTree (const SubTrees&... subTrees)
      : Super{std::any(subTrees)...}
    {}
  };


  //! Non-uniform size-tree with all sub-tree of the same type and static size.
  template<class SubTree, std::size_t n>
  struct StaticTypeUniformSizeTree
      : std::array<SubTree, n>
  {
    using Super = std::array<SubTree, n>;
    using State = typeUniform_t;
    static constexpr bool isStatic = true;

    StaticTypeUniformSizeTree (SubTree subTree)
      : Super{Dune::filledArray<n>(std::move(subTree))}
    {}

    StaticTypeUniformSizeTree (std::integral_constant<std::size_t,n>, SubTree subTree)
      : Super{Dune::filledArray<n>(std::move(subTree))}
    {}

    template <class... SubTrees,
      std::enable_if_t<(std::is_same_v<SubTrees,SubTree> &&...), int> = 0>
    StaticTypeUniformSizeTree (SubTrees... subTrees)
      : Super{std::move(subTrees)...}
    {}
  };


  //! Non-uniform size-tree with all sub-trees of the same type and dynamic size.
  template<class SubTree>
  struct TypeUniformSizeTree
      : std::vector<SubTree>
  {
    using Super = std::vector<SubTree>;
    using State = typeUniform_t;
    static constexpr bool isStatic = false;

    using Super::Super;
  };


  //! Uniform size-tree with static size.
  template<class SubTree, std::size_t n>
  struct StaticUniformSizeTree
  {
    using State = uniform_t;
    static constexpr bool isStatic = true;

    StaticUniformSizeTree (SubTree subTree)
      : subTree_{std::move(subTree)}
    {}

    StaticUniformSizeTree (std::integral_constant<std::size_t,n>, SubTree subTree)
      : subTree_{std::move(subTree)}
    {}

    template <class Index>
    SubTree const& operator[] (const Index&) const { return subTree_; }

    static constexpr std::size_t size () { return n; }

  private:
    SubTree subTree_;
  };


  //! Uniform size-tree with dynamic size.
  template<class SubTree>
  struct UniformSizeTree
  {
    using State = uniform_t;
    static constexpr bool isStatic = false;

    UniformSizeTree (std::size_t size, SubTree subTree)
      : size_{size}
      , subTree_{std::move(subTree)}
    {}

    template <class Index>
    SubTree const& operator[] (const Index&) const { return subTree_; }

    std::size_t size () const { return size_; }

  private:
    std::size_t size_;
    SubTree subTree_;
  };


  namespace Impl
  {
    template<class PreBasis>
    auto sizeTreeImpl (const PreBasis& preBasis, PriorityTag<2>)
      -> decltype(preBasis.sizeTree())
    {
      return preBasis.sizeTree();
    }

    template<class PreBasis,
      class = decltype(std::declval<PreBasis>().dimension())>
    auto sizeTreeImpl (const PreBasis& preBasis, PriorityTag<1>)
    {
      return FlatSizeTree{preBasis.dimension()};
    }
  }

  //! Generate a SizeTree associated to a PreBasis
  template<class PreBasis>
  auto sizeTree (const PreBasis& preBasis)
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
      static auto create (const ST& sizeTree)
      {
        return UnknownSizeTree{};
      }

      static auto create (const ST& sizeTree, std::size_t s)
      {
        return UnknownSizeTree{};
      }
    };

  } // end namespace Impl

  //! Generate a sum of a SizeTree consisting of n identical `sizeTree`s
  template<std::size_t n, class ST>
  auto sumSizeTrees (const ST& sizeTree)
  {
    return Impl::SumSizeTrees<ST>::template create<n>(sizeTree);
  }

  //! Generate a sum of a SizeTree consisting of n identical `sizeTree`s
  template<class ST>
  auto sumSizeTrees (const ST& sizeTree, std::size_t n)
  {
    return Impl::SumSizeTrees<ST>::create(sizeTree, n);
  }


  namespace Impl
  {
    template<>
    struct SumSizeTrees<FlatSizeTree>
    {
      using ST = FlatSizeTree;

      template<std::size_t s>
      static auto create (const ST& sizeTree)
      {
        return FlatSizeTree{s*sizeTree.size()};
      }

      static auto create (const ST& sizeTree, std::size_t s)
      {
        return FlatSizeTree{s*sizeTree.size()};
      }
    };

    template<std::size_t n>
    struct SumSizeTrees<StaticFlatSizeTree<n>>
    {
      using ST = StaticFlatSizeTree<n>;

      template<std::size_t s>
      static auto create (const ST& sizeTree)
      {
        return StaticFlatSizeTree<s*n>{};
      }

      static auto create (const ST& sizeTree, std::size_t s)
      {
        return FlatSizeTree{s*n};
      }
    };

    template<class SubTree, std::size_t n>
    struct SumSizeTrees<StaticUniformSizeTree<SubTree,n>>
    {
      using ST = StaticUniformSizeTree<SubTree,n>;

      template<std::size_t s>
      static auto create (const ST& sizeTree)
      {
        return StaticUniformSizeTree<SubTree,s*n>{sizeTree[0]};
      }

      static auto create (const ST& sizeTree, std::size_t s)
      {
        return UniformSizeTree<SubTree>{s*n, sizeTree[0]};
      }
    };

    template<class SubTree>
    struct SumSizeTrees<UniformSizeTree<SubTree>>
    {
      using ST = UniformSizeTree<SubTree>;

      template<std::size_t s>
      static auto create (const ST& sizeTree)
      {
        return UniformSizeTree<SubTree>{s*sizeTree.size(), sizeTree[0]};
      }

      static auto create (const ST& sizeTree, std::size_t s)
      {
        return UniformSizeTree<SubTree>{s*sizeTree.size(), sizeTree[0]};
      }
    };

  } // end namespace Impl



  namespace Impl
  {
    template<class ST1, class ST2>
    struct SumNonUniformSubTrees
    {
      static auto create (const ST1&, const ST2&)
      {
        return UnknownSizeTree{};
      }
    };

  } // end namespace Impl

  //! Overload for zero sizeTrees, return an unknown tree.
  inline auto sumNonUniformSubTrees ()
  {
    UnknownSizeTree{};
  }

  //! Overload for one sizeTrees, return the tree itself.
  template<class ST>
  auto const& sumNonUniformSubTrees (const ST& sizeTree)
  {
    return sizeTree;
  }

  template<class ST0, class ST1>
  auto sumNonUniformSubTrees (const ST0& sizeTree0, const ST1& sizeTree1)
  {
    return Impl::SumNonUniformSubTrees<ST0,ST1>::create(sizeTree0, sizeTree1);
  }

  //! Merge size trees.
  template<class ST0, class... ST,
    std::enable_if_t<(sizeof...(ST) > 1), int> = 0>
  auto sumNonUniformSubTrees (const ST0& sizeTree0, const ST&... sizeTrees)
  {
    return sumNonUniformSubTrees(sizeTree0, sumNonUniformSubTrees(sizeTrees...));
  }


  namespace Impl
  {
    template<>
    struct SumNonUniformSubTrees<FlatSizeTree,FlatSizeTree>
    {
      using ST1 = FlatSizeTree;
      using ST2 = FlatSizeTree;

      static auto create (const ST1& sizeTree1, const ST2& sizeTree2)
      {
        return FlatSizeTree{sizeTree1.size() + sizeTree2.size()};
      }
    };

    template<std::size_t n, std::size_t m>
    struct SumNonUniformSubTrees<StaticFlatSizeTree<n>,StaticFlatSizeTree<m>>
    {
      using ST1 = StaticFlatSizeTree<n>;
      using ST2 = StaticFlatSizeTree<m>;

      static auto create (const ST1& sizeTree1, const ST2& sizeTree2)
      {
        return StaticFlatSizeTree<n+m>{};
      }
    };

    template<std::size_t n>
    struct SumNonUniformSubTrees<StaticFlatSizeTree<n>,FlatSizeTree>
    {
      using ST1 = StaticFlatSizeTree<n>;
      using ST2 = FlatSizeTree;

      static auto create (const ST1& sizeTree1, const ST2& sizeTree2)
      {
        return FlatSizeTree{n + sizeTree2.size()};
      }
    };

    template<std::size_t n>
    struct SumNonUniformSubTrees<FlatSizeTree,StaticFlatSizeTree<n>>
    {
      using ST1 = FlatSizeTree;
      using ST2 = StaticFlatSizeTree<n>;

      static auto create (const ST1& sizeTree1, const ST2& sizeTree2)
      {
        return FlatSizeTree{sizeTree1.size() + n};
      }
    };

  } // end namespace Impl


  namespace Impl
  {
    template<class ST>
    struct AppendToSizeTree
    {
      static auto create (const ST&, std::size_t) { return UnknownSizeTree{}; }

      template<std::size_t>
      static auto create (const ST&) { return UnknownSizeTree{}; }
    };

  } // end namespace Impl

  //! Append the size `s` at the inner-most node of the tree
  template<class ST>
  auto appendToSizeTree (const ST& sizeTree, std::size_t s)
  {
    return Impl::AppendToSizeTree<ST>::create(sizeTree, s);
  }

  //! Append the size `s` at the inner-most node of the tree
  template<std::size_t s, class ST>
  auto appendToSizeTree (const ST& sizeTree)
  {
    return Impl::AppendToSizeTree<ST>::template create<s>(sizeTree);
  }


  namespace Impl
  {
    template<std::size_t n>
    struct AppendToSizeTree<StaticFlatSizeTree<n>>
    {
      using ST = StaticFlatSizeTree<n>;

      static auto create (const ST& sizeTree, std::size_t s)
      {
        using SubTree = FlatSizeTree;
        return StaticUniformSizeTree<SubTree,n>{SubTree{s}};
      }

      template<std::size_t s>
      static auto create (const ST& sizeTree)
      {
        using SubTree = StaticFlatSizeTree<s>;
        return StaticUniformSizeTree<SubTree,n>{SubTree{}};
      }
    };

    template<>
    struct AppendToSizeTree<FlatSizeTree>
    {
      using ST = FlatSizeTree;

      static auto create (const ST& sizeTree, std::size_t s)
      {
        using SubTree = FlatSizeTree;
        return UniformSizeTree<SubTree>{sizeTree.size(), SubTree{s}};
      }

      template<std::size_t s>
      static auto create (const ST& sizeTree)
      {
        using SubTree = StaticFlatSizeTree<s>;
        return UniformSizeTree<SubTree>{sizeTree.size(), SubTree{}};
      }
    };

    template<class... SubTrees>
    struct AppendToSizeTree<StaticNonUniformSizeTree<SubTrees...>>
    {
      using ST = StaticNonUniformSizeTree<SubTrees...>;

      static auto create (const ST& sizeTree, std::size_t s)
      {
        return std::apply([&](auto const&... subTree) {
          return StaticNonUniformSizeTree<decltype(appendToSizeTree(subTree,s))...>
            {appendToSizeTree(subTree,s)...};
        }, sizeTree);
      }

      template<std::size_t s>
      static auto create (const ST& sizeTree)
      {
        return std::apply([&](auto const&... subTree) {
          return StaticNonUniformSizeTree<decltype(appendToSizeTree<s>(subTree))...>
            {appendToSizeTree<s>(subTree)...};
        }, sizeTree);
      }
    };

    template<class SubTree, std::size_t n>
    struct AppendToSizeTree<StaticUniformSizeTree<SubTree,n>>
    {
      using ST = StaticUniformSizeTree<SubTree,n>;

      static auto create (const ST& sizeTree, std::size_t s)
      {
        auto subTree = appendToSizeTree(sizeTree[0], s);
        return StaticUniformSizeTree<decltype(subTree),n>{std::move(subTree)};
      }

      template<std::size_t s>
      static auto create (const ST& sizeTree)
      {
        auto subTree = appendToSizeTree<s>(sizeTree[0]);
        return StaticUniformSizeTree<decltype(subTree),n>{std::move(subTree)};
      }
    };

    template<class SubTree>
    struct AppendToSizeTree<UniformSizeTree<SubTree>>
    {
      using ST = UniformSizeTree<SubTree>;

      static auto create (const ST& sizeTree, std::size_t s)
      {
        auto subTree = appendToSizeTree(sizeTree[0], s);
        return UniformSizeTree<decltype(subTree)>{sizeTree.size(), std::move(subTree)};
      }

      template<std::size_t s>
      static auto create (const ST& sizeTree)
      {
        auto subTree = appendToSizeTree<s>(sizeTree[0]);
        return UniformSizeTree<decltype(subTree)>{sizeTree.size(), std::move(subTree)};
      }
    };

  } // end namespace Impl

}} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIZETREE_HH
