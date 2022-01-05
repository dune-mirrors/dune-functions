// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_TYPETREE2_TYPETREE_HH
#define DUNE_TYPETREE2_TYPETREE_HH

#include <array>
#include <vector>

#include <dune/common/indices.hh>
#include <dune/common/tuplevector.hh>
#include <dune/common/typelist.hh>
#include <dune/common/typeutilities.hh>

namespace Dune {
namespace TypeTree2 {

  template <bool b> struct IsLeaf : std::bool_constant<b> {};
  template <bool b> struct IsUniform : std::bool_constant<b> {};
  template <bool b> struct IsTypeUniform : std::bool_constant<b> {};
  template <bool b> struct IsStatic : std::bool_constant<b> {};

  template<class, class, class, class>
  struct BaseTypeTree;

  template<bool leaf, bool uniform, bool typeUniform, bool staticSize>
  struct BaseTypeTree<IsLeaf<leaf>, IsUniform<uniform>,
                      IsTypeUniform<typeUniform>, IsStatic<staticSize>>
  {
    inline static constexpr bool isLeaf = leaf;
    inline static constexpr bool isUniform = uniform;
    inline static constexpr bool isTypeUniform = typeUniform;
    inline static constexpr bool isStatic = staticSize;
  };



  //! Nested size-info that cannot be assigned one of the other type-tree types.
  struct UnknownTypeTree
      : public BaseTypeTree<IsLeaf<false>,IsUniform<false>,IsTypeUniform<false>,IsStatic<false>>
  {};

  //! Leaf type-tree with degree 0
  struct LeafTypeTree
      : public BaseTypeTree<IsLeaf<true>,IsUniform<false>,IsTypeUniform<false>,IsStatic<true>>
  {
    static constexpr index_constant<0> degree() { return {}; }
  };


  //! Non-uniform type-tree with all sub-trees of different type
  template<class... SubTrees>
  struct NonUniformTypeTree
      : public BaseTypeTree<IsLeaf<false>,IsUniform<false>,IsTypeUniform<false>,IsStatic<true>>
      , public std::tuple<SubTrees...>
  {
    using Super = std::tuple<SubTrees...>;
    using Super::Super;

    using ChildTypes = Super;

    template <std::size_t i>
    using Child = std::tuple_element_t<i,Super>;

    template<std::size_t i>
    auto& child(index_constant<i> ii)
    {
      return std::get<i>(static_cast<Super&>(*this));
    }

    template<std::size_t i>
    auto& child(index_constant<i> ii) const
    {
      return std::get<i>(static_cast<Super const&>(*this));
    }

    static constexpr index_constant<(sizeof...(SubTrees))> degree() { return {}; }
  };

  //! Non-uniform type-tree with all sub-tree of the same type and static size.
  template<class SubTree, std::size_t n>
  struct StaticNonUniformTypeTree
      : public BaseTypeTree<IsLeaf<false>,IsUniform<false>,IsTypeUniform<true>,IsStatic<true>>
      , public std::array<SubTree,n>
  {
    using Child = SubTree;
    using Super = std::array<SubTree, n>;

    //! Explicit definition of default constructor
    StaticNonUniformTypeTree()
      : Super{}
    {}

    //! Forward all arguments to the array
    template<class... Args,
      Dune::disableCopyMove<StaticNonUniformTypeTree, Args...> = 0,
      std::enable_if_t<(sizeof...(Args) == n || n == 1), int> = 0>
    explicit StaticNonUniformTypeTree(Args&&... args)
      : Super{std::forward<Args>(args)...}
    {}

    //! Repeat the single argument to fill the array
    explicit StaticNonUniformTypeTree(SubTree const& subTree)
      : Super{Dune::unpackIntegerSequence(
          [&](auto... ii) {
            return Super{(void(ii), subTree)...};
          }, std::make_index_sequence<n>{})
        }
    {}

    auto& child(std::size_t i)       { return Super::operator[](i); }
    auto& child(std::size_t i) const { return Super::operator[](i); }

    static constexpr index_constant<n> degree() { return {}; }
  };

  //! Non-uniform type-tree with all sub-trees of the same type and dynamic size.
  template<class SubTree>
  struct DynamicNonUniformTypeTree
      : public BaseTypeTree<IsLeaf<false>,IsUniform<false>,IsTypeUniform<true>,IsStatic<false>>
      , public std::vector<SubTree>
  {
    using Child = SubTree;
    using Super = std::vector<SubTree>;
    using Super::Super;

    auto& child(std::size_t i)       { return Super::operator[](i); }
    auto& child(std::size_t i) const { return Super::operator[](i); }

    std::size_t degree() { return Super::size(); }
  };

  //! Uniform type-tree with static size.
  template<class SubTree, std::size_t n>
  struct StaticUniformTypeTree
      : public BaseTypeTree<IsLeaf<false>,IsUniform<true>,IsTypeUniform<true>,IsStatic<true>>
  {
    using Child = SubTree;

    StaticUniformTypeTree(std::integral_constant<std::size_t,n>, SubTree subTree)
      : subTree_{std::move(subTree)}
    {}

    template <class Index>
    SubTree const& child(const Index& i) const
    {
      assert(std::size_t(i) < n);
      return subTree_;
    }

    template <class Index>
    SubTree& child(const Index& i)
    {
      assert(std::size_t(i) < n);
      return subTree_;
    }

    static constexpr index_constant<n> degree() { return {}; }

    SubTree subTree_;
  };

  //! Uniform type-tree with dynamic size.
  template<class SubTree>
  struct DynamicUniformTypeTree
      : public BaseTypeTree<IsLeaf<false>,IsUniform<true>,IsTypeUniform<true>,IsStatic<false>>
  {
    using Child = SubTree;

    DynamicUniformTypeTree(std::size_t size, SubTree subTree)
      : size_{size}
      , subTree_{std::move(subTree)}
    {}

    template <class Index>
    SubTree const& child(const Index& i) const
    {
      assert(std::size_t(i) < size_);
      return subTree_;
    }

    template <class Index>
    SubTree& child(const Index& i)
    {
      assert(std::size_t(i) < size_);
      return subTree_;
    }

    std::size_t degree() const { return size_; }

    std::size_t size_;
    SubTree subTree_;
  };

}} // end namespace Dune::TypeTree2

#endif // DUNE_TYPETREE2_TYPETREETREE_HH
