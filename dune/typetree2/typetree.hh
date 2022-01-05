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
  struct TreeProperties;

  //! Define basic properties of a tree
  template<bool leaf, bool uniform, bool typeUniform, bool staticSize>
  struct TreeProperties<
    IsLeaf<leaf>, IsUniform<uniform>, IsTypeUniform<typeUniform>, IsStatic<staticSize> >
  {
    //! The tree has no children
    inline static constexpr bool isLeaf = leaf;

    //! All nodes of the tree can be represented by one object
    inline static constexpr bool isUniform = uniform;

    //! All node type are identical
    inline static constexpr bool isTypeUniform = typeUniform;

    //! The tree has static degree
    inline static constexpr bool isStatic = staticSize;

    // for backwards compatibility
    inline static constexpr bool isPower = !isLeaf && isTypeUniform && !isUniform;
    inline static constexpr bool isComposite = !isLeaf && !isTypeUniform && !isUniform && isStatic;
  };


  //! A type-tree that cannot be represented with the other tree types
  struct UnknownTypeTree
  {};


  //! Leaf tree node with degree 0
  struct LeafNode
      : public TreeProperties<IsLeaf<true>,IsUniform<false>,IsTypeUniform<false>,IsStatic<true>>
  {
    static constexpr index_constant<0> degree() { return {}; }
  };


  //! Non-uniform type-tree with all sub-trees of different type
  template<class... SubTrees>
  struct VariadicNonUniformTypeTree
      : public TreeProperties<IsLeaf<false>,IsUniform<false>,IsTypeUniform<false>,IsStatic<true>>
      , public std::tuple<SubTrees...>
  {
    using Super = std::tuple<SubTrees...>;

    //! The type of the childs tuple
    using ChildTypes = Super;

    //! The type of the i'th child
    template <std::size_t i>
    using Child = std::tuple_element_t<i,Super>;

    //! Inherit the constructors from std::tuple
    using Super::Super;

    //! Return a reference to the i'th child of the tree
    template<std::size_t i>
    auto& child(index_constant<i> ii)
    {
      return std::get<i>(static_cast<Super&>(*this));
    }

    //! Return a const reference to the i'th child of the tree
    template<std::size_t i>
    auto& child(index_constant<i> ii) const
    {
      return std::get<i>(static_cast<Super const&>(*this));
    }

    //! Return the number of nodes
    static constexpr index_constant<(sizeof...(SubTrees))> degree() { return {}; }
  };


  //! Non-uniform type-tree with all sub-tree of the same type and static size.
  template<class SubTree, std::size_t n>
  struct StaticNonUniformTypeTree
      : public TreeProperties<IsLeaf<false>,IsUniform<false>,IsTypeUniform<true>,IsStatic<true>>
      , public std::array<SubTree,n>
  {
    using Super = std::array<SubTree, n>;

    //! The type of the childs
    using ChildType = SubTree;

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

    //! Return a reference to the i'th child of the tree
    auto& child(std::size_t i)       { return Super::operator[](i); }

    //! Return a const reference to the i'th child of the tree
    auto& child(std::size_t i) const { return Super::operator[](i); }

    //! Return the number of nodes
    static constexpr index_constant<n> degree() { return {}; }
  };


  //! Non-uniform type-tree with all sub-trees of the same type and dynamic size.
  template<class SubTree>
  struct DynamicNonUniformTypeTree
      : public TreeProperties<IsLeaf<false>,IsUniform<false>,IsTypeUniform<true>,IsStatic<false>>
      , public std::vector<SubTree>
  {
    using Super = std::vector<SubTree>;

    //! The type of the childs
    using ChildType = SubTree;

    //! Inherit the constructors from std::vector
    using Super::Super;

    //! Return a reference to the i'th child of the tree
    auto& child(std::size_t i)       { return Super::operator[](i); }

    //! Return a const reference to the i'th child of the tree
    auto& child(std::size_t i) const { return Super::operator[](i); }

    //! Return the number of nodes
    std::size_t degree() { return Super::size(); }
  };


  //! Uniform type-tree with static size.
  template<class SubTree, std::size_t n>
  struct StaticUniformTypeTree
      : public TreeProperties<IsLeaf<false>,IsUniform<true>,IsTypeUniform<true>,IsStatic<true>>
  {
    //! The type of the childs
    using ChildType = SubTree;

    //! Constructor that stores the `subTree`. Can be used for class-template-argument deduction.
    StaticUniformTypeTree(std::integral_constant<std::size_t,n>, SubTree subTree)
      : subTree_{std::move(subTree)}
    {}

    /**
     * \brief Return a reference to the i'th child of the tree
     * \param i Any integral constant or index type convertible to `std::size_t`
     **/
    template <class Index>
    SubTree& child(const Index& i)
    {
      assert(std::size_t(i) < n);
      return subTree_;
    }

    /**
     * \brief Return a const reference to the i'th child of the tree
     * \param i Any integral constant or index type convertible to `std::size_t`
     **/
    template <class Index>
    SubTree const& child(const Index& i) const
    {
      assert(std::size_t(i) < n);
      return subTree_;
    }

    //! Return the number of nodes
    static constexpr index_constant<n> degree() { return {}; }

    SubTree subTree_;
  };


  //! Uniform type-tree with dynamic size.
  template<class SubTree>
  struct DynamicUniformTypeTree
      : public TreeProperties<IsLeaf<false>,IsUniform<true>,IsTypeUniform<true>,IsStatic<false>>
  {
    //! The type of the childs
    using ChildType = SubTree;

    //! Constructor that stores the `degree` and the `subTree`.
    //! Can be used for class-template-argument deduction.
    DynamicUniformTypeTree(std::size_t degree, SubTree subTree)
      : degree_{degree}
      , subTree_{std::move(subTree)}
    {}

    /**
     * \brief Return a reference to the i'th child of the tree
     * \param i Any integral constant or index type convertible to `std::size_t`
     **/
    template <class Index>
    SubTree& child(const Index& i)
    {
      assert(std::size_t(i) < degree_);
      return subTree_;
    }

    /**
     * \brief Return a reference to the i'th child of the tree
     * \param i Any integral constant or index type convertible to `std::size_t`
     **/
    template <class Index>
    SubTree const& child(const Index& i) const
    {
      assert(std::size_t(i) < degree_);
      return subTree_;
    }

    //! Return the number of nodes
    std::size_t degree() const { return degree_; }

    std::size_t degree_;
    SubTree subTree_;
  };

}} // end namespace Dune::TypeTree2

#endif // DUNE_TYPETREE2_TYPETREETREE_HH
