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
#include <dune/common/typetraits.hh>

namespace Dune {
namespace TypeTree2 {

#ifndef DOXYGEN
  template <bool b> struct IsLeaf : std::bool_constant<b> {};
  template <bool b> struct IsUniform : std::bool_constant<b> {};
  template <bool b> struct IsTypeUniform : std::bool_constant<b> {};
  template <bool b> struct IsStatic : std::bool_constant<b> {};
#endif

#ifndef DOXYGEN
  template<class, class, class, class>
  struct TreeProperties;
#endif

  //! Define basic properties of a tree
  template<bool leaf, bool uniform, bool typeUniform, bool staticSize>
  struct TreeProperties<
    IsLeaf<leaf>, IsUniform<uniform>, IsTypeUniform<typeUniform>, IsStatic<staticSize> >
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
      , public Dune::TupleVector<SubTrees...>
  {
    using Super = Dune::TupleVector<SubTrees...>;

    //! The type of the childs tuple
    using ChildTypes = std::tuple<SubTrees...>;

    //! The type of the i'th child
    template <std::size_t i>
    using Child = std::tuple_element_t<i,ChildTypes>;

    //! Inherit the constructors from std::tuple
    using Super::Super;

    //! Return a reference to the i'th child of the tree
    template<std::size_t i>
    Child<i>& child(index_constant<i> ii)
    {
      return Super::operator[](ii);
    }

    //! Return a const reference to the i'th child of the tree
    template<std::size_t i>
    const Child<i>& child(index_constant<i> ii) const
    {
      return Super::operator[](ii);
    }

    //! Return the number of nodes
    static constexpr index_constant<(sizeof...(SubTrees))> degree() { return {}; }
  };

  // deduction guides
  template<class... SubTrees>
  VariadicNonUniformTypeTree(const SubTrees&...)
    -> VariadicNonUniformTypeTree<SubTrees...>;


  //! Non-uniform type-tree with all sub-trees of the same type and static size.
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
    template<class... SubTrees,
      std::enable_if_t<(1+sizeof...(SubTrees) == n), int> = 0>
    explicit StaticNonUniformTypeTree(const SubTree& subTree, SubTrees&&... subTrees)
      : Super{subTree, std::forward<SubTrees>(subTrees)...}
    {}

    //! Repeat the single argument to fill the array
    explicit StaticNonUniformTypeTree(index_constant<n>, const SubTree& subTree)
      : Super{Dune::unpackIntegerSequence(
          [&](auto... ii) {
            return Super{(void(ii), subTree)...};
          }, std::make_index_sequence<n>{})
        }
    {}

    //! Return a reference to the i'th child of the tree
    SubTree& child(std::size_t i) { return Super::operator[](i); }

    //! Return a const reference to the i'th child of the tree
    const SubTree& child(std::size_t i) const { return Super::operator[](i); }

    //! Return the number of nodes
    static constexpr index_constant<n> degree() { return {}; }
  };

  // deduction guides
  template<std::size_t n, class SubTree>
  StaticNonUniformTypeTree(index_constant<n>, const SubTree&)
    -> StaticNonUniformTypeTree<SubTree, n>;

  template<class SubTree, class... SubTrees,
    std::enable_if_t<not Dune::IsIntegralConstant<SubTree>::value, int> = 0,
    std::enable_if_t<(std::is_convertible_v<SubTrees,SubTree> &&...), int> = 0>
  StaticNonUniformTypeTree(const SubTree&, const SubTrees&...)
    -> StaticNonUniformTypeTree<SubTree, 1+sizeof...(SubTrees)>;



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
    SubTree& child(std::size_t i) { return Super::operator[](i); }

    //! Return a const reference to the i'th child of the tree
    const SubTree& child(std::size_t i) const { return Super::operator[](i); }

    //! Return the number of nodes
    std::size_t degree() { return Super::size(); }
  };

  // deduction guides
  template<class SubTree>
  DynamicNonUniformTypeTree(std::size_t, const SubTree&)
    -> DynamicNonUniformTypeTree<SubTree>;


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
     * \brief Return a reference to the stored `subTree`
     * \param i Any integral constant or index type convertible to `std::size_t`
     **/
    template <class Index>
    SubTree& child(const Index& i)
    {
      assert(std::size_t(i) < n);
      return subTree_;
    }

    /**
     * \brief Return a const reference to the stored `subTree`
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
     * \brief Return a reference to the stored `subTree`
     * \param i Any integral constant or index type convertible to `std::size_t`
     **/
    template <class Index>
    SubTree& child(const Index& i)
    {
      assert(std::size_t(i) < degree_);
      return subTree_;
    }

    /**
     * \brief Return a reference to the stored `subTree`
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
