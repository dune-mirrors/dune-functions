// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONTAINERDESCRIPTORS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONTAINERDESCRIPTORS_HH

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
 * \file containerdescriptor.hh
 * \brief Lightweight representation of (hierarchical) size and block structure extracted
 * from a bsis, to describe data structures like containers that can be accessed by
 * the multi-indices provided by the basis.
 *
 * The structure of an container-descriptor is a reduced container interface:
 * \code
  struct [Container]Descriptor
  {
    template<class Index>
    [SubContainerDescriptor] operator[](Index i) const;  // return the i-th sub-container-descriptor

    [static constexpr] std::size_t size() [const];  // return the number of children
  };
 * \endcode
 *
 * With the `operator[]` you can access the children. The `Index` type is either
 * an integral value or an `integral_constant` for tuple nodes.
 *
 * Size is either a static property, or a runtime value.
 **/

namespace Dune::Functions::ContainerDescriptors {


  template <class D>
  struct IsTypeUniform : std::false_type {};

  template <class D>
  inline constexpr bool isTypeUniform = IsTypeUniform<D>::value;

  template <class D>
  struct IsUniform : std::false_type {};

  template <class D>
  inline constexpr bool isUniform = IsUniform<D>::value;

  template <class D>
  struct IsFlat : std::false_type {};

  template <class D>
  inline constexpr bool isFlat = IsFlat<D>::value;


  struct Unknown {};

  //! The node in the descriptor tree representing a value placeholder
  struct Value
  {
    template <class Index>
    Value operator[] (const Index&) const { return {}; }

    static constexpr std::size_t size () { return 0; }
  };

  template <>
  struct IsTypeUniform<Value> : std::true_type {};

  template <>
  struct IsUniform<Value> : std::true_type {};


  //! Descriptor with all children of different type
  template<class... Children>
  struct Tuple
      : private Dune::TupleVector<Children...>
  {
    using Super = Dune::TupleVector<Children...>;

    template<class C = std::tuple<Children...>,
      std::enable_if_t<std::is_default_constructible_v<C>, int> = 0>
    Tuple ()
      : Tuple{Children{}...}
    {}

    explicit Tuple (Children... children)
      : Super{std::move(children)...}
    {}

    using Super::operator[];
    using Super::size;
  };

  //! Generate a descriptor in case the children are not all of the same type.
  template<class Child0, class... Children,
    std::enable_if_t<(sizeof...(Children) > 0), int> = 0,
    std::enable_if_t<(...|| (not std::is_same_v<Child0, Children>)), int> = 0>
  auto makeDescriptor (Child0 child0, Children... children)
  {
    using Descriptor = Tuple<Child0,Children...>;
    return Descriptor{std::move(child0),std::move(children)...};
  }


  //! Descriptor for arrays with all children of the same type and static size.
  template<class Child, std::size_t n>
  struct Array
      : private std::array<Child, n>
  {
    using Super = std::array<Child, n>;

    template<class C = Child,
      std::enable_if_t<std::is_default_constructible_v<C>, int> = 0>
    Array ()
      : Array{C{}}
    {}

    explicit Array (Child child)
      : Super{Dune::filledArray<n>(std::move(child))}
    {}

    Array (std::integral_constant<std::size_t,n>, Child child)
      : Super{Dune::filledArray<n>(std::move(child))}
    {}

    template <class... Children,
      std::enable_if_t<(std::is_same_v<Children,Child> &&...), int> = 0>
    explicit Array (Children... children)
      : Super{std::move(children)...}
    {}

    using Super::operator[];

    static constexpr std::size_t size () { return n; }
  };

  template <class C, std::size_t n>
  struct IsTypeUniform<Array<C,n>>
    : std::true_type {};


  //! Generate a descriptor in case the children are all of the same type.
  template<class Child0, class... Children,
    std::enable_if_t<(...&& std::is_same_v<Child0, Children>), int> = 0>
  auto makeDescriptor (Child0 child, Children... children)
  {
    using Descriptor = Array<Child0,1+sizeof...(Children)>;
    return Descriptor{std::move(child),std::move(children)...};
  }


  //! Descriptor for a vector with all children of the same type and dynamic size.
  template<class Child>
  struct Vector
      : private std::vector<Child>
  {
    using Super = std::vector<Child>;
    using Super::Super;
    using Super::operator[];
    using Super::size;
  };

  template <class C>
  struct IsTypeUniform<Vector<C>>
    : std::true_type {};


  //! Descriptor for an array with all children identical and the number of children a static size.
  template<class Child, std::size_t n>
  struct UniformArray
  {
    template<class C = Child,
      std::enable_if_t<std::is_default_constructible_v<C>, int> = 0>
    UniformArray ()
      : child_{}
    {}

    explicit UniformArray (Child child)
      : child_{std::move(child)}
    {}

    UniformArray (std::integral_constant<std::size_t,n>, Child child)
      : child_{std::move(child)}
    {}

    template <class Index>
    const Child& operator[] (const Index&) const { return child_; }

    static constexpr std::size_t size () { return n; }

  private:
    Child child_;
  };

  template <class C, std::size_t n>
  struct IsTypeUniform<UniformArray<C,n>>
    : std::true_type {};

  template <class C, std::size_t n>
  struct IsUniform<UniformArray<C,n>>
    : std::true_type {};


  template <std::size_t n>
  using FlatArray = UniformArray<Value,n>;

  template <std::size_t n>
  struct IsFlat<FlatArray<n>> : std::true_type {};


  //! Generate a uniform descriptor in case the size is a static constant
  template<class Child, std::size_t n>
  auto makeUniformDescriptor (std::integral_constant<std::size_t,n>, Child child)
  {
    return UniformArray<Child,n>{std::move(child)};
  }


  //! Uniform descriptor with dynamic size.
  template<class Child>
  struct UniformVector
  {
    template<class C = Child,
      std::enable_if_t<std::is_default_constructible_v<C>, int> = 0>
    UniformVector ()
      : size_{0}
      , child_{}
    {}

    template<class C = Child,
      std::enable_if_t<std::is_default_constructible_v<C>, int> = 0>
    explicit UniformVector (std::size_t size)
      : size_{size}
      , child_{}
    {}

    UniformVector (std::size_t size, Child child)
      : size_{size}
      , child_{std::move(child)}
    {}

    template <class Index>
    const Child& operator[] (const Index&) const { return child_; }

    std::size_t size () const { return size_; }

  private:
    std::size_t size_;
    Child child_;
  };

  template <class C>
  struct IsTypeUniform<UniformVector<C>>
    : std::true_type {};

  template <class Child>
  struct IsUniform<UniformVector<Child>>
    : std::true_type {};


  using FlatVector = UniformVector<Value>;

  template <>
  struct IsFlat<FlatVector> : std::true_type {};

  //! Generate a uniform descriptor in case the size is a dynamic value
  template<class Child>
  auto makeUniformDescriptor (std::size_t n, Child child)
  {
    return UniformVector<Child>{n,std::move(child)};
  }


  namespace Impl {

  // ---------------------------------------------------------------------------
  // Some utilities for generating container descriptors
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

    // overload for Unknown node. Cannot access any component.
    template<class FlatIndex, class OuterOffsetIndex, class IMS>
    static void getEntry (Unknown, FlatIndex, OuterOffsetIndex, IMS) {}

    // overload for Value nodes. Cannot access any component.
    template<class FlatIndex, class OuterOffsetIndex, class IMS>
    static void getEntry (Value, FlatIndex, OuterOffsetIndex, IMS) {}

    // overload for flat-lexicographic index-merging strategy
    template<class Tree, class FlatIndex, class OuterOffsetIndex>
    static auto getEntry (const Tree& tree, FlatIndex i, OuterOffsetIndex o,
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
    template<class Tree, class FlatIndex, class InnerIndex>
    static auto getEntry (const Tree& tree, FlatIndex i, InnerIndex o,
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
    // Return the sub-node in a given tree associated to a flat index
    // after merging the nodes with a given index-merging strategy
    template<class IMS, class Tree, class FlatIndex>
    static auto getEntry (const Tree& tree, FlatIndex i)
    {
      return getEntry(tree,i,Indices::_0,IMS{});
    }
  };


  /*
   * Generic implementation of the tree merging. The trees to merge are
   * collected into an artificial super tree `tree`. The number of nodes
   * after merging the sub-trees of `tree` is given by the parameter `size`. And
   * the tree properties as `isUniform` and `isTypeUniform`.
   *
   * \tparam IMS  The index-merging strategy used as strategy for merging the tree
   *    nodes, either BasisFactory::FlatLexicographic or BasisFactory::FlatInterleaved.
   *
   * \param tree  Collection of nodes to be merged
   * \param size  The size of the resulting container
   * \param isUniform  The resulting descriptor is uniform (all children identical)
   * \param isTypeUniform  The resulting descriptor is a container (all children have the same type)
   */
  template<class IMS, class Tree, class Size, bool isUniform, bool isTypeUniform>
  auto mergeTreesImpl (const Tree& tree, Size size,
                       std::bool_constant<isUniform>,
                       std::bool_constant<isTypeUniform>)
  {
    // access one of the sub-nodes of the nodes in `tree` by a flat index `ii`
    auto child = [&](auto ii) { return FlatIndexAccess::getEntry<IMS>(tree,ii); };
    using Child00 = std::decay_t<decltype(tree[Indices::_0][Indices::_0])>;

    if constexpr(isUniform)
      return makeUniformDescriptor(size, tree[Indices::_0][Indices::_0]);
    else if constexpr(IsIntegralConstant<Size>::value)
      return unpackIntegerSequence([&](auto... ii) {
        return makeDescriptor(child(ii)...);
      }, std::make_index_sequence<Size::value>{});
    else if constexpr(isTypeUniform) {
      Vector<Child00> result(size);
      for (std::size_t i = 0; i < size; ++i)
        result[i] = std::move(child(i));
      return result;
    }
    else {
      DUNE_THROW(Dune::NotImplemented,
        "Merging of dynamic non-uniform trees not implemented");
      return Unknown{};
    }
  }


  // Overload for zero trees, return an unknown tree.
  template<class IMS>
  inline void mergeTrees () {}

  // Overload for one tree, return the tree itself.
  template<class IMS, class T>
  auto mergeTrees (const T& tree)
  {
    return tree;
  }

  // Merge a variadic list of trees
  template<class IMS, class T0, class... Ts,
    std::enable_if_t<(sizeof...(Ts) > 0), int> = 0>
  auto mergeTrees (const T0& tree0, const Ts&... trees)
  {
    auto and_ = Hybrid::hybridFunctor(LogicalAnd{});
    auto or_ = Hybrid::hybridFunctor(std::logical_or<>{});
    auto equal_ = Hybrid::hybridFunctor(std::equal_to<>{});

    using T00 = std::decay_t<decltype(tree0[Indices::_0])>;

    // collect some properties of the passed trees
    const bool allUniform = (isUniform<T0> &&...&& isUniform<Ts>);
    const bool allTypeUniform = (isTypeUniform<T0> &&...&& isTypeUniform<Ts>);
    const bool allSubTypeUniform = (...&&
      std::is_same_v<T00, std::decay_t<decltype(trees[Indices::_0])>>);
    const auto allSameSize = and_(equal_(Hybrid::size(tree0), Hybrid::size(trees))...);

    // The resulting tree can only be uniform if we can deduce at compile-time that
    // all nodes of all nodes are identical. This is only possible in some cases,
    // e.g. if all nodes are `Value`s.
    const bool allAncestorsAreValues = (std::is_same_v<T00, Value> &&...&&
      std::is_same_v<std::decay_t<decltype(trees[Indices::_0])>, Value>);

    // The merge of variadic trees currently is only implemented with lexicographic merge
    // or interleaved merge if all nodes have the same size
    const bool isFlatLexicographic = std::is_same_v<IMS,BasisFactory::FlatLexicographic>;
    hybridAssert(or_(std::bool_constant<isFlatLexicographic>{}, allSameSize));

    auto sumSizes = Hybrid::plus(Hybrid::size(tree0),Hybrid::size(trees)...);
    return mergeTreesImpl<IMS>(makeDescriptor(tree0,trees...),
      sumSizes,
      std::bool_constant<allUniform && allAncestorsAreValues>{},
      std::bool_constant<allTypeUniform && allSubTypeUniform>{}
    );
  }

  // Merge `n` identical `rees
  template<class IMS, class Size, class T,
    std::enable_if_t<std::is_convertible_v<Size,std::size_t>, int> = 0>
  auto mergeIdenticalTrees (Size n, const T& tree)
  {
    auto multiplies = Hybrid::hybridFunctor(std::multiplies<>{});
    auto sumSizes = multiplies(Hybrid::size(tree), n);
    return mergeTreesImpl<IMS>(makeUniformDescriptor(n,tree),
      sumSizes,
      std::bool_constant<isUniform<T>>{},
      std::bool_constant<isTypeUniform<T>>{}
    );
  }

  // Merge `n` identical trees, redirects to the function with 2 arguments
  template<std::size_t n, class IMS, class T>
  auto mergeIdenticalTrees (const T& tree)
  {
    return mergeIdenticalTrees<IMS>(index_constant<n>{}, tree);
  }



  // Append the size `s` at the inner-most node of the tree
  template<class Size>
  auto appendToTree (Value, Size s)
  {
    return makeUniformDescriptor(s, Value{});
  }

  /*
   * Append a size to the inner-most node of the tree
   *
   * This transforming of the given tree is used to implement
   * a blocked-interleaved index-merging strategy in a power-basis.
   *
   * Examples:
   * append( Flat[Container] it, size ) -> Uniform[Container]( it.size(), Flat[Container](size) )
   * append( Descriptor(child...), size ) -> Descriptor( append(child, size)... )
   */
  template<class T, class Size>
  auto appendToTree (const T& tree, Size s)
  {
    auto child = [&](auto ii) { return appendToTree(tree[ii], s); };

    if constexpr(isUniform<T>)
      return makeUniformDescriptor(Hybrid::size(tree), child(Indices::_0));
    else if constexpr(HasStaticSize_v<T>)
      return unpackIntegerSequence([&](auto... ii) {
        return makeDescriptor(child(ii)...);
      }, std::make_index_sequence<std::size_t(T::size())>());
    else if constexpr(isTypeUniform<T>) {
      Vector<decltype(child(0))> result(tree.size());
      for (std::size_t i = 0; i < tree.size(); ++i)
        result[i] = std::move(child(i));
      return result;
    }
    else {
      DUNE_THROW(Dune::NotImplemented,
        "Merging of dynamic non-uniform trees not implemented");
      return Unknown{};
    }
  }

  // Append the size `s` at the inner-most node of the tree. Redirects to the
  // two-argument functions
  template<std::size_t s, class T>
  auto appendToTree (const T& tree)
  {
    return appendToTree(tree, index_constant<s>{});
  }

  } // end namespace Impl
} // end namespace Dune::Functions::ContainerDescriptors

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONTAINERDESCRIPTORS_HH
