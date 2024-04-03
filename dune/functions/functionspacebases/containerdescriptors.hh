// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONTAINERDESCRIPTORS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONTAINERDESCRIPTORS_HH

#include <array>
#include <cassert>
#include <functional>
#include <type_traits>
#include <vector>

#include <dune/common/filledarray.hh>
#include <dune/common/tuplevector.hh>
#include <dune/common/typeutilities.hh>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/basistags.hh>

/**
 * \file containerdescriptor.hh
 * \brief Lightweight representation of (hierarchical) size and block structure extracted
 * from a bsis to describe data structures like containers that can be accessed by
 * multi-indices provided by the basis.
 *
 * The structure of a container-descriptor is a reduced container interface:
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

namespace Dune::Functions {
namespace ContainerDescriptors {

//! Fallback container descriptor if nothing else fits
struct Unknown {};

} // end namespace ContainerDescriptors

namespace Impl {

template<class PreBasis>
auto containerDescriptorImpl(const PreBasis& preBasis, Dune::PriorityTag<1>)
  -> decltype(preBasis.containerDescriptor())
{
  return preBasis.containerDescriptor();
}

template<class PreBasis>
auto containerDescriptorImpl(const PreBasis& preBasis, Dune::PriorityTag<0>)
{
  return ContainerDescriptors::Unknown{};
}

} // end namespace Impl

//! Return the container descriptor of the pre-basis, if defined, otherwise ContainerDescriptor::Unknown
template<class PreBasis>
auto containerDescriptor(const PreBasis& preBasis)
{
  return Impl::containerDescriptorImpl(preBasis, Dune::PriorityTag<2>{});
}


namespace ContainerDescriptors {

//! The node in the descriptor tree representing a value placeholder
struct Value
{
  //! The child access method is only available for the interface, but should not be called.
  template<class Index>
  Value operator[] (const Index&) const { return {}; }

  //! A value placeholder does not have any sub-descriptors, thus its size is zero.
  static constexpr std::size_t size () { return 0; }
};

//! Descriptor with all children of possibly different type
template<class... Children>
using Tuple = Dune::TupleVector<Children...>;

//! Generate a descriptor in case the children are not all of the same type.
//! \relates Tuple
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
using Array = std::array<Child, n>;

//! Generate a descriptor in case the children are all of the same type.
template<class Child0, class... Children,
  std::enable_if_t<(std::is_same_v<Child0, Children> &&...), int> = 0>
auto makeDescriptor (Child0 child, Children... children)
{
  using Descriptor = Array<Child0,1+sizeof...(Children)>;
  return Descriptor{std::move(child),std::move(children)...};
}


//! Descriptor for vectors with all children of the same type and dynamic size.
template<class Child>
using Vector = std::vector<Child>;

//! Descriptor for arrays with all children identical and the number of children a static size.
template<class Child, std::size_t n>
struct UniformArray
{
  //! Default constructor. Is enable if the child-type is default constructible.
  template<class C = Child,
    std::enable_if_t<std::is_default_constructible_v<C>, int> = 0>
  UniformArray ()
    : child_{}
  {}

  //! Constructor that stores a single child only.
  explicit UniformArray (Child child)
    : child_{std::move(child)}
  {}

  //! Access the i'th child that is always the same, i.e., `child_`.
  template<class Index>
  const Child& operator[] (const Index& /*i*/) const { return child_; }

  //! The static size information, i.e., number of children.
  static constexpr std::size_t size () { return n; }

private:
  Child child_;
};

//! Alias for a uniform array storing value placeholders
template<std::size_t n>
using FlatArray = UniformArray<Value,n>;

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
  //! Default constructor with size. Is enable if the child-type is default constructible.
  template<class C = Child,
    std::enable_if_t<std::is_default_constructible_v<C>, int> = 0>
  explicit UniformVector (std::size_t size)
    : size_{size}
    , child_{}
  {}

  //! Constructor that stores the size and a single child only.
  UniformVector (std::size_t size, Child child)
    : size_{size}
    , child_{std::move(child)}
  {}

  //! Access the i'th child that is always the same, i.e., `child_`.
  template<class Index>
  const Child& operator[] (const Index& /*i*/) const { return child_; }

  //! The dynamic size information, i.e., number of children.
  std::size_t size () const { return size_; }

private:
  std::size_t size_;
  Child child_;
};

//! Alias for a uniform vector storing value placeholders
using FlatVector = UniformVector<Value>;

//! Generate a uniform descriptor in case the size is a dynamic value
template<class Child>
auto makeUniformDescriptor (std::size_t n, Child child)
{
  return UniformVector<Child>{n,std::move(child)};
}


namespace Impl {

template<class InnerFunc, class LeafFunc>
struct TreeTransform
{
  TreeTransform (const InnerFunc& innerFunc, const LeafFunc& leafFunc)
    : innerFunc_(innerFunc)
    , leafFunc_(leafFunc)
  {}

  Unknown operator() (const Unknown& tree) const
  {
    return tree;
  }

  auto operator() (const Value& tree) const
  {
    return leafFunc_(tree);
  }

  template<class... V>
  auto operator() (const Tuple<V...>& tree) const
  {
    return unpackIntegerSequence([&](auto... ii) {
      return Tuple{innerFunc_(tree,ii)...};
    }, std::make_index_sequence<sizeof...(V)>());
  }

  template<class V, std::size_t n>
  auto operator() (const Array<V,n>& tree) const
  {
    return unpackIntegerSequence([&](auto... ii) {
      return Array{innerFunc_(tree,ii)...};
    }, std::make_index_sequence<n>());
  }

  template<class V>
  auto operator() (const Vector<V>& tree) const
  {
    using W = decltype(innerFunc_(tree,0));
    Vector<W> result;
    result.reserve(tree.size());
    for (std::size_t i = 0; i < tree.size(); ++i)
      result.emplace_back(innerFunc_(tree,i));
    return result;
  }

  template<class V, std::size_t n>
  auto operator() (const UniformArray<V,n>& tree) const
  {
    using W = decltype(innerFunc_(tree,Indices::_0));
    return UniformArray<W,n>(innerFunc_(tree,Indices::_0));
  }

  template<class V>
  auto operator() (const UniformVector<V>& tree) const
  {
    using W = decltype(innerFunc_(tree,0));
    return UniformVector<W>(tree.size(), innerFunc_(tree,0));
  }

private:
  InnerFunc innerFunc_;
  LeafFunc leafFunc_;
};


/**
  * Append a size to the inner-most node of the tree
  *
  * This transf of the given tree is used to implement
  * a blocked-interleaved index-merging strategy in a power-basis.
  *
  * Examples:
  * append( Flat[Container] it, size ) -> Uniform[Container]( it.size(), Flat[Container](size) )
  * append( Descriptor(child...), size ) -> Descriptor( append(child, size)... )
  */
template<class T, class Size>
auto appendToTree (const T& tree, Size s)
{
  auto transform = TreeTransform(
    [s](auto&& node, auto i) { return appendToTree(node[i], s); },
    [s](auto&& node) { return makeUniformDescriptor(s, node); });
  return transform(tree);
}


template<class AccessFunc>
struct MergeTree
{
  MergeTree (const AccessFunc& innerFunc)
    : innerFunc_(innerFunc)
  {}

  template<class... V, std::size_t m>
  auto operator() (const UniformArray<Tuple<V...>,m>& tree) const
  {
    return unpackIntegerSequence([&](auto... ii) {
      return Tuple{innerFunc_(tree,ii)...};
    }, std::make_index_sequence<(sizeof...(V)*m)>());
  }

  template<class V, std::size_t n, std::size_t m>
  auto operator() (const UniformArray<Array<V,n>,m>& tree) const
  {
    return unpackIntegerSequence([&](auto... ii) {
      return Array{innerFunc_(tree,ii)...};
    }, std::make_index_sequence<n*m>());
  }

  template<class V, std::size_t n>
  auto operator() (const UniformVector<Array<V,n>>& tree) const
  {
    using W = decltype(innerFunc_(tree,Indices::_0));
    Vector<W> result;
    result.reserve(n*tree.size());
    for (std::size_t i = 0; i < n*tree.size(); ++i)
      result.emplace_back(innerFunc_(tree,i));
    return result;
  }

  template<class V, std::size_t m>
  auto operator() (const UniformArray<Vector<V>,m>& tree) const
  {
    using W = decltype(innerFunc_(tree,0));
    Vector<W> result;
    result.reserve(tree[0].size()*m);
    for (std::size_t i = 0; i < tree[0].size()*m; ++i)
      result.emplace_back(innerFunc_(tree,i));
    return result;
  }

  template<class V>
  auto operator() (const UniformVector<Vector<V>>& tree) const
  {
    using W = decltype(innerFunc_(tree,0));
    Vector<W> result;
    result.reserve(tree[0].size()*tree.size());
    for (std::size_t i = 0; i < tree[0].size()*tree.size(); ++i)
      result.emplace_back(innerFunc_(tree,i));
    return result;
  }

  template<class V, std::size_t n, std::size_t m>
  auto operator() (const UniformArray<UniformArray<V,n>,m>& tree) const
  {
    using W = decltype(innerFunc_(tree,Indices::_0));
    return makeUniformDescriptor(std::integral_constant<std::size_t,n*m>{},
      innerFunc_(tree,Indices::_0));
  }

  template<class V, std::size_t n>
  auto operator() (const UniformVector<UniformArray<V,n>>& tree) const
  {
    using W = decltype(innerFunc_(tree,Indices::_0));
    return makeUniformDescriptor(n*tree.size(), innerFunc_(tree,Indices::_0));
  }

  template<class V, std::size_t m>
  auto operator() (const UniformArray<UniformVector<V>,m>& tree) const
  {
    using W = decltype(innerFunc_(tree,0));
    return makeUniformDescriptor(tree[0].size()*m,innerFunc_(tree,Indices::_0));
  }

  template<class V>
  auto operator() (const UniformVector<UniformVector<V>>& tree) const
  {
    using W = decltype(innerFunc_(tree,0));
    return makeUniformDescriptor(tree[0].size()*tree.size(),innerFunc_(tree,Indices::_0));
  }

private:
  AccessFunc innerFunc_;
};


// Merge `n` identical `rees
template<class IMS, class T, class Size>
auto mergeIdenticalTrees (const T& tree, Size n)
{
  auto transform = MergeTree(
    [&](auto&&, auto i) { return FlatIndexAccess::getEntry<IMS>(tree,i); });
  return transform(makeUniformDescriptor(n,tree));
}



} // end namespace Impl
} // end namespace ContainerDescriptors
} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONTAINERDESCRIPTORS_HH
