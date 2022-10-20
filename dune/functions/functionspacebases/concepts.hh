// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH

#include <concepts>
#include <iterator>

#include <dune/common/concept.hh>
#include <dune/common/reservedvector.hh>

#include <dune/functions/common/utility.hh>

#include <dune/functions/functionspacebases/nodes.hh>


namespace Dune {
namespace Functions {
namespace Concept {

template <class C>
concept HasResize = requires(C& c) {
    c.resize(0);
};


template<class C>
concept HasSizeMethod = requires(const C& c) {
  c.size();
};


template<class C, class I>
concept HasIndexAccess = requires(C c, I i) {
  c[i];
};


// Concept for a BasisNode in a local ansatz tree
template <class N>
concept BasisNode = std::derived_from<N, BasisNodeMixin>
&& requires(const N& node, typename N::size_type idx)
{
  typename N::size_type;

  { node.size() } -> std::convertible_to<typename N::size_type>;
  { node.localIndex(idx) } -> std::convertible_to<typename N::size_type>;
  { node.treeIndex() } -> std::convertible_to<typename N::size_type>;
};



// Concept for a LeafBasisNode in a local ansatz tree
template<class N, class GridView>
concept LeafBasisNode = N::isLeaf && BasisNode<N>
&& requires(const N& node)
{
  typename N::Element;
  typename N::FiniteElement;

  { node.element() } -> std::convertible_to<typename N::Element>;
  { node.finiteElement() } -> std::same_as<const typename N::FiniteElement&>;

  requires std::same_as<typename N::Element, typename GridView::template Codim<0>::Entity>;
};

// Concept for a PowerBasisNode in a local ansatz tree
template<class N, class GridView>
concept PowerBasisNode = N::isPower && BasisNode<N>;

// Concept for a CompositeBasisNode in a local ansatz tree
template<class N, class GridView>
concept CompositeBasisNode = N::isComposite && BasisNode<N>;

// Concept for a DynamicPowerBasisNode in a local ansatz tree
template<class GridView>
struct DynamicPowerBasisNode : Refines<BasisNode>
{
  template<class N>
  auto require(const N& node) -> decltype(
    requireBaseOf<Dune::Functions::DynamicPowerBasisNode<typename N::ChildType>, N>(),
    requireConcept<BasisTree<GridView>, typename N::ChildType>()
  );
};

namespace Impl {

template <class GridView, LeafBasisNode<GridView> N>
constexpr void checkBasisTree(const N& node) {}

template <class GridView, CompositeBasisNode<GridView> N>
constexpr void checkBasisTree(const N& node);

template <class GridView, PowerBasisNode<GridView> N>
constexpr void checkBasisTree(const N& node)
  requires requires { checkBasisTree<GridView>(node.child(0)); }{}

template <class GridView, class N, std::size_t... I>
constexpr void checkCompositeBasisTree(const N& node, std::index_sequence<I...>)
  requires requires { (checkBasisTree<GridView>(node.child(Dune::index_constant<I>{})),...); }{}

template <class GridView, CompositeBasisNode<GridView> N>
constexpr void checkBasisTree(const N& node)
  requires requires {
    checkCompositeBasisTree<GridView>(node, std::make_index_sequence<N::degree()>{});
  }{}

} // end namespace Impl


// Concept for a full local BasisTree
template<class N, class GridView>
concept BasisTree = BasisNode<N>
&& ((N::isLeaf && LeafBasisNode<N,GridView>) ||
    (N::isPower && PowerBasisNode<N,GridView>) ||
    (N::isComposite && CompositeBasisNode<N,GridView>))
&& requires(N tree)
{
  // additionally check by function recursion
  Impl::checkBasisTree<GridView>(tree);
};


template<class PB>
using MultiIndex = Dune::ReservedVector<typename PB::size_type, PB::multiIndexBufferSize>;

// Concept for a PreBasis
template<class PB, class GridView = typename PB::GridView>
concept PreBasis = requires(const PB& preBasis)
{
  typename PB::GridView;
  typename PB::size_type;
  typename PB::Node;

  requires std::same_as<typename PB::GridView, GridView>;
  requires std::convertible_to<decltype(PB::maxMultiIndexSize), typename PB::size_type>;
  requires std::convertible_to<decltype(PB::maxMultiIndexSize), typename PB::size_type>;
  requires std::convertible_to<decltype(PB::maxMultiIndexSize), typename PB::size_type>;
  requires std::convertible_to<decltype(PB::multiIndexBufferSize), typename PB::size_type>;
  requires (PB::minMultiIndexSize <= PB::maxMultiIndexSize);
  requires (PB::maxMultiIndexSize <= PB::multiIndexBufferSize);

  { preBasis.gridView() } -> std::convertible_to<typename PB::GridView>;
  { preBasis.makeNode() } -> std::convertible_to<typename PB::Node>;
  { preBasis.size() } -> std::convertible_to<typename PB::size_type>;
  { preBasis.dimension() } -> std::convertible_to<typename PB::size_type>;
  { preBasis.maxNodeSize() } -> std::convertible_to<typename PB::size_type>;

  requires BasisTree<typename PB::Node,typename PB::GridView>;
  requires requires(MultiIndex<PB> mi) {
    { preBasis.size(mi) } -> std::convertible_to<typename PB::size_type>;
  };
  requires requires(typename PB::Node node, std::vector<MultiIndex<PB>>& vec) {
    { preBasis.indices(node, vec.begin()) } -> std::forward_iterator;
  };
} && requires(PB& preBasis, GridView gridView) {
  preBasis.initializeIndices();
  preBasis.update(gridView);
};



// Concept for a LocalView
template<class V, class GlobalBasis>
concept LocalView = requires(const V& localView)
{
  typename V::size_type;
  typename V::MultiIndex;
  typename V::GlobalBasis;
  typename V::Tree;
  typename V::GridView;
  typename V::Element;

  requires std::same_as<typename V::GlobalBasis, GlobalBasis>;
  requires std::same_as<typename V::GridView, typename GlobalBasis::GridView>;
  requires std::same_as<typename V::size_type, typename GlobalBasis::size_type>;
  requires std::same_as<typename V::Element, typename GlobalBasis::GridView::template Codim<0>::Entity>;

  { localView.bound() } -> std::convertible_to<bool>;
  { localView.tree() } -> std::convertible_to<typename V::Tree>;
  { localView.size() } -> std::convertible_to<typename V::size_type>;
  { localView.maxSize() } -> std::convertible_to<typename V::size_type>;
  { localView.globalBasis() } -> std::convertible_to<const typename V::GlobalBasis&>;

  requires BasisTree<typename V::Tree, typename V::GridView>;
  requires requires(typename V::size_type idx) {
    { localView.index(idx) } -> std::convertible_to<typename V::MultiIndex>;
  };
} && requires(V& localView, typename V::Element element) {
  localView.bind(element);
  localView.unbind();
};



// Concept for a GlobalBasis
template<class B, class GridView>
concept GlobalBasis = requires(const B& basis) {
  typename B::GridView;
  typename B::size_type;
  typename B::MultiIndex;
  typename B::SizePrefix;
  typename B::LocalView;
  requires std::same_as<typename B::GridView, GridView>;

  { basis.gridView() }  -> std::convertible_to<typename B::GridView>;
  { basis.localView() } -> std::convertible_to<typename B::LocalView>;
  { basis.size() }      -> std::convertible_to<typename B::size_type>;
  { basis.dimension() } -> std::convertible_to<typename B::size_type>;

  requires LocalView<typename B::LocalView,B>;
  requires requires(typename B::SizePrefix sizePrefix) {
    { basis.size(sizePrefix) } -> std::convertible_to<typename B::size_type>;
  };

} && requires(B& basis, GridView gridView) {
  basis.update(gridView);
};



} // namespace Dune::Functions::Concept
} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CONCEPTS_HH
