// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMEDGRIDVIEWPREBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMEDGRIDVIEWPREBASIS_HH

#include <concepts>
#include <functional>
#include <optional>
#include <utility>

namespace Dune{
namespace Functions {

// *****************************************************************************
// *****************************************************************************

/**
 * \brief A pre-basis with a transformed grid view.
 *
 * This class provides a mechanism to wrap a pre-basis on a transformed grid view.
 * It maps entities of the transformed grid view to entities of the underlying
 * basis grid view. This is particularly useful for extending or restricting
 * the domain of a pre-basis constructed on a sub-grid or super-grid, such as in
 * multi-domain grids (e.g., using dune-multidomaingrid).
 *
 * The transformation of the entity must be consistent with the domain of the underlying
 * pre-basis. This means that geometry transformations are possible, but care must be taken
 * on the underlying pre-basis so that its functions still make sense. For instance, this class
 * is sufficient to transform the geometry of affine equivalent finite element basis (e.g. Lagrange),
 * but not of affine-interpolation equivalent finite elements (e.g. Hermite, Argyris).
 *
 * \tparam RPB Raw Basis to be wrapped
 * \tparam TGV Transformed grid view type
 * \tparam GVIT Grid view inverse transformation type
 * \tparam EIT Element inverse transformation type
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 */
template <class RPB,
          class TGV = typename RPB::GridView,
          class GVIT = std::function<typename RPB::GridView(const TGV &)>,
          class EIT = std::function<typename RPB::Node::Element(const typename RPB::GridView &, const TGV &, const typename TGV::template Codim<0>::Entity &)>>
class TransformedGridViewPreBasis
{
  using GridViewInverseTransformation = GVIT;
  using ElementInverseTransformation = EIT;
  static_assert(std::is_invocable_v<GridViewInverseTransformation, const TGV &>,
              "GridViewInverseTransformation must be invocable with a TGV");
  static_assert(std::convertible_to<std::invoke_result_t<GridViewInverseTransformation, const TGV &>, const typename RPB::GridView &>,
              "GridViewInverseTransformation must return a GridView from the underlying basis");
  static_assert(std::is_invocable_v<ElementInverseTransformation, const typename RPB::GridView &, const TGV &, const typename TGV::template Codim<0>::Entity &>,
              "ElementInverseTransformation must be invocable with a GridView, TGV and an entity of the transformed grid view");
  static_assert(std::convertible_to<std::invoke_result_t<ElementInverseTransformation, const typename RPB::GridView &, const TGV &, const typename TGV::template Codim<0>::Entity &>, const typename RPB::Node::Element &>,
              "ElementInverseTransformation must return an element of the underlying basis");

public:
  //! The raw pre-basis that is being transformed
  using RawPreBasis = RPB;

  //! The grid view that the FE basis is defined on
  using GridView = TGV;

  //! Type used for indices and size information
  using size_type = typename RawPreBasis::size_type;

  static constexpr size_type maxMultiIndexSize = RawPreBasis::maxMultiIndexSize;
  static constexpr size_type minMultiIndexSize = RawPreBasis::minMultiIndexSize;
  static constexpr size_type multiIndexBufferSize = RawPreBasis::multiIndexBufferSize;

  //! Template mapping root tree path to type of created tree node
  class Node;

  /**
   * \brief Constructor from another pre-basis and a new grid view
   *
    * \param basis The raw pre-basis to be transformed
    * \param gridView The new grid view to be used
    * \param gridViewInverseTransformation A transformation from the transformed grid view to the underlying grid view
    * \param elementInverseTransformation A transformation of entities from the transformed grid view to the underlying grid view
   */
  TransformedGridViewPreBasis(
    const RawPreBasis &pre_basis,
    const GridView &gridView,
    GridViewInverseTransformation gridViewTransformation = {},
    ElementInverseTransformation elementTransformation = {})
  : rawPreBasis_{pre_basis},
    gridViewInverseTransformation_{std::move(gridViewTransformation)},
    elementInverseTransformation_{std::move(elementTransformation)},
    gridView_{gridView}
  {}

  //! Obtain the underlying untransformed pre-basis
  const RawPreBasis &rawPreBasis() const
  {
    return rawPreBasis_;
  }

  //! Obtain the underlying untransformed pre-basis
  RawPreBasis &rawPreBasis()
  {
    return rawPreBasis_;
  }

  //! Initialize the global indices
  const TransformedGridViewPreBasis& initializeIndices()
  {
    rawPreBasis().initializeIndices();
    return *this;
  }

  //! Return number of possible values for next position in empty multi index
  size_type size() const
  {
    return rawPreBasis().size();
  }

  //! Return number possible values for next position in multi index
  template<typename SizePrefix>
  size_type size(const SizePrefix &prefix) const
  {
    return rawPreBasis().size(prefix);
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return rawPreBasis().dimension();
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    return rawPreBasis_.maxNodeSize();
  }

  //! Obtain the transformed grid view that the basis is defined on
  const GridView &gridView() const
  {
    return gridView_;
  }

  //! Fill the multi-index container with the indices of the node
  template <typename It>
  It indices(const Node &node, It it) const
  {
    if (node.size() != 0)
      it = rawPreBasis().indices(node, it);
    return it;
  }

   //! Create tree node
  Node makeNode() const
  {
    return Node{rawPreBasis().makeNode(), *this};
  }

  //! Update the stored grid view
  void update(const GridView &gv)
  {
    gridView_ = gv;
    rawPreBasis().update(gridViewInverseTransformation_(gridView()));
  }

  //! Transform an entity from the transformed grid view to the underlying grid view
  std::convertible_to<typename RawPreBasis::Node::Element>
  decltype(auto) elementTransformation(const typename GridView::template Codim<0>::Entity &element) const
  {
    return elementInverseTransformation_(rawPreBasis().gridView(), gridView(), element);
  }

private:
  RawPreBasis rawPreBasis_;
  GridViewInverseTransformation gridViewInverseTransformation_;
  ElementInverseTransformation elementInverseTransformation_;
  GridView gridView_;
};


template <class RPB, class TGV, class GVIT, class EIT>
class TransformedGridViewPreBasis<RPB, TGV, GVIT, EIT>::Node : public RawPreBasis::Node
{
public:
  using RawNode = typename RawPreBasis::Node;
  using Element = typename TGV::template Codim<0>::Entity;

  Node(typename RawPreBasis::Node &&node, const TransformedGridViewPreBasis &preBasis)
    : RawPreBasis::Node(std::move(node)), preBasis_{&preBasis}
  {
  }

  //! Obtain the underlying untransformed node
  const RawNode &rawNode() const
  {
    return *this;
  }

  //! Obtain the underlying untransformed node
  RawNode &rawNode()
  {
    return *this;
  }


  const Element &element() const
  {
    return *element_;
  }

  //! Bind the node to a transformed grid element
  void bind(const Element &e)
  {
    element_ = &e;
    decltype(auto) raw_element = preBasis_->elementTransformation(element());
    typename RawNode::Element const *raw_element_view = nullptr;
    if constexpr (std::is_rvalue_reference_v<decltype(raw_element) &&>)
      raw_element_view = &raw_element_.emplace(std::move(raw_element));
    else
      raw_element_view = &raw_element;
    rawNode().bind(*raw_element_view);
  }

  //! Unbind from the current element
  void unbind()
  {
    rawNode().unbind();
    element_ = nullptr;
    raw_element_.reset();
  }

private:
  const TransformedGridViewPreBasis *preBasis_;
  std::optional<const typename RawNode::Element> raw_element_;
  Element const *element_;
};

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMEDGRIDVIEWPREBASIS_HH
