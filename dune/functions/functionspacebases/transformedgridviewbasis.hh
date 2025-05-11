// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMEDGRIDVIEWBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMEDGRIDVIEWBASIS_HH

#include <concepts>
#include <functional>
#include <optional>
#include <utility>

namespace Dune {
namespace Functions {

// *****************************************************************************
// *****************************************************************************

/**
 * \brief A basis with a transformed grid view.
 *
 * This class provides a mechanism to wrap a basis on a transformed grid view.
 * It maps entities of the transformed grid view to entities of the underlying
 * basis grid view. This is particularly useful for extending or restricting
 * the domain of a basis constructed on a sub-grid or super-grid, such as in
 * multi-domain grids (e.g., using dune-multidomaingrid).
 *
 * The transformation of the entity must be consistent with the domain of the underlying
 * basis. This means that geometry transformations are possible, but care must be taken
 * on the underlying basis so that its functions still make sense. For instance, this class
 * is sufficient to transform the geometry of affine equivalent finite element basis (e.g. Lagrange),
 * but not of affine-interpolation equivalent finite elements (e.g. Hermite, Argyris).
 *
 * \tparam RB Raw Basis to be wrapped
 * \tparam TGV Transformed grid view type
 * \tparam GVIT Grid view inverse transformation type
 * \tparam EIT Element inverse transformation type
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 */
template <class RB,
          class TGV = typename RB::GridView,
          class GVIT = std::function<typename RB::GridView(const TGV &)>,
          class EIT = std::function<typename RB::LocalView::Element(const typename RB::GridView &, const TGV &, const typename TGV::template Codim<0>::Entity &)>>
class TransformedGridViewBasis
{
  using GridViewInverseTransformation = GVIT;
  using ElementInverseTransformation = EIT;
  static_assert(std::is_invocable_v<GridViewInverseTransformation, const TGV &>,
              "GridViewInverseTransformation must be invocable with a TGV");
  static_assert(std::convertible_to<std::invoke_result_t<GridViewInverseTransformation, const TGV &>, const typename RB::GridView &>,
              "GridViewInverseTransformation must return a GridView from the underlying basis");
  static_assert(std::is_invocable_v<ElementInverseTransformation, const typename RB::GridView &, const TGV &, const typename TGV::template Codim<0>::Entity &>,
              "ElementInverseTransformation must be invocable with a GridView, TGV and an entity of the transformed grid view");
  static_assert(std::convertible_to<std::invoke_result_t<ElementInverseTransformation, const typename RB::GridView &, const TGV &, const typename TGV::template Codim<0>::Entity &>, const typename RB::LocalView::Element &>,
              "ElementInverseTransformation must return an element of the underlying basis");
public:
  //! The raw basis that is being transformed
  using RawBasis = RB;

  //! The transformed grid view that the FE space is defined on
  using GridView = TGV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Type of the local view on the restriction of the basis to a single element
  class LocalView;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = typename RawBasis::MultiIndex;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = typename RawBasis::SizePrefix;

  /**
   * \brief Constructor from another basis and a new grid view
   *
    * \param basis The raw basis to be transformed
    * \param gridView The new grid view to be used
    * \param gridViewInverseTransformation A transformation from the transformed grid view to the underlying grid view
    * \param elementInverseTransformation A transformation of entities from the transformed grid view to the underlying grid view
   */
  TransformedGridViewBasis(
    const RawBasis &basis,
    const GridView &gridView,
    GridViewInverseTransformation gridViewInverseTransformation = {},
    ElementInverseTransformation elementInverseTransformation = {})
  : rawBasis_(basis),
    gridViewInverseTransformation_{std::move(gridViewInverseTransformation)},
    elementInverseTransformation_{std::move(elementInverseTransformation)},
    gridView_{gridView}
  {}

  //! Obtain the underlying untransformed basis
  const RawBasis &rawBasis() const
  {
    return rawBasis_;
  }

  //! Obtain the underlying untransformed basis
  RawBasis &rawBasis()
  {
    return rawBasis_;
  }

  //! Return number of possible values for next position in empty multi index
  size_type size() const
  {
    return rawBasis().size();
  }

  //! Return number possible values for next position in multi index
  size_type size(const SizePrefix &prefix) const
  {
    return rawBasis().size(prefix);
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return rawBasis().dimension();
  }

  //! Obtain the transformed grid view that the basis is defined on
  const GridView &gridView() const
  {
    return gridView_;
  }

  //! \brief Return local view for basis
  LocalView localView() const
  {
    return LocalView(rawBasis().localView(), *this);
  }

  //! Update the stored grid view
  void update(const GridView &gv)
  {
    gridView_ = gv;
    rawBasis().update(gridViewInverseTransformation_(gridView()));
  }

  //! Transform an entity from the transformed grid view to the underlying grid view
  std::convertible_to<typename RawBasis::LocalView::Element>
  decltype(auto) elementInverseTransformation(const typename GridView::template Codim<0>::Entity &element) const
  {
    return elementInverseTransformation_(rawBasis().gridView(), gridView(), element);
  }

private:
  RawBasis rawBasis_;
  GridViewInverseTransformation gridViewInverseTransformation_;
  ElementInverseTransformation elementInverseTransformation_;
  GridView gridView_;
};

//! Type of the local view on the restriction of the basis to a single element
template <class RB, class TGV, class GVIT, class EIT>
class TransformedGridViewBasis<RB,TGV,GVIT,EIT>::LocalView : public RB::LocalView
{
  using RawBasis = RB;
public:
  //! The raw basis that is being transformed
  using RawLocalView = typename RawBasis::LocalView;
  //! The transformed grid view that the FE space is defined on
  using GridView = TGV;
  //! Transformed type of the grid element we are bound to
  using Element = typename GridView::template Codim<0>::Entity;

  //! Obtain the underlying untransformed local view
  const RawLocalView &rawLocalView() const
  {
    return *this;
  }

  //! Obtain the underlying untransformed local view
  RawLocalView &rawLocalView()
  {
    return *this;
  }

  //! Constructor from the underlying local view and the transformed global basis
  LocalView(RawBasis::LocalView &&localView, const TransformedGridViewBasis &globalBasis)
      : RawBasis::LocalView(std::move(localView)), globalBasis_(&globalBasis)
  {}

  //! return the grid element that the view is bound to
  const Element &element() const
  {
    return *element_;
  }

  //! Bind the view to a transformed grid element
  void bind(const typename GridView::template Codim<0>::Entity &e)
  {
    element_ = &e;
    decltype(auto) raw_element = globalBasis_->elementInverseTransformation(element());
    typename RawLocalView::Element const *raw_element_view = nullptr;
    if constexpr (std::is_rvalue_reference_v<decltype(raw_element) &&>)
      raw_element_view = &raw_element_.emplace(std::move(raw_element));
    else
      raw_element_view = &raw_element;
    rawLocalView().bind(*raw_element_view);
  }

  //! Unbind from the current element
  void unbind()
  {
    rawLocalView().unbind();
    element_ = nullptr;
    raw_element_.reset();
  }

private:
  const TransformedGridViewBasis *globalBasis_;
  std::optional<const typename RawLocalView::Element> raw_element_;
  Element const *element_;
};

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMEDGRIDVIEWBASIS_HH
