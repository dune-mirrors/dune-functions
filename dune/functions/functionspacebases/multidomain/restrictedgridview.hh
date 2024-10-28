// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_RESTRICTEDGRIDVIEW_HH
#define DUNE_FUNCTIONS_RESTRICTEDGRIDVIEW_HH

#include <typeinfo>

#include <dune/common/std/type_traits.hh>
#include <dune/common/iteratorrange.hh>
#include <dune/common/parallel/future.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/rangegenerators.hh>

// implementation adopted from Simon
// #include "subdomainindexset.hh"
#include "restrictedindexset.hh"

namespace Dune::Functions::MultiDomain
{

/**
   \brief restrict a GridView to a subdomain

   Given a GridView we keep the original iterators and everything,
   but we restrict the IndexSet to only consider those entities with
   support in the subdomain.
 */
template<typename GridView>
class RestrictedGridView : public GridView
{
public:
  /** \brief type of the index set */
  using IndexSet = RestrictedIndexSet<GridView>;

  RestrictedGridView(const GridView & gridView, std::shared_ptr<IndexSet> indexSet) :
    GridView(gridView), _restrictedIndexSet(*indexSet)
  {}

  RestrictedGridView(const RestrictedGridView&) = default;

  RestrictedGridView& operator=(const RestrictedGridView& other)
  {
    (GridView&)*this = other;
    // the indexSet should be implicitly updated
    return *this;
  }

  /** \brief obtain the index set
   *
   * The lifetime of the returned index set is bound to the lifetime of the
   * grid view. Keep a copy of the grid view to prevent the index set from
   * becoming a dangling reference.
   */
  const IndexSet &indexSet () const
  {
    return _restrictedIndexSet;
  }

  /** \brief obtain number of entities in a given codimension */
  int size ( int codim ) const
  {
    return _restrictedIndexSet.size( codim );
  }

  /** \brief obtain number of entities with a given geometry type */
  int size ( const GeometryType &type ) const
  {
    return _restrictedIndexSet.size( type );
  }

  /** @brief Return true if the given entity is contained in this grid view
   * @todo Currently we call the implementation on the IndexSet.  This may lead to suboptimal efficiency.
   *
   * \note If the input element e is not an element of the grid, then
   *       the result of contains() is undefined.
   */
  template<class EntityType>
  bool contains (const EntityType& e) const
  {
    return _restrictedIndexSet.contains(e);
  }

private:
  RestrictedGridView() = delete;
  IndexSet _restrictedIndexSet;
};

} // namespace Dune::Functions::MultiDomain

#endif // DUNE_FUNCTIONS_RESTRICTEDGRIDVIEW_HH
