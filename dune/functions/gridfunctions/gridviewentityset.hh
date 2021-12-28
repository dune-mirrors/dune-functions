// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDVIEWENTITYSET_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDVIEWENTITYSET_HH

#include <cstddef>


namespace Dune {
namespace Functions {


/**
 * \brief An entity set for all entities of given codim in a grid view.
 *
 * \ingroup FunctionUtility
 *
 * This implements the \ref Concept::EntitySet concept.
 */
template<class GV, int cd>
class GridViewEntitySet
{
public:

  using GridView = GV;
  enum {
    codim = cd
  };

  //! Type of Elements contained in this EntitySet
  using Element = typename GridView::template Codim<codim>::Entity;

  //! Type of local coordinates with respect to the Element
  using LocalCoordinate = typename Element::Geometry::LocalCoordinate;
  using GlobalCoordinate = typename Element::Geometry::GlobalCoordinate;

  using value_type = Element;

  //! A forward iterator
  using const_iterator = typename GridView::template Codim<codim>::Iterator;

  //! Same as const_iterator
  using iterator = const_iterator;

  //! Construct GridViewEntitySet for a GridView.
  GridViewEntitySet(const GridView& gv) :
    gv_(gv)
  {}

  //! Return true if `e` is contained in the EntitySet.
  bool contains(const Element& e) const
  {
    return gv_.contains(e);
  }

  //! Return number of Elements visited by an iterator.
  std::size_t size() const
  {
    return gv_.size(codim);
  }

  //! Create a begin iterator.
  const_iterator begin() const
  {
    return gv_.template begin<codim>();
  }

  //! Create an end iterator.
  const_iterator end() const
  {
    return gv_.template end<codim>();
  }

  //! Return the associated GridView.
  const GridView& gridView() const
  {
    return gv_;
  }

private:
  GridView gv_;
};


} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDVIEWENTITYSET_HH
