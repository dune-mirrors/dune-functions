// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDVIEWENTITYSET_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_GRIDVIEWENTITYSET_HH

#include <memory>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune {

namespace Functions {


/**
 * \brief An entity set for all entities of given codim in a grid view
 *
 * \ingroup FunctionUtility
 *
 * This implements the \ref Concept::EntitySet concept.
 */
template<class GV, int cd, PartitionIteratorType pitype = PartitionIteratorType::All_Partition>
class GridViewEntitySet
{
public:

  typedef GV GridView;
  enum {
    codim = cd
  };

  //! Type of Elements contained in this EntitySet
  typedef typename GridView::template Codim<codim>::Entity Element;

  //! Type of local coordinates with respect to the Element
  typedef typename Element::Geometry::LocalCoordinate LocalCoordinate;
  typedef typename Element::Geometry::GlobalCoordinate GlobalCoordinate;

  typedef Element value_type;

  //! A forward iterator
  typedef typename GridView::template Codim<codim>::template Partition<pitype>::Iterator const_iterator;

  //! Same as const_iterator
  typedef const_iterator iterator;

  /**
   * \brief Construct GridViewEntitySet for a GridView
   *
   */
  GridViewEntitySet(const GridView& gv) :
    gv_(gv)
  {}

  //! Returns true if e is contained in the EntitySet
  bool contains(const Element& e) const
  {
    return gv_.contains(e) && partitionSet<pitype>().contains(e.partitionType());
  }

  //! Number of Elements visited by an iterator
  size_t size() const
  {
    if constexpr(pitype == PartitionIteratorType::All_Partition)
      return gv_.size(codim);
    else
      return std::count_if(begin(), end(),
        [](auto const& e) { return partitionSet<pitype>().contains(e.partitionType()); });
  }

  //! Create a begin iterator
  const_iterator begin() const
  {
    return gv_.template begin<codim,pitype>();
  }

  //! Create an end iterator
  const_iterator end() const
  {
    return gv_.template end<codim,pitype>();
  }

  //! Return the associated partition-iterator
  static constexpr PartitionIteratorType partitionIterator()
  {
    return pitype;
  }

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
