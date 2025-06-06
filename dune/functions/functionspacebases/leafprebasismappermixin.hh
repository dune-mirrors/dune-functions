// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LEAFPREBASISMAPPERMIXIN_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LEAFPREBASISMAPPERMIXIN_HH

#include <dune/common/rangeutilities.hh>

#include <dune/functions/functionspacebases/leafprebasismixin.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/mcmgmapper.hh>



namespace Dune::Functions {
namespace Impl {

struct NoTwist
{
  template <class Element>
  unsigned int operator() (const Element& element, unsigned int s, unsigned int c, unsigned int i) const
  {
    return i;
  }
};

template <class IndexSet>
struct EdgeTwist
{
  static constexpr int dim = IndexSet::dimension;

  /**
   * \param indexSet The gridView indexSet containing the elements and vertices
   * \param n  The number of DOFs per edge.
   */
  EdgeTwist (const IndexSet& indexSet, unsigned int n)
    : indexSet_(&indexSet)
    , n_(n)
  {}

  template <class Element>
  unsigned int operator() (const Element& element, unsigned int s, unsigned int c, unsigned int i) const
  {
    if (c == 1) {
      auto refElem = referenceElement(element);
      auto v0 = indexSet_->subIndex(element,refElem.subEntity(s,c,0,dim),dim);
      auto v1 = indexSet_->subIndex(element,refElem.subEntity(s,c,1,dim),dim);
      bool flip = (v0 > v1);
      return flip ? n_ - i : i;
    } else {
      return i;
    }
  }

private:
  const IndexSet* indexSet_;
  unsigned int n_;
};

} // end namespace Impl


// Helper function returning a random access range
// of global indices associated to the element. The ordering
// is derived from the LocalCoefficients object.
// Having this as a member of MCMGMapper would be nice.
// But this would introduce the LocalCoefficient in dune-grid.
// This would introduce at least a weak 'conceptual' dependency problem.
template<class GridView, class LocalCoefficients, class Twist = Impl::NoTwist>
auto subIndexRange(const Dune::MultipleCodimMultipleGeomTypeMapper<GridView>& mapper, const typename GridView::template Codim<0>::Entity& element, const LocalCoefficients& localCoefficients, Twist const& twist = {})
{
  // Here we make use of the 'hidden' (poorly documented) MCMGMapper feature to support
  // multiple DOFs per subentity. However, we do not take care for any reordering.
  return Dune::transformedRangeView(Dune::range(localCoefficients.size()), [&](auto localIndex) {
    auto localKey = localCoefficients.localKey(localIndex);
    return mapper.subIndex(element, localKey.subEntity(), localKey.codim())
      + twist(element, localKey.subEntity(), localKey.codim(), localKey.index());
  });
}



/**
 * \brief A generic MixIn class for PreBasis with flat indices computed from a mapper.
 *
 * This abstracts all index computations that can be implemented using a
 * MultipleCodimMultipleGeomTypeMapper with appropriate MCMGLayout.
 * In order to use this, you need to derive from this class and
 * pass the layout in the constructor. Then the mixin takes care
 * for all the index and size computation and the derived class
 * only needs to add the node creation.
 *
 * Be careful: This does not do any reordering of indices
 * if multiple basis functions are associated to the same
 * subentity.
 *
 * \tparam GV The grid view the basis is defined on.
 */
template<typename GV, typename Twist = Impl::NoTwist>
class LeafPreBasisMapperMixin
    : public LeafPreBasisMixin<LeafPreBasisMapperMixin<GV>>
{
  static const int gridDim = GV::dimension;

public:

  //! Type of the associated GridView
  using GridView = GV;

  //! Type used for index digits
  using size_type = std::size_t;

  //! Construct from GridView and local DOF layout
  LeafPreBasisMapperMixin(const GridView& gv, Dune::MCMGLayout layout, const Twist& twist = {}) :
    gridView_(gv),
    mapper_(gridView_, std::move(layout)),
    twist_(twist)
  {}

  //! Initialize the global index information
  void initializeIndices()
  {
    // Determine upper bound for node size by traversing
    // the layout for all element types in the grid
    maxNodeSize_ = 0;
    for(const GeometryType& elementType : gridView_.indexSet().types(0))
    {
      auto referenceElement = Dune::referenceElement<double, gridDim>(elementType);
      for(auto codim : Dune::range(gridDim+1))
        for(auto i : Dune::range(referenceElement.size(codim)))
          maxNodeSize_ += mapper_.layout()(referenceElement.type(i, codim), gridDim);
    }
  }

  //! Export the stored GridView
  const GridView& gridView() const
  {
    return gridView_;
  }

  //! Update the stored GridView
  void update(const GridView& gv)
  {
    gridView_ = gv;
    mapper_.update(gridView_);
  }

  //! Return total number of basis functions
  size_type dimension() const
  {
    return mapper_.size();
  }

  //! Return maximal number of basis functions per element
  size_type maxNodeSize() const
  {
    return maxNodeSize_;
  }

  //! Fill cache with global indices of DOFs associated to the given bound node
  template<class Node, class It>
  It indices(const Node& node, It it) const
  {
    for(const auto& globalIndex : subIndexRange(mapper_, node.element(), node.finiteElement().localCoefficients(), twist_))
    {
      *it = {{ (size_type)globalIndex }};
      ++it;
    }
    return it;
  }

protected:
  GridView gridView_;
  Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper_;
  Twist twist_;
  std::size_t maxNodeSize_;
};



} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LEAFPREBASISMAPPERMIXIN_HH
