// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASIXBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASIXBASIS_HH

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <span>
#include <type_traits>

#include <dune/localfunctions/basix/globalbasix.hh>
#include <dune/localfunctions/basix/utility.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>

#include <dune/geometry/type.hh>

namespace Dune::Functions {

/**
  * \brief A pre-basis based on the basix framework
  *
  * \ingroup FunctionSpaceBasesImplementations
  *
  * \tparam GV        The grid view that the FE basis is defined on.
  * \tparam rangeClass  The class of the basis function range, i.e., RangeClass:scalar, vector or matrix.
  */
template<class GV, RangeClass rangeClass>
class BasixPreBasis
    : public LeafPreBasisMapperMixin<GV>
{
  using Base = LeafPreBasisMapperMixin<GV>;
  using Geometry = typename GV::template Codim<0>::Entity::Geometry;
  using FiniteElement = BasixFiniteElement<Geometry,rangeClass,typename GV::ctype>;
  using Basix = typename FiniteElement::Basix;

  static auto makeLayout (const Basix& b)
  {
    auto& entity_dofs = b.entity_dofs();
    auto subentity_types =  basix::cell::subentity_types(b.cell_type());

    std::array<std::size_t, Dune::Impl::number_of_cell_types> sizes{};
    for (std::size_t d = 0; d < entity_dofs.size(); ++d) {
      for (std::size_t s = 0; s < entity_dofs[d].size(); ++s) {
        int i = d+1 < entity_dofs.size() ? (int)(subentity_types[d][s]) : (int)(b.cell_type());
        assert(i >= 0 && i < int(sizes.size()));
        assert(sizes[i] == 0 || sizes[i] == entity_dofs[d][s].size());
        sizes[i] = entity_dofs[d][s].size();
      }
    }
    return [sizes](GeometryType gt, int) -> std::size_t { return sizes[(int)(Dune::Impl::cellType(gt))]; };
  }


public:
  using GridView = GV;
  class Node;
  using size_type = typename Base::size_type;

public:
  BasixPreBasis (const GV& gridView, const Basix& basix)
    : Base(gridView, makeLayout(basix))
    , basix_{basix}
  {
    computeCellInfos(gridView);
  }

  using Base::gridView;

  //! Create tree node
  Node makeNode () const
  {
    return Node(basix_, gridView().indexSet(), cellInfos_);
  }

  //! Update the stored GridView
  void update (const GridView& gv)
  {
    Base::update(gv);
    computeCellInfos(gv);
  }

  //! Fill cache with global indices of DOFs associated to the given bound node
  template<class Node, class It>
  It indices (const Node& node, It it) const
  {
    static_assert(std::is_same_v<typename std::iterator_traits<It>::iterator_category, std::random_access_iterator_tag>);

    auto first = it;
    auto last = Base::indices(node, it);
    assert(last - first == std::ptrdiff_t(node.size()));

    if (basix_.dof_transformations_are_permutations())
    {
      // basix_.permute(std::span{&*first, node.size()}, fe_.cellInfo());
      // NOTE: cannot pass a span over the indices into permute, since indices are not int32_t
    }

    return last;
  }

protected:
  template <class IndexSet, class Element, class RefElem>
  std::uint32_t computeEdgeInfo (const IndexSet& indexSet, const Element& e, const RefElem& r, int i, int c)
  {
    int ii0 = r.subEntity(i,c,0,Element::dimension);
    int ii1 = r.subEntity(i,c,1,Element::dimension);

    int jj0 = indexSet.subIndex(e,ii0,Element::dimension);
    int jj1 = indexSet.subIndex(e,ii1,Element::dimension);
    int flipOrientation = (ii0 < ii1) != (jj0 < jj1);

    return flipOrientation;
  }

  template <class IndexSet, class Element, class RefElem>
  std::uint32_t computeFaceInfo (const IndexSet& indexSet, const Element& e, const RefElem& r, int i, int c)
  {
    std::array<int,4> ii{99,99,99,99};
    std::array<std::size_t,4> jj{std::size_t(-1),std::size_t(-1),std::size_t(-1),std::size_t(-1)};

    for (int j = 0; j < e.subEntities(i,c,Element::dimension); ++j) {
      ii[j] = r.subEntity(i,c,j,Element::dimension);
      jj[j] = indexSet.subIndex(e,ii[j],Element::dimension);
    }
    int flipOrientation = (ii[0] < ii[1]) != (jj[0] < jj[1]);

    int iimin = std::distance(ii.begin(),std::min_element(ii.begin(), ii.end()));
    int jjmin = std::distance(jj.begin(),std::min_element(jj.begin(), jj.end()));
    int rotations = std::abs(iimin - jjmin);

    return flipOrientation | rotations<<1;
  }

  void computeCellInfos (const GridView& gv)
  {
    auto& indexSet = gv.indexSet();
    cellInfos_.resize(indexSet.size(0));
    for (auto const& e : elements(gv))
    {
      std::uint32_t cellInfo = 0;
      auto r = referenceElement(e);
      int pos = 0;

      if constexpr(GridView::dimension > 1) {
        for (int i = 0; i < r.size(1); ++i)
        {
          if constexpr(GridView::dimension == 2)
            cellInfo |= computeEdgeInfo(indexSet,e,r,i,1)<<(pos++);
          else if constexpr(GridView::dimension == 3) {
            cellInfo |= computeFaceInfo(indexSet,e,r,i,1)<<pos;
            pos += 3;
          }
        }
      }

      if constexpr(GridView::dimension > 2) {
        for (int i = 0; i < r.size(2); ++i)
        {
          if constexpr(GridView::dimension == 3)
            cellInfo |= computeEdgeInfo(indexSet,e,r,i,2)<<(pos++);
        }
      }

      cellInfos_[indexSet.index(e)] = cellInfo;
    }
  }

protected:
  Basix basix_;
  std::vector<std::uint32_t> cellInfos_;
};


template <class GV, RangeClass rangeClass>
class BasixPreBasis<GV,rangeClass>::Node
  : public LeafBasisNode
{
  static constexpr int dim = GV::dimension;

public:
  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using IndexSet = typename GV::IndexSet;
  using Basix = typename BasixPreBasis::Basix;
  using FiniteElement = typename BasixPreBasis::FiniteElement;

  //! Constructor; stores a pointer to the passed local finite-element `fe`.
  explicit Node (const Basix& basix, const IndexSet& indexSet, const std::vector<std::uint32_t>& cellInfos)
    : fe_{basix}
    , indexSet_(&indexSet)
    , cellInfos_(&cellInfos)
    , element_{nullptr}
  {}

  //! Return current element; might raise an error if unbound
  const Element& element () const
  {
    assert(!!element_);
    return *element_;
  }

  /**
   * \brief Return the LocalFiniteElement for the element we are bound to; might raise an error if unbound.
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement () const
  {
    return fe_;
  }

  //! Bind to element. Stores a pointer to the passed element reference.
  void bind (const Element& e)
  {
    element_ = &e;
    fe_.bind(element_->geometry(), (*cellInfos_)[indexSet_->index(e)]);
    this->setSize(fe_.size());
  }

protected:
  FiniteElement fe_;
  const IndexSet* indexSet_;
  const std::vector<std::uint32_t>* cellInfos_;
  const Element* element_;
};


namespace BasisFactory {

/**
  * \brief A factory that can create a Basix pre-basis
  *
  * \ingroup FunctionSpaceBasesImplementations
  *
  * \tparam Basix    The type of the basix finite-element
  * \tparam dimRange Dimension of the range of the basis functions
  */
template<RangeClass rangeClass, class Basix>
auto basix (Basix basix)
{
  return [basix=std::move(basix)](const auto& gridView) {
    return BasixPreBasis<std::decay_t<decltype(gridView)>, rangeClass>(gridView, std::move(basix));
  };
}

} // end namespace BasisFactory
} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASIXBASIS_HH
