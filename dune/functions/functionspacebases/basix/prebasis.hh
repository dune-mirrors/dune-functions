// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASIX_PREBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASIX_PREBASIS_HH

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <map>
#include <vector>
#include <type_traits>

#include <basix/cell.h>

#include <dune/common/copyableoptional.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/basix/finiteelement.hh>
#include <dune/functions/functionspacebases/basix/utility.hh>

#include <dune/geometry/type.hh>

namespace Dune::Functions {

/**
  * \brief A pre-basis based on the basix framework
  *
  * \ingroup FunctionSpaceBasesImplementations
  *
  * \tparam GV          The grid view that the FE basis is defined on.
  * \tparam rangeClass  The class of the basis function range, i.e., RangeClass:scalar, vector or matrix.
  * \tparam Factory     A basix-finite-element factory that gets a basix cell-type as input.
  */
template<class GV, RangeClass rangeClass, class Factory>
class BasixPreBasis
    : public LeafPreBasisMapperMixin<GV>
{
  using Base = LeafPreBasisMapperMixin<GV>;

  static_assert(std::is_invocable_v<Factory, ::basix::cell::type>);

  // Define a mapper layout by counting the entity-dofs per (sub-)entity type.
  template <class Types>
  static auto makeLayout (const Factory& factory, const Types& types)
  {
    std::array<std::size_t, Basix::number_of_cell_types> sizes{};
    for (const GeometryType& type : types)
    {
      ::basix::cell::type cell_type = Basix::cellType(type);
      auto b = factory(cell_type); // TODO: we construct the basix object multiple times. Better: use a cache.

      auto& entity_dofs = b.entity_dofs();
      auto subentity_types =  ::basix::cell::subentity_types(cell_type);

      for (std::size_t d = 0; d < entity_dofs.size(); ++d) {
        for (std::size_t s = 0; s < entity_dofs[d].size(); ++s) {
          int i = d+1 < entity_dofs.size() ? (int)(subentity_types[d][s]) : (int)(cell_type);
          assert(i >= 0 && i < int(sizes.size()));
          assert(sizes[i] == 0 || sizes[i] == entity_dofs[d][s].size());
          sizes[i] = entity_dofs[d][s].size();
        }
      }
    }
    return [sizes](GeometryType gt, int) -> std::size_t { return sizes[(int)(Basix::cellType(gt))]; };
  }


public:
  using GridView = GV;
  class Node;
  using size_type = typename Base::size_type;

public:
  /** \brief Construct the pre-basis and index-mapping based on a mcmgmapper. */
  BasixPreBasis (const GV& gridView, const Factory& factory)
    : Base(gridView, makeLayout(factory,gridView.indexSet().types(0)))
    , factory_{factory}
  {
    computeCellInfos(gridView);
  }

  /** \brief Create tree node and the node finite-element. */
  Node makeNode () const
  {
    return Node(*factory_, Base::gridView().indexSet(), cellInfos_);
  }

  /** \brief Update the stored GridView and cell infos. */
  void update (const GridView& gv)
  {
    Base::update(gv);
    computeCellInfos(gv);
  }

  /** \brief Fill cache with global indices of DOFs associated to the given bound node. */
  template<class Node, class It>
  It indices(const Node& node, It it) const
  {
    auto globalIndexRange = subIndexRange(Base::mapper_, node.element(), node.finiteElement().localCoefficients());
    if (node.finiteElement().basix().dof_transformations_are_identity()) {
      for(const auto& globalIndex : globalIndexRange)
      {
        *it = {{ (size_type)globalIndex }};
        ++it;
      }
    } else {
      // global indices must be permuted to match the neighboring elements
      std::vector<std::int32_t> globalIndices(globalIndexRange.begin(), globalIndexRange.end());
      if (node.finiteElement().basix().dof_transformations_are_permutations() && node.finiteElement().cellInfo())
        node.finiteElement().basix().permute(globalIndices, node.finiteElement().cellInfo());

      for(const auto& globalIndex : globalIndices)
      {
        *it = {{ (size_type)globalIndex }};
        ++it;
      }
    }
    return it;
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

    for (int j = 0; j < r.size(i,c,Element::dimension); ++j) {
      ii[j] = r.subEntity(i,c,j,Element::dimension);
      jj[j] = indexSet.subIndex(e,ii[j],Element::dimension);
    }

    // The face must be reflected if the number of flipped edges is odd
    int flipOrientation =(((ii[0] < ii[1]) != (jj[0] < jj[1]))
                        + ((ii[0] < ii[2]) != (jj[0] < jj[2]))
                        + ((ii[1] < ii[2]) != (jj[1] < jj[2]))) % 2;

    // Compute how many rotations after the reflection are necessary
    int iimin = std::distance(ii.begin(),std::min_element(ii.begin(), ii.end()));
    int jjmin = std::distance(jj.begin(),std::min_element(jj.begin(), jj.end()));
    int rotations = std::abs(iimin - jjmin);

    return flipOrientation | rotations<<1;
  }

  // Fill the internal `cellInfos_` vector with an encoded edge/face permutation
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
          int j = Basix::entityIndex(Basix::cellType(e.type()),GridView::dimension-1,i);
          if constexpr(GridView::dimension == 2)
            cellInfo |= computeEdgeInfo(indexSet,e,r,j,1)<<(pos++);
          else if constexpr(GridView::dimension == 3) {
            cellInfo |= computeFaceInfo(indexSet,e,r,j,1)<<pos;
            pos += 3;
          }
        }
      }

      if constexpr(GridView::dimension > 2) {
        for (int i = 0; i < r.size(2); ++i)
        {
          int j = Basix::entityIndex(Basix::cellType(e.type()),GridView::dimension-2,i);
          if constexpr(GridView::dimension == 3)
            cellInfo |= computeEdgeInfo(indexSet,e,r,j,2)<<(pos++);
        }
      }

      cellInfos_[indexSet.index(e)] = cellInfo;
    }
  }

protected:
  CopyableOptional<Factory> factory_; // TODO: maybe store the fe cache directly here in the prebasis?
  std::vector<std::uint32_t> cellInfos_;
};


template <class GV, RangeClass rangeClass, class Factory>
class BasixPreBasis<GV,rangeClass,Factory>::Node
  : public LeafBasisNode
{
  static constexpr int dim = GV::dimension;

  class FiniteElementCache;

public:
  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using IndexSet = typename GV::IndexSet;
  using FiniteElement = typename FiniteElementCache::FiniteElement;

public:
  /** \brief Constructor; stores a pointer to the passed local finite-element `fe`. */
  explicit Node (const Factory& factory, const IndexSet& indexSet, const std::vector<std::uint32_t>& cellInfos)
    : cache_(factory, indexSet.types(0))
    , indexSet_(&indexSet)
    , cellInfos_(&cellInfos)
  {}

  /** \brief Return current element; might raise an error if unbound */
  const Element& element () const
  {
    assert(!!element_);
    return *element_;
  }

  /**
   * \brief Return the FiniteElement for the element we are bound to; might raise an error if unbound.
   *
   * The FiniteElement implements the corresponding interfaces of the dune-localfunctions module.
   */
  const FiniteElement& finiteElement () const
  {
    assert(!!fe_);
    return *fe_;
  }

  /** \brief Bind to element. Stores a pointer to the passed element reference. */
  void bind (const Element& e)
  {
    element_ = &e;
    fe_ = &(cache_.get(element_->type()));
    fe_->bind(element_->geometry(), (*cellInfos_)[indexSet_->index(e)]);
    this->setSize(fe_->size());
  }

protected:
  FiniteElementCache cache_;
  const IndexSet* indexSet_;
  const std::vector<std::uint32_t>* cellInfos_;

  FiniteElement* fe_ = nullptr;
  const Element* element_ = nullptr;
};


template <class GV, RangeClass rangeClass, class Factory>
class BasixPreBasis<GV,rangeClass,Factory>::Node::FiniteElementCache
{
  using Geometry = typename GV::template Codim<0>::Entity::Geometry;

public:
  using FiniteElement = BasixFiniteElement<Geometry,rangeClass,typename GV::ctype>;

public:
  /** \brief Construct the cache by instantiating all finite-elements on the GeometryType in the type-list `types`. */
  template <class Types>
  FiniteElementCache (const Factory& factory, const Types& types)
  {
    for (const GeometryType& type : types)
      map_.emplace(type, factory(Basix::cellType(type)));
  }

  /** \brief Return a reference to the finite-element on the given GeometryType `type`. */
  FiniteElement& get (const GeometryType& type)
  {
    auto it = map_.find(type);
    if (it == map_.end())
      DUNE_THROW(Dune::RangeError,"There is no FiniteElement of the requested GeometryType.");
    return it->second;
  }

private:
  std::map<GeometryType::Id, FiniteElement> map_ = {};
};

} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASIX_PREBASIS_HH
