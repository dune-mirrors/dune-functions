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


#include <dune/common/copyableoptional.hh>
#include <dune/localfunctions/basix/globalbasix.hh>
#include <dune/localfunctions/basix/utility.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>

#include <dune/geometry/type.hh>

#include <basix/e-lagrange.h>

namespace Dune::Functions {

/**
  * \brief A pre-basis based on the basix framework
  *
  * \ingroup FunctionSpaceBasesImplementations
  *
  * \tparam GV        The grid view that the FE basis is defined on.
  * \tparam rangeClass  The class of the basis function range, i.e., RangeClass:scalar, vector or matrix.
  */
template<class GV, RangeClass rangeClass, class Factory>
class BasixPreBasis
    : public LeafPreBasisMapperMixin<GV>
{
  using Base = LeafPreBasisMapperMixin<GV>;

  template <class Types>
  static auto makeLayout (const Factory& factory, const Types& types)
  {
    std::array<std::size_t, Dune::Impl::number_of_cell_types> sizes{};
    for (const GeometryType& type : types)
    {
      ::basix::cell::type cell_type = Dune::Impl::cellType(type);
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
    return [sizes](GeometryType gt, int) -> std::size_t { return sizes[(int)(Dune::Impl::cellType(gt))]; };
  }


public:
  using GridView = GV;
  class Node;
  using size_type = typename Base::size_type;

public:
  BasixPreBasis (const GV& gridView, const Factory& factory)
    : Base(gridView, makeLayout(factory,gridView.indexSet().types(0)))
    , factory_{factory}
  {
    computeCellInfos(gridView);
  }

  using Base::gridView;

  //! Create tree node
  Node makeNode () const
  {
    return Node(*factory_, gridView().indexSet(), cellInfos_);
  }

  //! Update the stored GridView and cell infos
  void update (const GridView& gv)
  {
    Base::update(gv);
    computeCellInfos(gv);
  }

  //! Fill cache with global indices of DOFs associated to the given bound node
  template<class Node, class It>
  It indices(const Node& node, It it) const
  {
    auto globalIndexRange = subIndexRange(Base::mapper_, node.element(), node.finiteElement().localCoefficients());
    std::vector<std::int32_t> globalIndices(globalIndexRange.begin(), globalIndexRange.end());

    if (node.finiteElement().basix().dof_transformations_are_permutations() && node.finiteElement().cellInfo()) {
      node.finiteElement().basix().permute(globalIndices, node.finiteElement().cellInfo());
    }

    for(const auto& globalIndex : globalIndices)
    {
      *it = {{ (size_type)globalIndex }};
      ++it;
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
    int flipOrientation =(((ii[0] < ii[1]) != (jj[0] < jj[1]))
                        + ((ii[0] < ii[2]) != (jj[0] < jj[2]))
                        + ((ii[1] < ii[2]) != (jj[1] < jj[2]))) % 2;

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
          int j = Dune::Impl::entityIndex(Dune::Impl::cellType(e.type()),GridView::dimension-1,i);
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
          int j = Dune::Impl::entityIndex(Dune::Impl::cellType(e.type()),GridView::dimension-2,i);
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

  //! Constructor; stores a pointer to the passed local finite-element `fe`.
  explicit Node (const Factory& factory, const IndexSet& indexSet, const std::vector<std::uint32_t>& cellInfos)
    : cache_(factory, indexSet.types(0))
    , indexSet_(&indexSet)
    , cellInfos_(&cellInfos)
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
    return *fe_;
  }

  //! Bind to element. Stores a pointer to the passed element reference.
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
  template <class Types>
  FiniteElementCache (const Factory& factory, const Types& types)
  {
    for (const GeometryType& type : types)
      map_.emplace(type, factory(Dune::Impl::cellType(type)));
  }

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
  return [basix=std::move(basix)]<class GridView>(const GridView& gridView) {
    auto factory = [basix=std::move(basix)](::basix::cell::type cell_type) {
      assert(basix.cell_type() == cell_type);
      return basix;
    };
    return BasixPreBasis<GridView, rangeClass, decltype(factory)>(gridView, std::move(factory));
  };
}

template<class F = double>
auto basix_lagrange (int degree,
  ::basix::element::lagrange_variant variant = ::basix::element::lagrange_variant::equispaced)
{
  return [degree,variant]<class GridView>(const GridView& gridView) {
    auto factory = [degree,variant](::basix::cell::type cell_type) {
      return ::basix::element::create_lagrange<F>(cell_type, degree, variant, false);
    };
    return BasixPreBasis<GridView, RangeClass::scalar, decltype(factory)>(gridView, std::move(factory));
  };
}

template<class F = double>
auto basix_lagrangedg (int degree,
  ::basix::element::lagrange_variant variant = ::basix::element::lagrange_variant::equispaced)
{
  return [degree,variant]<class GridView>(const GridView& gridView) {
    auto factory = [degree,variant](::basix::cell::type cell_type) {
      return ::basix::element::create_lagrange<F>(cell_type, degree, variant, true);
    };
    return BasixPreBasis<GridView, RangeClass::scalar, decltype(factory)>(gridView, std::move(factory));
  };
}

} // end namespace BasisFactory
} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASIXBASIS_HH
