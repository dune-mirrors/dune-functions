// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_MULTIDOMAIN_MULTIDOMAINGRIDVIEW_HH
#define DUNE_FUNCTIONS_MULTIDOMAIN_MULTIDOMAINGRIDVIEW_HH

#include <memory>
#include <type_traits>
#include <utility>

#include "restrictedgridview.hh"
#include "restrictedindexset.hh"

namespace Dune::Functions::MultiDomain {

template <class GridView>
class MultiDomainGridView : public GridView
{
public:
  using HostGridView = GridView;
  using SubdomainGridView = RestrictedGridView<GridView>;
  using SubdomainIndexSet = typename SubdomainGridView::IndexSet;

  MultiDomainGridView(const GridView & gridView, const std::shared_ptr<PartitionedDomainInfo> & domainInfo) :
    GridView(gridView),
    // _cellMapper(gridView, mcmgElementLayout()),
    // _entityMapper(gridView, [](GeometryType, int) { return true; }),
    _cellMapper((const GridView&)(*this), mcmgElementLayout()),
    // here we might optimize and setup only required indices
    _entityMapper((const GridView&)(*this), [](GeometryType, int) { return true; }),
    _domainInfo(domainInfo)
  {
    unsigned int subdomains = _domainInfo->subdomains().size();
    _sizes.resize(subdomains);
    _indices.resize(subdomains);
    _indexSets.resize(subdomains);
    _gridViews.resize(subdomains);

    updateRestrictedgridViews();

    // setup indexSets
    for (int d=0; d<subdomains; d++)
      _indexSets[d] = std::make_shared<SubdomainIndexSet>(
          *this, _entityMapper,
          _indices[d], _sizes[d]);

    // seup gridViews
    for (int d=0; d<subdomains; d++)
      _gridViews[d] = std::make_shared<SubdomainGridView>(
          *this, _indexSets[d]);
  }

  /** \brief obtain the index set
   *
   * The lifetime of the returned index set is bound to the lifetime of the
   * grid view. Keep a copy of the grid view to prevent the index set from
   * becoming a dangling reference.
   */
  const SubdomainGridView & subdomainGridView (int i) const
  {
    assert(i < _gridViews.size());
    return *_gridViews[i];
  }

  void update(const GridView & backgroundView)
  {
    (GridView &)*this = backgroundView;
    _domainInfo->update(backgroundView);
    _cellMapper.update(backgroundView);
    _entityMapper.update(backgroundView);
    updateRestrictedgridViews();
  }

private:
  void updateRestrictedgridViews()
  {
    constexpr int dim = GridView::dimension;
    std::size_t subdomains = _domainInfo->subdomains().size();

    // create a temporary map from partition index to subdomain index
    std::map<int,std::vector<int>> parts;
    for (int d = 0; d<subdomains; d++)
      for (int p : _domainInfo->subdomain(d))
        parts[p].push_back(d);

    // clear and initialize indices and sizes
    for (int dom=0; dom<subdomains; dom++)
    {
      _indices[dom].resize(_entityMapper.size());
      std::fill(_indices[dom].begin(), _indices[dom].end(), -1);
      _sizes[dom].clear();
      for (unsigned int d = 0; d <= dim; d++)
      {
        for (int t=0; t<LocalGeometryTypeIndex::size(d); t++)
        {
          auto gt = LocalGeometryTypeIndex::type(d,t);
          _sizes[dom][gt] = 0;
        }
      }
    }

    // loop through mesh, check sub-domain and update indices
    for (auto && e : elements(*this))
    {
      auto i = _cellMapper.index(e);
      auto refcell = referenceElement<double,dim>(e.type());
      auto p = _domainInfo->partition(i);
      for (int subdomain : parts[p])
      {
        // every subentity is part of the subdomain
        for (int codim = 0; codim<=dim; codim++)
          // n'th subentity of codim
          for (int n = 0; n<refcell.size(codim); n++)
          {
            // use the next unused index within this geometry type
            auto gt = refcell.type(n, codim);
            auto j = _entityMapper.subIndex(e,n,codim);
            // set new index, if not already set
            if (_indices[subdomain][j] == -1)
              _indices[subdomain][j] = _sizes[subdomain][gt]++;
          }
      }
    }
  }

  MultipleCodimMultipleGeomTypeMapper<GridView> _cellMapper;
  MultipleCodimMultipleGeomTypeMapper<GridView> _entityMapper;
  std::shared_ptr<PartitionedDomainInfo> _domainInfo;
  std::vector<std::map<GeometryType,int>> _sizes;
  std::vector<std::vector<int>> _indices;
  std::vector<std::shared_ptr<SubdomainIndexSet>> _indexSets;
  std::vector<std::shared_ptr<SubdomainGridView>> _gridViews;
};

}

#endif // DUNE_FUNCTIONS_MULTIDOMAIN_MULTIDOMAINGRIDVIEW_HH
