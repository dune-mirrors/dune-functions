// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_COMMON_MULTIDOMAIN_HH
#define DUNE_FUNCTIONS_COMMON_MULTIDOMAIN_HH

#include <array>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/iteratorrange.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/functions/common/subdomain.hh>

namespace Dune::Functions {

  class PartitionedDomainInfo
  {
    std::vector<int> _partitions;
    std::vector<std::vector<int>> _subdomains;
  public:
    PartitionedDomainInfo() = delete;

    PartitionedDomainInfo(std::vector<int> && partitions,
      std::vector<std::vector<int>> && domains) :
      _partitions ( std::move(partitions) ),
      _subdomains (std::move( domains) )
    {}

    const std::vector<std::vector<int>>& subdomains() const {
      return _subdomains;
    }

    const std::vector<int>& subdomain(std::size_t i) const {
      return _subdomains[i];
    }

    const std::vector<int>& partitions() const {
      return _partitions;
    }

    const int& partition(std::size_t i) const {
      return _partitions[i];
    }
  };

  auto createPartitionedDomainInfo(std::vector<int> && partitions,
    std::vector<std::vector<int>> && domains)
  {
    return std::make_shared<PartitionedDomainInfo>(
      std::forward<std::vector<int>>(partitions),
      std::forward<std::vector<std::vector<int>>>(domains));
  }

  /**
   * \brief Class representing a multi-domain of a GridView
   *
   * \ingroup Utility
   *
   * A MultiDomain is a multiset of grid elements from a given
   * underlying grid view together with their multi-entities.
   * It allows to create a `MultiDomainGridView` which implements
   * a reasonable multiset of the grid view interface defined
   * in dune-grid. In particular the `MultiDomainGridView`
   * implements an index set, a `contains()` methods, an
   * element iterator and intersection iterators.
   *
   * \tparam HGV The underlying host grid view type.
   */
  template<class GV>
  class MultiDomain
  {
  public:

    using GridView = GV;
    using Grid = typename GridView::Grid;
    using SubDomain = Dune::Functions::SubDomain<GridView>;

    //! Construct MultiDomain for underlying host grid view
    template<class Partition>
    MultiDomain(const GridView& gridView,
      const Partition & partition, // maps e -> partition number
      const std::vector<std::vector<int>> & subdomains)
      : gridView_(gridView)
    {
      subdomains_.resize(subdomains.size(), SubDomain(gridView_));

      // create a temporary map from partition index to subdomain index
      std::map<std::size_t,std::vector<int>> parts;
      for (int d = 0; d<subdomains.size(); d++)
        for (int p : subdomains[d])
          parts[p].push_back(d);

      // loop through elements and insert these according to their partition
      for(auto && e : elements(gridView))
      {
        int p = partition(e);
        for (const auto & d : parts[p])
          subdomains_[d].insertElement(e);
      }
    }

    std::size_t subDomains() const {
      return subdomains_.size();
    }

    const SubDomain & subDomain(int i) const {
      return subdomains_[i];
    }

  private:
    GridView gridView_;
    std::vector<SubDomain> subdomains_;
  };



} // namespace Dune::Functions

#endif// DUNE_FUNCTIONS_COMMON_MULTIDOMAIN_HH
