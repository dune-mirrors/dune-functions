// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_MULTIDOMAIN_DOMININFO_HH
#define DUNE_FUNCTIONS_MULTIDOMAIN_DOMININFO_HH

namespace Dune::Functions::MultiDomain {

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

  template<typename GridView>
  void update(const GridView &) {
    std::cout << "WARNING ! not updating the PartitionedDomainInfo" << std::endl;
  }
};

auto createPartitionedDomainInfo(std::vector<int> && partitions,
  std::vector<std::vector<int>> && domains)
{
  return std::make_shared<PartitionedDomainInfo>(
    std::forward<std::vector<int>>(partitions),
    std::forward<std::vector<std::vector<int>>>(domains));
}

} // end namespace Dune::Functions::MultiDomain

#endif // DUNE_FUNCTIONS_MULTIDOMAIN_DOMININFO_HH
