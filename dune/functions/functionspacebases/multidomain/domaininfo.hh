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
