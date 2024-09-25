// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_SUBDOMAINBASIS_HH
#define DUNE_SUBDOMAINBASIS_HH

#include <vector>

#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/indices.hh>

#include <dune/functions/common/utility.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/concepts.hh>

#include <dune/subdomains/subdomaingridview.hh>
#include <dune/subdomains/subdomainindexset.hh>


namespace Dune {
namespace Subdomains {

template <class GV, class PreBasis>
class SubdomainPreBasis
    : public PreBasis
{
  using PGridView = typename PreBasis::GridView;
  using PIndexSet = typename PGridView::IndexSet;

public:
  //! The grid view that the FE basis is defined on
  using GridView = GV;

public:
  //! Constructor for given pre-basis objects
  template <class PB>
  SubdomainPreBasis (PB&& preBasis, std::shared_ptr<PIndexSet> indexSet)
    : PreBasis(std::forward<PB>(preBasis))
    , indexSet_(std::move(indexSet))
  {}

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView () const
  {
    return PreBasis::gridView().impl();
  }

  //! Initialize the global indices
  void initializeIndices ()
  {
    indexSet_->update(gridView());
    PreBasis::initializeIndices();
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update (const GridView& gv)
  {
    indexSet_ = std::make_shared<PIndexSet>(gv.indexSet());
    PreBasis::update(PGridView{gv, indexSet_});
  }

  template <class Entity>
  void setSubdomain (Entity const& entity, int partition)
  {
    static_assert(Entity::codimension == 0);
    indexSet_->setSubdomain(entity, partition);
  }

private:
  std::shared_ptr<PIndexSet> indexSet_;
};

} // end namespace Subdomains


namespace Functions {
namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can build a SubdomainPreBasis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam numPartitions   Number of partitions
 * \param preBasisFactory  Pre-basis factory to wrap
 */
template <class PreBasisFactory>
auto subdomains(PreBasisFactory const& preBasisFactory, std::vector<int> partitions = {})
{
  return [preBasisFactory, partitions](auto const& gridView)
  {
    using GridView = std::decay_t<decltype(gridView)>;
    using PIndexSet = Subdomains::SubdomainIndexSet<typename GridView::IndexSet>;
    auto pIndexSet = std::make_shared<PIndexSet>(gridView.indexSet(), partitions);

    using PGridView = Subdomains::SubdomainGridView<GridView>;
    auto preBasis = preBasisFactory(PGridView{gridView, pIndexSet});

    using PreBasis = Subdomains::SubdomainPreBasis<GridView, decltype(preBasis)>;
    return PreBasis{std::move(preBasis), pIndexSet};
  };
}

} // end namespace BasisFactory
} // end namespace Functions
} // end namespace Dune


#endif // DUNE_SUBDOMAINBASIS_HH