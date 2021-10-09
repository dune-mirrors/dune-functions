// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_BACKENDS_ISTLCOMMUNICATION_HH
#define DUNE_FUNCTIONS_BACKENDS_ISTLCOMMUNICATION_HH

#include <memory>
#include <vector>

#include <dune/functions/parallel/uniqueborderpartition.hh>
#include <dune/functions/parallel/subpartitionprovider.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/typetree/traversal.hh>

namespace Dune {
namespace Functions {

/// Implementation class for ISTL-specific communication to be used in parallel solvers
template <class GlobalIndex, class LocalIndex>
class ISTLCommunication
{
  using Impl = Dune::OwnerOverlapCopyCommunication<GlobalIndex, LocalIndex>;

public:
  using IndexSet = typename Impl::ParallelIndexSet;
  using RemoteIndices = typename Impl::RemoteIndices;

public:
  template <class LocalBasis, class GlobalBasis>
  ISTLCommunication(const LocalBasis& localBasis, const GlobalBasis& globalBasis,
                    MPI_Comm comm, typename Dune::SolverCategory::Category cat)
    : comm_(comm)
    , cat_(cat)
  {
    update(localBasis, globalBasis);
  }

  const IndexSet& indexSet() const
  {
    assert(bool(impl_));
    return impl_->indexSet();
  }

  const RemoteIndices& remoteIndices() const
  {
    assert(bool(impl_));
    return impl_->remoteIndices();
  }

  typename Dune::SolverCategory::Category category() const
  {
    return cat_;
  }

  const Impl& get() const
  {
    assert(bool(impl_));
    return *impl_;
  }

  template <class LocalBasis, class GlobalBasis>
  void update(const LocalBasis& localBasis, const GlobalBasis& globalBasis);

private:
  template <class LocalBasis, class GlobalBasis>
  void updateIndexSet(const LocalBasis& localBasis, const GlobalBasis& globalBasis,
                      IndexSet& indexSet);

private:
  MPI_Comm comm_;
  typename Dune::SolverCategory::Category cat_;
  std::unique_ptr<Impl> impl_ = nullptr;
};


} // end namespace Functions
} // end namespace Dune


// ---------- implementation details ----------


template <class GlobalIndex, class LocalIndex>
  template <class LocalBasis, class GlobalBasis>
void Dune::Functions::ISTLCommunication<GlobalIndex, LocalIndex>
  :: update(const LocalBasis& localBasis, const GlobalBasis& globalBasis)
{
  impl_ = std::make_unique<Impl>(comm_, cat_);

  auto& pis = impl_->indexSet();
  updateIndexSet(localBasis, globalBasis, pis);

  auto& ris = impl_->remoteIndices();
  ris.setIndexSets(pis, pis, impl_->communicator());
  ris.template rebuild<true>();
}


template <class GlobalIndex, class LocalIndex>
  template <class LocalBasis, class GlobalBasis>
void Dune::Functions::ISTLCommunication<GlobalIndex, LocalIndex>
  ::updateIndexSet(const LocalBasis& localBasis, const GlobalBasis& globalBasis, IndexSet& indexSet)
{
  using GV = typename LocalBasis::GridView;
  using GI = typename IndexSet::GlobalIndex;
  using LI = typename IndexSet::LocalIndex;
  using Attribute = typename LI::Attribute;

  auto gv = localBasis.gridView();
  auto localMap = localBasis.localView();
  auto globalMap = globalBasis.localView();

  // make disjoint partition of border entities
  auto const& idSet = gv.grid().globalIdSet();
  UniqueBorderPartition borderEntities{gv.comm().rank(), idSet};
  gv.communicate(borderEntities,
    Dune::InterfaceType::InteriorBorder_All_Interface,
    Dune::CommunicationDirection::ForwardCommunication);
  // NOTE: this should actually be available in the GlobalGridView, already

  // insert index-pair only once
  std::vector<bool> visited(localBasis.dimension(), false);

  indexSet.beginResize();
  for (auto const& e : elements(gv))
  {
    localMap.bind(e);
    globalMap.bind(e);
    TypeTree::forEachLeafNode(localMap.tree(), [&](auto&& node, auto&&)
    {
      auto const& fe = node.finiteElement();
      for (std::size_t i = 0; i < node.size(); ++i)
      {
        const auto key = fe.localCoefficients().localKey(i);
        const auto localIndex = localMap.index(node.localIndex(i));

        if (!visited[localIndex]) {
          const GI globalIndex = globalMap.index(node.localIndex(i));

          using PType = Dune::PartitionType;
          PType pt = Dune::Hybrid::switchCases(std::make_index_sequence<GV::dimension+1>{}, key.codim(),
            [&](auto codim) {return e.template subEntity<codim>(key.subEntity()).partitionType();},
            [&]() {          return Dune::PartitionType{}; });

          switch (pt)
          {
          case PType::InteriorEntity:
            // interior entities are always owned
            indexSet.add(globalIndex, LI(localIndex, Attribute::owner, true));
            break;
          case PType::BorderEntity:
            // for border entities we need to communicate the owner status, using the
            // UniqueBorderPartition utility
            if (borderEntities.contains(idSet.subId(e, key.codim(), key.subEntity())))
              indexSet.add(globalIndex, LI(localIndex, Attribute::owner, true));
            else
              indexSet.add(globalIndex, LI(localIndex, Attribute::overlap, true));
            break;
          case PType::OverlapEntity:
            // overlap entities are never owned
            indexSet.add(globalIndex, LI(localIndex, Attribute::overlap, true));
            break;
          case PType::FrontEntity:
          case PType::GhostEntity:
            // front and ghost entities have the special status "copy"
            indexSet.add(globalIndex, LI(localIndex, Attribute::copy, true));
            break;
          default:
            assert(false && "Unknown partition type.");
            std::abort();
          }

          visited[localIndex] = true;
        }
      }
    });
    globalMap.unbind();
    localMap.unbind();
  }
  indexSet.endResize();

  // test that all indices are inserted into the indexset
  assert(indexSet.size() == localBasis.dimension());
}

#endif // DUNE_FUNCTIONS_BACKENDS_ISTLCOMMUNICATION_HH
