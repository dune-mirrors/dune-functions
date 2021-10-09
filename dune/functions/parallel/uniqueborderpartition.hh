// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_PARALLEL_UNIQUEBORDERPARTITION_HH
#define DUNE_FUNCTIONS_PARALLEL_UNIQUEBORDERPARTITION_HH

#include <cassert>
#include <set>

#include <dune/geometry/dimension.hh>
#include <dune/grid/common/datahandleif.hh>

namespace Dune {
namespace Functions {

/// \brief Determine for each border entity which processor owns it
/**
 * All entities must be uniquely owned by exactly one processor, but they can
 * exist on multiple processors. For interior, overlap, and ghost entities the
 * assignment is trivial: interior: owner, otherwise not owner.
 *
 * For border entities (codim != 0) the ownership is not known a priori and must
 * be communicated. Here we assign the entity to the processor with the lowest rank.
 **/
template <class IdSet>
class UniqueBorderPartition
    : public Dune::CommDataHandleIF<UniqueBorderPartition<IdSet>, int>
{
  using IdType = typename IdSet::IdType;

public:
  using EntitySet = std::set<IdType>;

public:
  /// \brief Construct a UniqueBorderPartition DataHandle to be used in a GridView
  /// communicator.
  /**
   * \param rank            The own processor rank
   * \param idSet           The id set of entity ids to store in borderEntities,
   *                        typically the grid globalIdSet.
   *
   * NOTE: Since idSet is stored by reference it must not go out of scope
   * until all calls to \ref gather and \ref scatter are finished.
   **/
  UniqueBorderPartition(int rank, IdSet const& idSet)
    : myrank_(rank)
    , idSet_(idSet)
  {}

  /** \brief Communicate all entities */
  bool contains(int /*dim*/, int /*codim*/) const { return true; }

  /** \brief communicate exactly one integer, the rank */
  bool fixedSize(int /*dim*/, int /*codim*/) const { return true; }

  /** \brief Always contains one int, the rank */
  template <class Entity>
  std::size_t size(Entity const& e) const { return 1; }


  template <class MessageBuffer, class Entity>
  void gather(MessageBuffer& buff, Entity const& e) const
  {
    buff.write(myrank_);
  }

  template <class MessageBuffer, class Entity>
  void scatter(MessageBuffer& buff, Entity const& e, std::size_t n)
  {
    scatterImpl(buff, e, n, Codim<Entity::codimension>{});
  }

  /** \brief Returns whether id is owned by this rank */
  bool contains(IdType const& id) const
  {
    return notOwner_.count(id) == 0;
  }

private:

  template <class MessageBuffer, class Entity, int cd>
  void scatterImpl(MessageBuffer& buff, Entity const& e, [[maybe_unused]] std::size_t n, Codim<cd>)
  {
    assert(n == 1);

    int rank = 0;
    buff.read(rank);

    // insert only border entities that are owned by other processors, i.e. rank > myrank
    // Those entities are not owned by this rank.
    if (rank > myrank_)
      notOwner_.insert(idSet_.id(e));
  }

  template <class MessageBuffer, class Entity>
  void scatterImpl(MessageBuffer& buff, Entity const& e, [[maybe_unused]] std::size_t n, Codim<0>)
  {
    assert(n == 1);

    int rank = 0;
    buff.read(rank);

    // insert only border entities that are owned by other processors, i.e. rank > myrank
    // Those entities are not owned by this rank.
    if (rank > myrank_) {
      for (int codim = 1; codim <= IdSet::dimension; ++codim) {
        for (unsigned int i = 0; i < e.subEntities(codim); ++i) {
          notOwner_.insert(idSet_.subId(e, i, codim));
        }
      }
    }
  }

private:
  int myrank_;
  EntitySet notOwner_;
  IdSet const& idSet_;
};

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_PARALLEL_UNIQUEBORDERPARTITION_HH
