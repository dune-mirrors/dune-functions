// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_PARALLEL_GLOBALINDEXSET_HH
#define DUNE_FUNCTIONS_PARALLEL_GLOBALINDEXSET_HH

#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <map>
#include <utility>
#include <algorithm>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>

#if HAVE_MPI
  #include <dune/common/parallel/mpihelper.hh>
#endif

#include <dune/functions/parallel/entityowner.hh>

namespace Dune {
namespace Functions {

/** \brief Calculate globally unique index over all processes in a Dune grid
 */
template<class GridView>
class GlobalIndexSet
{
public:
  /** \brief The number type used for global indices  */
  using IndexType = int;

private:
  /** define data types */
  typedef typename GridView::Grid Grid;

  typedef typename GridView::Grid::GlobalIdSet GlobalIdSet;
  typedef typename GridView::Grid::GlobalIdSet::IdType IdType;
  typedef typename GridView::Traits::template Codim<0>::Iterator Iterator;

  typedef typename Grid::CollectiveCommunication CollectiveCommunication;

  using MapId2Index = std::map<IdType,IndexType>;
  using IndexMap = std::map<IndexType,IndexType>;

private:
  /* A DataHandle class to communicate the global index from the
    * owning to the non-owning entity; the class is based on the MinimumExchange
    * class in the parallelsolver.hh header file.
    */
  class IndexExchange
    : public Dune::CommDataHandleIF<IndexExchange,Index>
  {
  public:
    //! returns true if data for this codim should be communicated
    bool contains (int dim, int codim) const
    {
      return codim_ < 0 || codim == codim_;
    }

    //! returns true if size per entity of given dim and codim is a constant
    bool fixedSize (int dim, int codim) const
    {
      return true;
    }

    /** \brief How many objects of type DataType have to be sent for a given entity
     *
     * \note Only the sender side needs to know this size.
     */
    template<class Entity>
    size_t size (Entity& e) const
    {
      return 1;
    }

    /*! pack data from user to message buffer */
    template <class MessageBuffer, class Entity>
    void gather (MessageBuffer& buff, const Entity& e) const
    {
      IdType id = globalIdSet_.id(e);

      if (codim_ == 0)
        buff.write(mapid2entity_[id]);
      else
        buff.write((*mapid2entity_.find(id)).second);
    }

    /** \brief Unpack data from message buffer to user
     *
     * \param n The number of objects sent by the sender
     */
    template <class MessageBuffer, class Entity>
    void scatter (MessageBuffer& buff, const Entity& entity, std::size_t n)
    {
      assert(n == 1);

      IndexType idx;
      buff.read(idx);

      /** only if the incoming index is a valid one,
       *  i.e. if it is greater than zero, will it be
       *  inserted as the global index; it is made
       *  sure in the upper class, i.e. GlobalIndexSet,
       *  that non-owning processes use -1 to mark an entity
       *  that they do not own.
       */
      if (idx != IndexType(-1)) {
        const IdType id = globalIdSet_.id(entity);

        if (codim_ == 0)
          mapid2entity_[id] = idx;
        else
        {
          mapid2entity_.erase(id);
          mapid2entity_.insert(std::make_pair(id,x));

          const IndexType lindex = indexSet_.index(entity);
          localGlobalMap_[lindex] = x;
        }
      }
    }

    //! constructor
    IndexExchange (const GlobalIdSet& globalIdSet,
                   MapId2Index& mapid2entity,
                   const typename GridView::IndexSet& localIndexSet,
                   IndexMap& localGlobal,
                   int codim = -1)
      : globalIdSet_(globalIdSet)
      , mapid2entity_(mapid2entity)
      , indexSet_(localIndexSet)
      , localGlobalMap_(localGlobal)
      , codim_(codim)
    {}

  private:
    const GlobalIdSet& globalIdSet_;
    MapId2Index& mapid2entity_;

    const typename GridView::IndexSet& indexSet_;
    IndexMap& localGlobalMap_;
    int codim_;
  };

public:
  /** \brief Constructor for a given GridView
   *
   * This constructor calculates the complete set of global unique indices so that we can then
   *  later query the global index, by directly passing the entity in question.
   */
  GlobalIndexSet (const GridView& gridview, int codim = -1)
    : gridview_(gridview)
    , codim_(codim)
  {
    const int rank = gridview.comm().rank();
    const int size = gridview.comm().size();

    const typename GridView::IndexSet& indexSet = gridview.indexSet();

    std::unique_ptr<UniqueEntityPartition> uniqueEntityPartition;
    if (codim_!=0)
      uniqueEntityPartition = std::make_unique<UniqueEntityPartition>(gridview,codim_);

    int nLocalEntity = (codim_ == 0)
      ? std::distance(gridview.template begin<0, Dune::Interior_Partition>(),
                      gridview.template end<0, Dune::Interior_Partition>())
      : uniqueEntityPartition->numOwners(rank);

    // Compute the global, non-redundant number of entities, i.e. the number of entities in the set
    // without double, aka. redundant entities, on the interprocessor boundary via global reduce. */
    nGlobalEntity_ = gridview.comm().template sum<int>(nLocalEntity);

    /* communicate the number of locally owned entities to all other processes so that the respective offset
      * can be calculated on the respective processor; we use the Dune mpi collective communication facility
      * for this; first, we gather the number of locally owned entities on the root process and, second, we
      * broadcast the array to all processes where the respective offset can be calculated. */

    std::vector<int> offset(size);
    std::fill(offset.begin(), offset.end(), 0);

    /** Share number of locally owned entities */
    gridview_.comm().template allgather<int>(&nLocalEntity, 1, offset.data());

    int myoffset = 0;
    for (int i=1; i<rank+1; i++)
      myoffset += offset[i-1];

    /*  compute globally unique index over all processes; the idea of the algorithm is as follows: if
      *  an entity is owned by the process, it is assigned an index that is the addition of the offset
      *  specific for this process and a consecutively incremented counter; if the entity is not owned
      *  by the process, it is assigned -1, which signals that this specific entity will get its global
      *  unique index through communication afterwards;
      *
      *  thus, the calculation of the globally unique index is divided into 2 stages:
      *
      *  (1) we calculate the global index independently;
      *
      *  (2) we achieve parallel adjustment by communicating the index
      *      from the owning entity to the non-owning entity.
      *
      */

    // 1st stage of global index calculation: calculate global index for owned entities
    // initialize map that stores an entity's global index via it's globally unique id as key
    globalIndex_.clear();

    const GlobalIdSet& globalIdSet = gridview_.grid().globalIdSet();      /** retrieve globally unique Id set */

    Index globalcontrib = 0;      /** initialize contribution for the global index */

    if (codim_==0)  // This case is simpler
    {
      for (Iterator iter = gridview_.template begin<0>(); iter!=gridview_.template end<0>(); ++iter)
      {
        const IdType id = globalIdSet.id(*iter);      /** retrieve the entity's id */

        /** if the entity is owned by the process, go ahead with computing the global index */
        if (iter->partitionType() == Dune::InteriorEntity)
        {
          const IndexType gindex = myoffset + globalcontrib;    /** compute global index */

          globalIndex_[id] = gindex;                      /** insert pair (key, datum) into the map */
          globalcontrib++;                                /** increment contribution to global index */
        }

        /** if entity is not owned, insert -1 to signal not yet calculated global index */
        else
        {
          globalIndex_[id] = -1;     /** insert pair (key, datum) into the map */
        }
      }
    }
    else  // if (codim==0) else
    {
    std::vector<bool> firstTime(gridview_.size(codim_));
    std::fill(firstTime.begin(), firstTime.end(), true);

    for(Iterator iter = gridview_.template begin<0>();iter!=gridview_.template end<0>(); ++iter)
    {
      for (size_t i=0; i<iter->subEntities(codim_); i++)
      {
        IdType id=globalIdSet.subId(*iter,i,codim_);

        IndexType idx = gridview_.indexSet().subIndex(*iter,i,codim_);

        if (!firstTime[idx] )
          continue;

        firstTime[idx] = false;

        if (uniqueEntityPartition->owner(idx) == rank)  /** if the entity is owned by the process, go ahead with computing the global index */
        {
          const IndexType gindex = myoffset + globalcontrib;    /** compute global index */
          globalIndex_.insert(std::make_pair(id,gindex)); /** insert pair (key, value) into the map */

          const IndexType lindex = idx;
          localGlobalMap_[lindex] = gindex;

          globalcontrib++;                                /** increment contribution to global index */
        }
        else /** if entity is not owned, insert -1 to signal not yet calculated global index */
        {
          globalIndex_.insert(std::make_pair(id,-1));
        }
      }

    }
    }

    // 2nd stage of global index calculation: communicate global index for non-owned entities

    // Create the data handle and communicate.
    IndexExchange dataHandle(globalIdSet,globalIndex_,indexSet,localGlobalMap_,codim_);
    gridview_.communicate(dataHandle, Dune::All_All_Interface, Dune::ForwardCommunication);
  }

  /** \brief Return the global index of a given entity */
  template <class Entity>
  IndexType index(const Entity& entity) const
  {
    if (codim_==0)
    {
      /** global unique index is only applicable for inter or border type entities */
      const GlobalIdSet& globalIdSet = gridview_.grid().globalIdSet(); /** retrieve globally unique Id set */
      const IdType id = globalIdSet.id(entity);                        /** obtain the entity's id */
      const IndexType gindex = globalIndex_.find(id)->second;                /** retrieve the global index in the map with the id as key */

      return gindex;
    }
    else
      return localGlobalMap_.find(gridview_.indexSet().index(entity))->second;
  }

  /** \brief Return the global index of a subentity of a given entity
   *
   * \param i Number of the requested subentity among all subentities of the given codimension
   * \param codim Codimension of the requested subentity
   */
  template <class Entity>
  IndexType subIndex(const Entity& entity, unsigned int i, unsigned int codim) const
  {
    if (codim_==0)
    {
      /** global unique index is only applicable for inter or border type entities */
      const GlobalIdSet& globalIdSet = gridview_.grid().globalIdSet(); /** retrieve globally unique Id set */
      const IdType id = globalIdSet.subId(entity,i,codim);                        /** obtain the entity's id */
      const IndexType gindex = globalIndex_.find(id)->second;                /** retrieve the global index in the map with the id as key */

      return gindex;
    }
    else
      return localGlobalMap_.find(gridview_.indexSet().subIndex(entity,i,codim))->second;
  }

  /** \brief Return the total number of entities over all processes that we have indices for
   *
   * \param codim If this matches GlobalIndexSet codimension, the number of entities is returned.
   *              Otherwise, zero is returned.
   */
  unsigned int size(unsigned int codim) const
  {
    return (codim_==codim) ? nGlobalEntity_ : 0;
  }

protected:
  const GridView gridview_;

  /** \brief Codimension of the entities that we hold indices for */
  unsigned int codim_;

  //! Global number of entities, i.e. number of entities without redundant entities on interprocessor boundaries
  int nGlobalEntity_;

  IndexMap localGlobalMap_;

  /** \brief Stores global index of entities with entity's globally unique id as key
   */
  MapId2Index globalIndex_;
};

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_PARALLEL_GLOBALINDEXSET_HH
