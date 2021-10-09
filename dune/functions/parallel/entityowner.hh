// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_PARALLEL_ENTITYOWNER_HH
#define DUNE_FUNCTIONS_PARALLEL_ENTITYOWNER_HH

#include <algorithm>
#include <cassert>
#include <vector>

#include <dune/common/hybridutilities.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/mcmgmapper.hh>

namespace Dune {
namespace Functions {

/// Calculate unique owner rank for all entities (of a given codim) in a given GridView
template<class GridView>
class EntityOwner
{
private:
  /// A DataHandle class to calculate the minimum of a vector which is accompanied by an index set
  /**
   * \tparam IS  An indexSet of the GridView
   * \tparam V   The vector type to compute the elementwise minimum
   **/
  template<class IS, class Vec>
  class MinimumExchange
      : public Dune::CommDataHandleIF<MinimumExchange<IS,Vec>,typename Vec::value_type>
  {
  public:
    /** \brief export type of data for message buffer */
    using ValueType = typename Vec::value_type;

    /** \brief constructor. If `codim < 0` all codimensions are communicated */
    MinimumExchange (const IS& indexset, Vec& vec, int codim = -1)
      : indexset_(indexset)
      , vec_(vec)
      , codim_(codim)
    {}

    /** \brief returns true if data for this codim should be communicated */
    bool contains (int /*dim*/, int codim) const
    {
      return codim_ < 0 || codim == codim_;
    }

    /** \brief returns true if size per entity of given dim and codim is a constant */
    bool fixedSize (int /*dim*/, int /*codim*/) const
    {
      return true ;
    }

    /** \brief number of values to send */
    template <class Entity>
    std::size_t size (const Entity& e) const
    {
      return 1 ;
    }

    /** \brief pack data from user to message buffer */
    template <class MessageBuffer, class Entity>
    void gather (MessageBuffer& buff, const Entity& e) const
    {
      buff.write(vec_[indexset_.index(e)]);
    }

    /** \brief Unpack data from message buffer to user
     *
     * \param n The number of objects sent by the sender
     */
    template <class MessageBuffer, class Entity>
    void scatter (MessageBuffer& buff, const Entity& e, [[maybe_unused]] std::size_t n)
    {
      assert(n == 1);
      ValueType x;
      buff.read(x);
      if (x >= 0) // other is -1 means, he does not want it
        vec_[indexset_.index(e)] = std::min(x, vec_[indexset_.index(e)]);
    }

  private:
    const IS& indexset_;
    Vec& vec_;
    int codim_;
  };

  static MCMGLayout layout (int codim)
  {
    return [codim](GeometryType gt, int dim) {
      return codim < 0 || dim - int(gt.dim()) == codim;
    };
  }

public:

  using IndexSet = MultipleCodimMultipleGeomTypeMapper<GridView>;

public:
  /** \brief Constructor needs to know the grid function space */
  EntityOwner (const GridView& gridView, int codim = -1)
    : indexSet_{gridView, layout(codim)}
    , assignment_(indexSet_.size())
  {
    // assign own rank to entities that I might have
    for (auto const& e : elements(gridView))
    {
      Dune::Hybrid::forEach(Dune::StaticIntegralRange<int,GridView::dimension+1>{}, [&](auto cd)
      {
        if (codim >= 0 && codim != cd)
          return;

        for (unsigned int i = 0; i < e.subEntities(cd); ++i)
        {
          PartitionType subPartitionType = e.template subEntity<cd>(i).partitionType();

          assignment_[indexSet_.subIndex(e,i,cd)]
            = (subPartitionType == Dune::InteriorEntity or subPartitionType == Dune::BorderEntity)
            ? gridView.comm().rank()  // set to own rank
            : - 1;                    // it is a ghost entity, I will not possibly own it.
        }
      });
    }

    // exchange entity rank through communication
    MinimumExchange dh{indexSet_, assignment_, codim};
    gridView.communicate(dh,
      Dune::InterfaceType::InteriorBorder_All_Interface,
      Dune::CommunicationDirection::ForwardCommunication);
  }

  /** \brief Which rank is the i-th entity assigned to? */
  template <class Entity>
  int owner (const Entity& entity)
  {
    return assignment_[indexSet_.index(entity)];
  }

  template <class Entity>
  int owner (const Entity& entity, int i, unsigned int codim)
  {
    return assignment_[indexSet_.subIndex(entity,i,codim)];
  }

private:
  IndexSet indexSet_;
  std::vector<int> assignment_;
};

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_PARALLEL_ENTITYOWNER_HH