// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GLOBALGRIDVIEW_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GLOBALGRIDVIEW_HH

#include <array>

#include <dune/common/indices.hh>
#include <dune/common/typeutilities.hh>
#include <dune/grid/utility/globalindexset.hh>

namespace Dune {
namespace Functions {

/** \brief A GridView wrapper that provides the GlobalIndexSet
  *
  * \tparam HGV  The GridView to be wrapped.
  */
template<class HGV>
class GlobalGridView
    : public HGV
{
  using HostGridView = HGV;

public:
  /** \brief A GlobalIndexSet for each codimension */
  class IndexSet
  {
    using HostIndexSet = typename HostGridView::IndexSet;
    using CodimIndexSet = GlobalIndexSet<HostGridView>;

  public:
    /** \brief Export the type of the entity used as parameter in the index(...) method */
    template<int cc>
    struct Codim
    {
      using Entity = typename HostGridView::template Codim<cc>::Entity;
    };

    /** \brief The type used for the indices */
    using IndexType = typename CodimIndexSet::Index;

    /** \brief iterator range for geometry types in domain */
    using Types = typename HostIndexSet::Types;

    /** \brief maximum allowed codimension */
    static const int dimension = HostIndexSet::dimension;

  public:
    IndexSet (const HostGridView& hostGridView)
      : hostIndexSet_{&hostGridView.indexSet()}
      , codimIndexSets_{ Dune::unpackIntegerSequence(
          [&](auto... cc) { return std::array{CodimIndexSet{hostGridView, cc}...}; },
          std::make_index_sequence<dimension+1>{})
        }
    {}

    /** \brief Geometry types accessible by the indexSet */
    Types types (int codim) const
    {
      return hostIndexSet_->types(codim);
    }

    /** \brief Return the global index of a given entity of codim cc */
    template<int cc>
    IndexType index (const typename Codim<cc>::Entity& entity) const
    {
      return codimIndexSets_[cc].index(entity);
    }

    /** \brief Return the global index of a given entity */
    template<class Entity>
    IndexType index (const Entity& entity) const
    {
      return codimIndexSets_[Entity::codimension].index(entity);
    }

    /** \brief Return the global index of a subentity of a given entity
     *
     * \param i Number of the requested subentity among all subentities of the given codimension
     * \param codim Codimension of the requested subentity
     */
    template<int cc>
    IndexType subIndex (const typename Codim<cc>::Entity& entity,
                        int i, unsigned int codim) const
    {
      return codimIndexSets_[cc].subIndex(entity, i, codim);
    }

    /** \brief Return the global index of a subentity of a given entity
     *
     * \param i Number of the requested subentity among all subentities of the given codimension
     * \param codim Codimension of the requested subentity
     */
    template<class Entity>
    IndexType subIndex (const Entity& entity, unsigned int i, unsigned int codim) const
    {
      return codimIndexSets_[Entity::codimension].subIndex(entity, i, codim);
    }

    /** \brief Return the total number of entities over all processes that we have indices for */
    auto size (int codim) const
    {
      return codimIndexSets_[codim].size(codim);
    }

    /** \brief Return the total number of entities over all processes that we have indices for */
    auto size (GeometryType type) const
    {
      return size(dimension - type.dim());
    }

    /** \brief Return true if the given entity is contained in the indexSet */
    template<class Entity>
    bool contains (const Entity& e) const
    {
      DUNE_THROW(Dune::NotImplemented, "GlobalIndexSet::contains not implemented.");
      return true;
    }

  private:
    /** \brief This indexSet of the HostGridView */
    const HostIndexSet* hostIndexSet_;

    /** \brief For each codimension a global indexSet */
    std::array<CodimIndexSet, dimension+1> codimIndexSets_;
  };

public:
  GlobalGridView (const HostGridView& hostGridView)
    : HostGridView{hostGridView}
    , indexSet_{static_cast<const HostGridView&>(*this)}
  {
    static_assert(Dune::Capabilities::hasSingleGeometryType<typename HostGridView::Grid>::v,
      "Only grids with single GeometryType are supported.");
  }

  const IndexSet& indexSet () const
  {
    return indexSet_;
  }

  /** \brief obtain number of entities in a given codimension */
  int size (int codim) const
  {
    return indexSet_.size(codim);
  }

  /** \brief obtain number of entities with a given geometry type */
  int size (GeometryType type) const
  {
    return indexSet_.size(type);
  }

  /** \brief Return true if the given entity is contained in this grid view */
  template<class Entity>
  bool contains (const Entity& entity) const
  {
    return indexSet_.contains(entity);
  }

private:
  IndexSet indexSet_;
};


} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_GLOBALGRIDVIEW_HH
