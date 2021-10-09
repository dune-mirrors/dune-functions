// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_PARALLEL_SUBPARTITIONTYPEPROVIDER_HH
#define DUNE_FUNCTIONS_PARALLEL_SUBPARTITIONTYPEPROVIDER_HH

namespace Dune {
namespace Functions {

  /** \brief Helper class to provide access to subentity PartitionTypes with a run-time codimension
   *
   * This class can be removed if there is ever a method 'subPartitionType' similar to 'subIndex',
   * that takes a run-time codimension argument.
   */
  template <class Entity, int cd>
  struct SubPartitionTypeProvider
  {
    /** \brief Get PartitionType of the i-th subentity of codimension 'codim' of entity 'entity'
     */
    static PartitionType get (const Entity& entity, int codim, int i)
    {
      if (codim == cd)
        return entity.template subEntity<cd>(i).partitionType();
      else
        return SubPartitionTypeProvider<Entity,cd-1>::get(entity, codim, i);
    }
  };

  template <class Entity>
  struct SubPartitionTypeProvider<Entity,0>
  {
    static PartitionType get (const Entity& entity, int codim, int i)
    {
      return entity.template subEntity<0>(i).partitionType();
    }
  };

} // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_PARALLEL_SUBPARTITIONTYPEPROVIDER_HH
