// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_RESTRICTEDINDEXSET_HH
#define DUNE_FUNCTIONS_RESTRICTEDINDEXSET_HH

#include <vector>
#include <functional>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/typeindex.hh>

namespace Dune
{

  template<typename GridView>
  class RestrictedIndexSet
  {
    using HostGridView = GridView;
    using HostIndexSet = typename GridView::IndexSet;
    using HostIndexType = typename HostIndexSet::IndexType;
  public:
    /** \brief Export the type of the entity used as parameter in the index(...) method */
    template <int cc>
    using Codim = HostIndexSet :: template Codim<cc>;

    /** \brief The type used for the indices */
    using IndexType = int;

    /** \brief iterator range for geometry types in domain */
    using Types = typename HostIndexSet::Types;

    /** \brief dimension of the grid (maximum allowed codimension) */
    static const int dimension = HostIndexSet::dimension;

    RestrictedIndexSet(
      const HostGridView& hostGridView,
      const MultipleCodimMultipleGeomTypeMapper<HostGridView>& mapper,
      const std::vector<int>& indices,
      const std::map<GeometryType,int>& sizes
      ) :
      _hostIndexSet(hostGridView.indexSet()),
      _entityMapper(mapper),
      _indices(indices),
      _sizes(sizes)
    {}


    //===========================================================
    /** @name Index access from entity
     */
    //@{
    //===========================================================

    /** @brief Map entity to index. The result of calling this method with an entity that is not
            in the index set is undefined.

            \param e Reference to codim cc entity, where cc is the template parameter of the function.
            \return An index in the range 0 ... Max number of entities in set - 1.
     */
    template<int cc>
    IndexType index (const typename HostIndexSet::template Codim<cc>::Entity& e) const
    {
      return index(e);
    }

    /** @brief Map entity to index. Easier to use than the above because codimension template
            parameter need not be supplied explicitly.
            The result of calling this method with an entity that is not
            in the index set is undefined.

            \param e Reference to codim cc entity. Since
           entity knows its codimension, automatic extraction is possible.
            \return An index in the range 0 ... Max number of entities in set - 1.
     */
    template<class Entity>
    IndexType index (const Entity& e) const
    {
      auto j = _entityMapper.index(e);
      return _indices[j];
    }

    /** \brief Map a subentity to an index.
     *
     *  The result of calling this method with an entity that is not in the
     *  index set is undefined.
     *
     *  \note This method exists for convenience only.
     *        It extracts the codimension from the type of the entity, which can
     *        be guessed by the compiler.
     *
     *  \tparam  Entity  type of entity (must be GridImp::Codim< cc >::Entity
     *                   for some cc)
     *
     *  \param[in]  e      reference to entity
     *  \param[in]  i      number subentity of e within the codimension
     *  \param[in]  codim  codimension of the subentity we're interested in
     *
     *  \note The parameter <tt>codim</tt> denotes the codimension with respect
     *        to the grid, i.e., it must satisfy cc <= codim <= dimension.
     *
     *  \return An index in the range 0 ... Max number of entities in set - 1.
     */
    template< class Entity >
    IndexType subIndex ( const Entity &e, int i, unsigned int codim ) const
    {
      auto j = _entityMapper.subIndex(e,i,codim);
      return _indices[j];
    }
    //@}


    //===========================================================
    /** @name Access to entity set
     */
    //@{
    //===========================================================

    /**
     * \brief obtain all geometry types of entities in domain
     *
     * This method returns an iterator range (something that behaves like
     * Dune::IteratorRange) visiting all geometry types of codimension codim
     * in the domain of the index map exactly once.
     * The iterator must implement the concept of a forward iterator (in the
     * sense of the STL).
     * The elements in the iterator range are required to be of type
     * Dune::GeometryType.
     *
     * \param[in]  codim  a valid codimension
     *
     * \return iterator range over Const reference to a vector of geometry types.
     */
    Types types (int codim) const
    {
      return _hostIndexSet.types( codim );
    }

    /** @brief Return total number of entities of given geometry type in entity set \f$E\f$.

       \param[in] type A valid geometry type.
       \return    number of entities (type is auto determined by the
                  implementation. std::size_t is the expected return type).
     */
    auto size (GeometryType type) const
    {
#warning TODO switch to GeometryTypeIndex instead of a map, this avoids problems with uninitialized values
      /*
        we initialize all types reported be the grid, but the lagrange
        basis does not care about these types and also queries other
        types...

        If we switch to GeometryTypeIndex, we
         - avoid expensive std::map lookups
         - can simply initialize *every* type that is know to the GeometryTypeIndex utility
       */
      try {
        return _sizes.at(type);
      }
      catch(...) {
        return 0;
      }
    }

    /** @brief Return total number of entities of given codim in the entity set \f$E\f$. This
            is simply a sum over all geometry types.

       \param[in] codim A valid codimension
       \return    number of entities (type is auto determined by the
                  implementation. std::size_t is the expected return type).
     */
    auto size (int codim) const
    {
      std::size_t sz = 0;
      // for (auto && gt : types(codim))
      //   sz += size(gt);
      for (auto & [gt,s] : _sizes)
        if (gt.dim() == (dimension - codim))
          sz += s;
      return sz;
    }

    /** @brief Return true if the given entity is contained in \f$E\f$.
     *
     * \note If the input element e is not an element of the grid, then
     *       the result of contains() is undefined.
     */
    template<class Entity>
    bool contains (const Entity& e) const
    {
      return index(e) != -1;
    }

  private:
    const HostIndexSet& _hostIndexSet;
    const MultipleCodimMultipleGeomTypeMapper<HostGridView>& _entityMapper;
    const std::vector<int>& _indices;
    const std::map<GeometryType, int>& _sizes;
  };

} // namespace Dune

#endif // DUNE_FUNCTIONS_RESTRICTEDINDEXSET_HH
