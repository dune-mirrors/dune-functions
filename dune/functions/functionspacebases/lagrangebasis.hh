// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH

#include <type_traits>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>


namespace Dune {
namespace Functions {

// *****************************************************************************
// This is the reusable part of the LagrangeBasis. It contains
//
//   LagrangePreBasis
//   LagrangeNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

template<typename GV, int k, typename R=double>
class LagrangeNode;

template<typename GV, int k, typename R=double>
class LagrangePreBasis;



/**
 * \brief A pre-basis for a PQ-lagrange bases with given order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam k   The polynomial order of ansatz functions; -1 means 'order determined at run-time'
 * \tparam R   Range type used for shape function values
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * \warning For pyramid elements in 3d, the shape functions are different for
 *    run-time and compile-time order. While the former are defined using the
 *    Duffy-transformation the latter are continuous and piecewise polynomial
 *    with discontinuous gradients.
 */
template<typename GV, int k, typename R>
class LagrangePreBasis :
  public LeafPreBasisMapperMixin<GV, Impl::EdgeTwist<typename GV::IndexSet>>
{
  using Twist = Impl::EdgeTwist<typename GV::IndexSet>;
  using Mixin = LeafPreBasisMapperMixin<GV, Twist>;

  static const int dim = GV::dimension;
  static const bool useDynamicOrder = (k<0);

  static auto dofLayout(unsigned int order)
  {
    return [order](Dune::GeometryType type, int dim) -> std::size_t {
      if (order==0)
        return int(type.dim()) == dim;
      if (type.isSimplex())
        return Dune::binomial(order-1, type.dim());
      if (type.isCube())
        return Dune::power(order-1, type.dim());
      if (type.isPyramid() and (order>2))
        return (order-2)*(order-1)*(2*order-3)/6;
      if (type.isPrism() and (order>1))
        return (order-1)*(order-1)*(order-2)/2;
      return 0;
    };
  }

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Template mapping root tree path to type of created tree node
  using Node = LagrangeNode<GV, k, R>;

  //! Constructor for a given grid view object with compile-time order
  LagrangePreBasis(const GridView& gv) :
    LagrangePreBasis(gv, (unsigned int)(k))
  {}

  //! Constructor for a given grid view object and run-time order
  LagrangePreBasis(const GridView& gv, unsigned int order) :
    Mixin(gv, dofLayout(order), Twist{gv.indexSet(), (unsigned int)(dofLayout(order)(Dune::GeometryTypes::line,dim))}),
    order_(order)
  {
    if (!useDynamicOrder && order!=(unsigned int)(k))
      DUNE_THROW(RangeError, "Template argument k has to be -1 when supplying a run-time order!");
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    return Node{order_};
  }

  template<typename It>
  It indices(const Node& node, It it) const
  {
    if (dim == 3 && order() > 3)
      DUNE_THROW(Dune::NotImplemented, "LagrangeBasis for 3D grids is only implemented if k<=3");

    return Mixin::indices(node, it);
  }

  //! Polynomial order used in the local Lagrange finite-elements
  unsigned int order() const
  {
    return order_;
  }

protected:
  // Run-time order, only valid if k<0
  unsigned int order_;
};



template<typename GV, int k, typename R>
class LagrangeNode :
  public LeafBasisNode
{
  // Stores LocalFiniteElement implementations with run-time order as a function of GeometryType
  template<typename Domain, typename Range, int dim>
  class LagrangeRunTimeLFECache
  {
  public:
    using FiniteElementType = LagrangeLocalFiniteElement<EquidistantPointSet,dim,Domain,Range>;

    const FiniteElementType& get(GeometryType type)
    {
      auto i = data_.find(type);
      if (i==data_.end())
        i = data_.emplace(type,FiniteElementType(type,order_)).first;
      return (*i).second;
    }

    std::map<GeometryType, FiniteElementType> data_;
    unsigned int order_;
  };

  static constexpr int dim = GV::dimension;
  static constexpr bool useDynamicOrder = (k<0);

  using FiniteElementCache = std::conditional_t<(useDynamicOrder),
                                                       LagrangeRunTimeLFECache<typename GV::ctype, R, dim>,
                                                       LagrangeLocalFiniteElementCache<typename GV::ctype, R, dim, std::max(k,0)>
                                                      >;

public:

  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename FiniteElementCache::FiniteElementType;

  //! Constructor without order (uses the compile-time value)
  LagrangeNode() :
    finiteElement_(nullptr),
    element_(nullptr)
  {}

  //! Constructor with a run-time order
  LagrangeNode(unsigned int order) :
    order_(order),
    finiteElement_(nullptr),
    element_(nullptr)
  {
    // Only the cache for the run-time-order case (i.e., k<0), has the 'order_' member
    if constexpr (useDynamicOrder)
      cache_.order_ = order;
  }

  //! Return current element, throw if unbound
  const Element& element() const
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement() const
  {
    return *finiteElement_;
  }

  //! Bind to element.
  void bind(const Element& e)
  {
    element_ = &e;
    finiteElement_ = &(cache_.get(element_->type()));
    this->setSize(finiteElement_->size());
  }

protected:

  unsigned int order() const
  {
    return (useDynamicOrder) ? order_ : k;
  }

  // Run-time order, only valid if k<0
  unsigned int order_;

  FiniteElementCache cache_;
  const FiniteElement* finiteElement_;
  const Element* element_;
};



namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a  Lagrange pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam k   The polynomial order of the ansatz functions; -1 means 'order determined at run-time'
 * \tparam R   The range type of the local basis
 */
template<std::size_t k, typename R=double>
auto lagrange()
{
  return [](const auto& gridView) {
    return LagrangePreBasis<std::decay_t<decltype(gridView)>, k, R>(gridView);
  };
}

/**
 * \brief Create a pre-basis factory that can create a  Lagrange pre-basis with a run-time order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R   The range type of the local basis
 */
template<typename R=double>
auto lagrange(int order)
{
  return [=](const auto& gridView) {
    return LagrangePreBasis<std::decay_t<decltype(gridView)>, -1, R>(gridView, order);
  };
}

} // end namespace BasisFactory



/** \brief Nodal basis of a scalar k-th-order Lagrangean finite element space
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \note This only works for certain grids.  The following restrictions hold
 * - If k is no larger than 2, then the grids can have any dimension
 * - If k is larger than 3 then the grid must be two-dimensional
 * - If k is 3, then the grid can be 3d *if* it is a simplex grid
 *
 * All arguments passed to the constructor will be forwarded to the constructor
 * of LagrangePreBasis.
 *
 * \warning The implementation of the basis with run-time order order uses the
 *   LagrangeFiniteElement implementation of dune-localfunctions, which is known
 *   to violate strict-aliasing rules
 *   (see https://gitlab.dune-project.org/core/dune-localfunctions/issues/14)
 *   Keep this in mind if ever you experience difficult-to-explain crashes
 *   or wrong results.
 *
 * \warning For pyramid elements in 3d, the shape functions are different for
 *    run-time and compile-time order. While the former are defined using the
 *    Duffy-transformation the latter are continuous and piecewise polynomial
 *    with discontinuous gradients.
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam k The order of the basis; -1 means 'order determined at run-time'
 * \tparam R The range type of the local basis
 */
template<typename GV, int k=-1, typename R=double>
using LagrangeBasis = DefaultGlobalBasis<LagrangePreBasis<GV, k, R> >;





} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LAGRANGEBASIS_HH
