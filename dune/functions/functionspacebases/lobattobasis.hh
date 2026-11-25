// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LOBATTOBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LOBATTOBASIS_HH

#include <algorithm>
#include <array>
#include <concepts>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>
#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/localfunctions/common/faceorientations.hh>
#include <dune/localfunctions/common/localfiniteelementvariant.hh>
#include <dune/localfunctions/lobatto.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/leafprebasismixin.hh>
#include <dune/functions/functionspacebases/nodes.hh>

namespace Dune::Functions {

namespace Impl {

  // Maximal polynomial order supported by the Lobatto local finite elements
  inline constexpr unsigned int lobattoMaxOrder = 19;

  // Number of Lobatto basis functions associated to an entity of the given type
  // with (isotropic) polynomial order p.
  inline std::size_t lobattoEntitySize (GeometryType type, unsigned int p)
  {
    const unsigned int m = type.dim();
    if (type.isVertex())
      return 1;
    else if (type.isCube())
      return Dune::power(std::size_t(p-1), m);
    else if (type.isSimplex())
      return (p > m) ? Dune::binomial(std::size_t(p-1), std::size_t(m)) : 0;
    else
      DUNE_THROW(Dune::NotImplemented, "Lobatto basis functions on " << type);
  }

  // Check that all elements of the grid view are supported by the Lobatto basis
  template<class GridView>
  void checkLobattoElementTypes (const GridView& gridView)
  {
    for (auto type : gridView.indexSet().types(0))
      if (!type.isCube() && !type.isSimplex())
        DUNE_THROW(Dune::NotImplemented, "The Lobatto basis is implemented for cube and simplex elements only, not for " << type);
  }

  // Global orientation of the edges and faces of an element induced by the global vertex ids
  template<class IdSet, class Element>
  auto lobattoFaceOrientations (const IdSet& idSet, const Element& element)
  {
    constexpr int dim = Element::mydimension;
    std::array<typename IdSet::IdType, (1 << dim)> ids;
    for (unsigned int k = 0; k < element.subEntities(dim); ++k)
      ids[k] = idSet.subId(element, k, dim);
    return Dune::Experimental::FaceOrientations<dim>(element.type(), ids, Dune::Indices::_1);
  }

} // end namespace Impl



// *****************************************************************************
// LobattoNode
// *****************************************************************************

/**
 * \brief Leaf basis node of the Lobatto bases
 *
 * The node creates the local finite element for each element from the polynomial
 * orders provided by the pre-basis and the global orientation of the edges and faces
 * induced by the global vertex ids. The orientation-independent part of the local finite
 * elements is cached for each set of orders.
 *
 * \tparam PB  The pre-basis, either \ref LobattoPreBasis or \ref LobattoAdaptivePreBasis
 */
template<class PB>
class LobattoNode
  : public LeafBasisNode
{
  using GridView = typename PB::GridView;
  static constexpr int dim = GridView::dimension;
  using Orders = typename PB::Orders;
  using R = typename PB::RangeFieldType;

  using CubeFE = LobattoCubeLocalFiniteElement<typename GridView::ctype, R, dim>;
  using SimplexFE = LobattoSimplexLocalFiniteElement<typename GridView::ctype, R, dim>;

  using SGT = Capabilities::hasSingleGeometryType<typename GridView::Grid>;
  static constexpr bool hasSingleCubeType = SGT::v && GeometryType(SGT::topologyId, dim).isCube();
  static constexpr bool hasSingleSimplexType = SGT::v && GeometryType(SGT::topologyId, dim).isSimplex();

public:
  using size_type = std::size_t;
  using Element = typename GridView::template Codim<0>::Entity;

  //! The local finite element type
  using FiniteElement = std::conditional_t<hasSingleCubeType, CubeFE,
                        std::conditional_t<hasSingleSimplexType, SimplexFE,
                          LocalFiniteElementVariant<CubeFE, SimplexFE>>>;

  explicit LobattoNode (const PB& preBasis)
    : preBasis_(&preBasis)
  {}

  //! Return current element, throw if unbound
  const Element& element () const
  {
    return *element_;
  }

  //! Return the local finite element for the element we are bound to
  const FiniteElement& finiteElement () const
  {
    return *finiteElement_;
  }

  //! Bind to element
  void bind (const Element& e)
  {
    element_ = &e;
    const auto orientations = Impl::lobattoFaceOrientations(preBasis_->gridView().grid().globalIdSet(), e);
    const Orders orders = preBasis_->localOrders(e);
    if (e.type().isCube())
      setFiniteElement(cubeCache_, orders, orientations);
    else if (e.type().isSimplex())
      setFiniteElement(simplexCache_, orders, orientations);
    else
      DUNE_THROW(Dune::NotImplemented, "The Lobatto basis is implemented for cube and simplex elements only, not for " << e.type());
    this->setSize(finiteElement_->size());
  }

private:
  template<class FE, class Orientations>
  void setFiniteElement (std::map<Orders,FE>& cache, const Orders& orders, const Orientations& orientations)
  {
    auto it = cache.find(orders);
    if (it == cache.end())
      it = cache.emplace(orders, FE(orders)).first;

    if constexpr(std::is_same_v<FE, FiniteElement>)
      finiteElement_.emplace(it->second, orientations);
    else if constexpr(std::is_constructible_v<FiniteElement, FE>)
      finiteElement_.emplace(FE(it->second, orientations));
    else
      DUNE_THROW(Dune::Exception, "Element type does not match the grid capabilities");
  }

  const PB* preBasis_;
  const Element* element_ = nullptr;
  std::optional<FiniteElement> finiteElement_;

  // Finite elements with default orientation for all orders used so far
  std::map<Orders, CubeFE> cubeCache_;
  std::map<Orders, SimplexFE> simplexCache_;
};



// *****************************************************************************
// LobattoPreBasis
// *****************************************************************************

/**
 * \brief A pre-basis for a Lobatto basis with uniform polynomial order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * A continuous basis of hierarchic Lobatto shape functions of order `p` on grids
 * with simplex and cube elements, following
 *
 *   P. Šolín, K. Segeth, I. Doležel: "Higher-Order Finite Element Methods",
 *   Chapman & Hall/CRC, 2004.
 *
 * In contrast to Lagrange bases, the basis functions are not nodal. Instead, there are
 * vertex functions (the linear Lagrange functions) and edge, face, and interior functions
 * of increasing polynomial degree. The orientation of edge and face functions is determined
 * by the global vertex ids of the grid.
 *
 * See \ref LobattoAdaptivePreBasis for a basis with variable polynomial order per element.
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam R   Range field type used for shape function values
 */
template<class GV, class R = double>
class LobattoPreBasis
  : public LeafPreBasisMapperMixin<GV>
{
  using Base = LeafPreBasisMapperMixin<GV>;
  static constexpr int dim = GV::dimension;

  static Dune::MCMGLayout dofLayout (unsigned int p)
  {
    return [p](Dune::GeometryType type, int) -> std::size_t {
      return (type.isCube() || type.isSimplex()) ? Impl::lobattoEntitySize(type, p) : 0;
    };
  }

public:
  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Range field type of the shape functions
  using RangeFieldType = R;

  //! Type representing the polynomial orders on an element
  using Orders = LobattoHomogeneousOrders<dim>;

  //! Template mapping root tree path to type of created tree node
  using Node = LobattoNode<LobattoPreBasis>;

  //! Constructor for a given grid view object and polynomial order `p`
  LobattoPreBasis (const GridView& gv, unsigned int p)
    : Base(gv, dofLayout(p))
    , order_(p)
  {
    if (p < 1 || p > Impl::lobattoMaxOrder)
      DUNE_THROW(Dune::RangeError, "Lobatto order " << p << " not in the supported range 1,...," << Impl::lobattoMaxOrder);
    Impl::checkLobattoElementTypes(gv);
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update (const GridView& gv)
  {
    Impl::checkLobattoElementTypes(gv);
    Base::update(gv);
  }

  //! Create tree node
  Node makeNode () const
  {
    return Node{*this};
  }

  //! Polynomial order on all elements
  unsigned int order () const
  {
    return order_;
  }

  //! Polynomial orders of the local finite element on the given element
  template<class Element>
  Orders localOrders (const Element&) const
  {
    return Orders(order_);
  }

private:
  unsigned int order_;
};



// *****************************************************************************
// LobattoAdaptivePreBasis
// *****************************************************************************

/**
 * \brief A pre-basis for a Lobatto basis with variable polynomial order per element
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * The polynomial order of each element is given by a user-provided function
 * `elementOrder(element)` that returns either a single order or an array of orders
 * for the coordinate directions (anisotropic orders are supported on cube elements
 * only). This function is evaluated whenever the indices are initialized, i.e.,
 * on construction of the global basis and on `update()`.
 *
 * The orders of the edges and faces are determined by the minimum rule
 * (Šolín et al., Section 3.5.5): the order of a sub-entity is the minimum of the
 * orders of all adjacent elements. For cube elements with anisotropic orders, only
 * the orders in the directions tangential to the sub-entity are considered.
 * Optionally, a function `subEntityOrder(element, i, c)` can be given that
 * overrides the order of the `i`th sub-entity of codimension `c` (with
 * `0 < c < dim`) as seen from the element. The minimum over all adjacent elements
 * of these values is used.
 *
 * \note The minimum rule is evaluated on the local process only. Thus, the basis
 *       is not yet suitable for p-adaptivity on distributed grids.
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam R   Range field type used for shape function values
 */
template<class GV, class R = double>
class LobattoAdaptivePreBasis
  : public LeafPreBasisMixin<LobattoAdaptivePreBasis<GV,R>>
{
  static constexpr int dim = GV::dimension;

public:
  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Range field type of the shape functions
  using RangeFieldType = R;

  //! Type representing the polynomial orders on an element
  using Orders = LobattoOrders<dim>;

  //! Template mapping root tree path to type of created tree node
  using Node = LobattoNode<LobattoAdaptivePreBasis>;

  //! The grid element type
  using Element = typename GV::template Codim<0>::Entity;

  //! Function returning the (directional) orders of an element
  using ElementOrderFunction = std::function<std::array<unsigned int,dim>(const Element&)>;

  //! Function returning the order of the `i`th sub-entity of codim `c` of an element
  using SubEntityOrderFunction = std::function<unsigned int(const Element&, int, int)>;

  /**
   * \brief Constructor for a given grid view and element orders
   *
   * \param gv              The grid view
   * \param elementOrder    Function `element -> p` or `element -> std::array<unsigned int,dim>`
   * \param subEntityOrder  Optional function `(element, i, c) -> p` overriding the minimum rule
   **/
  template<class EO, class SO = SubEntityOrderFunction>
  LobattoAdaptivePreBasis (const GridView& gv, EO&& elementOrder, SO&& subEntityOrder = {})
    : gridView_(gv)
    , mapper_(gv, [](Dune::GeometryType, int) -> std::size_t { return 1; })
    , elementOrder_(makeElementOrderFunction(std::forward<EO>(elementOrder)))
    , subEntityOrder_(std::forward<SO>(subEntityOrder))
  {
    Impl::checkLobattoElementTypes(gv);
  }

  //! Evaluate the element orders, apply the minimum rule, and compute the global indices
  void initializeIndices ()
  {
    constexpr std::uint8_t unset = std::numeric_limits<std::uint8_t>::max();
    orders_.assign(mapper_.size(), unset);
    elementOrders_.assign(mapper_.size(), {});

    // element orders
    for (const auto& e : elements(gridView_)) {
      const auto p = elementOrder_(e);
      for (int k = 0; k < dim; ++k) {
        if (p[k] < 1 || p[k] > Impl::lobattoMaxOrder)
          DUNE_THROW(Dune::RangeError, "Lobatto order " << p[k] << " not in the supported range 1,...," << Impl::lobattoMaxOrder);
        if (e.type().isSimplex() && p[k] != p[0])
          DUNE_THROW(Dune::NotImplemented, "Anisotropic orders on simplex elements");
      }
      auto& pe = elementOrders_[mapper_.index(e)];
      for (int k = 0; k < dim; ++k)
        pe[k] = std::uint8_t(p[k]);
    }

    // minimum rule for the sub-entities
    for (const auto& e : elements(gridView_)) {
      const auto& re = referenceElement<double,dim>(e.type());
      const auto& pe = elementOrders_[mapper_.index(e)];
      for (int c = 1; c < dim; ++c) {
        for (int i = 0; i < re.size(c); ++i) {
          unsigned int q = subEntityOrder_ ? subEntityOrder_(e, i, c) : tangentialOrder(e.type(), pe, i, c);
          if (q < 1 || q > Impl::lobattoMaxOrder)
            DUNE_THROW(Dune::RangeError, "Lobatto order " << q << " not in the supported range 1,...," << Impl::lobattoMaxOrder);
          auto& ps = orders_[mapper_.subIndex(e, i, c)];
          ps = std::min(ps, std::uint8_t(q));
        }
      }
    }

    // number of basis functions per entity
    std::vector<size_type> sizes(mapper_.size(), 0);
    maxNodeSize_ = 0;
    for (const auto& e : elements(gridView_)) {
      const auto& re = referenceElement<double,dim>(e.type());
      size_type nodeSize = 0;
      for (int c = 0; c <= dim; ++c) {
        for (int i = 0; i < re.size(c); ++i) {
          const auto idx = mapper_.subIndex(e, i, c);
          sizes[idx] = entitySize(e.type(), re.type(i,c), idx, c);
          nodeSize += sizes[idx];
        }
      }
      maxNodeSize_ = std::max(maxNodeSize_, nodeSize);
    }

    // offsets as exclusive prefix sum of the sizes
    offsets_.resize(mapper_.size() + 1);
    offsets_[0] = 0;
    for (size_type idx = 0; idx < mapper_.size(); ++idx)
      offsets_[idx+1] = offsets_[idx] + sizes[idx];
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView () const
  {
    return gridView_;
  }

  //! Update the stored grid view, to be called if the grid or the orders have changed
  void update (const GridView& gv)
  {
    Impl::checkLobattoElementTypes(gv);
    gridView_ = gv;
    mapper_.update(gridView_);
  }

  //! Create tree node
  Node makeNode () const
  {
    return Node{*this};
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension () const
  {
    return offsets_.empty() ? 0 : offsets_.back();
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize () const
  {
    return maxNodeSize_;
  }

  //! Polynomial orders of the local finite element on the given element
  Orders localOrders (const Element& e) const
  {
    const auto& re = referenceElement<double,dim>(e.type());
    Orders orders;
    const auto& pe = elementOrders_[mapper_.index(e)];
    for (int k = 0; k < dim; ++k)
      orders.setInterior(k, pe[k]);
    for (int c = 1; c < dim; ++c)
      for (int i = 0; i < re.size(c); ++i)
        orders.set(i, c, orders_[mapper_.subIndex(e, i, c)]);
    return orders;
  }

  //! Polynomial order of the `i`th sub-entity of codimension `0 < c < dim` of the element (after the minimum rule)
  unsigned int subEntityOrder (const Element& e, int i, int c) const
  {
    return orders_[mapper_.subIndex(e, i, c)];
  }

  //! Maximal polynomial order over all elements
  unsigned int maxOrder () const
  {
    unsigned int p = 1;
    for (const auto& pe : elementOrders_)
      for (auto pk : pe)
        p = std::max(p, (unsigned int)(pk));
    return p;
  }

  //! Collect the global indices of the basis functions on a node
  template<class Node, class It>
  It indices (const Node& node, It it) const
  {
    const auto& localCoefficients = node.finiteElement().localCoefficients();
    for (size_type i = 0; i < localCoefficients.size(); ++i, ++it) {
      const auto& localKey = localCoefficients.localKey(i);
      *it = {{ offsets_[mapper_.subIndex(node.element(), localKey.subEntity(), localKey.codim())] + localKey.index() }};
    }
    return it;
  }

private:
  template<class EO>
  static ElementOrderFunction makeElementOrderFunction (EO&& elementOrder)
  {
    return [f=std::forward<EO>(elementOrder)](const Element& e) {
      if constexpr(std::is_convertible_v<std::invoke_result_t<const EO&, const Element&>, unsigned int>) {
        std::array<unsigned int,dim> p;
        p.fill(f(e));
        return p;
      }
      else
        return std::array<unsigned int,dim>(f(e));
    };
  }

  // Order of the `i`th sub-entity of codim `c` induced by the element orders `pe`:
  // the minimum of the orders in the tangential directions of the sub-entity.
  static unsigned int tangentialOrder (GeometryType type, const std::array<std::uint8_t,dim>& pe, int i, int c)
  {
    if (!type.isCube())
      return pe[0];
    const auto& re = referenceElement<double,dim>(type);
    const auto x0 = re.position(re.subEntity(i,c,0,dim), dim);
    unsigned int q = std::numeric_limits<unsigned int>::max();
    for (int k = 0; k < re.size(i,c,dim); ++k) {
      const auto xk = re.position(re.subEntity(i,c,k,dim), dim);
      for (int d = 0; d < dim; ++d)
        if (xk[d] != x0[d])
          q = std::min(q, (unsigned int)(pe[d]));
    }
    return q;
  }

  // Number of basis functions associated to the entity with mapper index `idx`
  size_type entitySize (GeometryType elementType, GeometryType type, size_type idx, int c) const
  {
    if (c == dim)
      return 1;
    else if (c == 0) {
      const auto& pe = elementOrders_[idx];
      if (elementType.isCube()) {
        size_type s = 1;
        for (int k = 0; k < dim; ++k)
          s *= pe[k] - 1;
        return s;
      }
      return Impl::lobattoEntitySize(type, pe[0]);
    }
    else
      return Impl::lobattoEntitySize(type, orders_[idx]);
  }

  GridView gridView_;
  Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper_;
  ElementOrderFunction elementOrder_;
  SubEntityOrderFunction subEntityOrder_;

  std::vector<std::array<std::uint8_t,dim>> elementOrders_;  // indexed by the mapper
  std::vector<std::uint8_t> orders_;                          // orders of sub-entities with 0 < codim < dim
  std::vector<size_type> offsets_;
  size_type maxNodeSize_ = 0;
};



namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a Lobatto pre-basis of uniform order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R  The range field type of the local basis
 * \param  p  The polynomial order
 */
template<class R = double>
auto lobatto (unsigned int p)
{
  return [p](const auto& gridView) {
    return LobattoPreBasis<std::decay_t<decltype(gridView)>, R>(gridView, p);
  };
}

/**
 * \brief Create a pre-basis factory that can create a Lobatto pre-basis with variable order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R               The range field type of the local basis
 * \param  elementOrder    Function `element -> p` or `element -> std::array<unsigned int,dim>`
 * \param  subEntityOrder  Optional function `(element, i, c) -> p` overriding the minimum rule
 *
 * See \ref LobattoAdaptivePreBasis for details.
 */
template<class R = double, class EO, class... SO>
  requires (!std::is_arithmetic_v<std::decay_t<EO>>) && (sizeof...(SO) <= 1)
auto lobatto (EO&& elementOrder, SO&&... subEntityOrder)
{
  return [=](const auto& gridView) {
    return LobattoAdaptivePreBasis<std::decay_t<decltype(gridView)>, R>(gridView, elementOrder, subEntityOrder...);
  };
}

} // end namespace BasisFactory



/**
 * \brief Basis of a scalar Lobatto finite element space with uniform order
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The GridView that the space is defined on
 * \tparam R   The range field type of the local basis
 */
template<class GV, class R = double>
using LobattoBasis = DefaultGlobalBasis<LobattoPreBasis<GV, R>>;

/**
 * \brief Basis of a scalar Lobatto finite element space with variable order per element
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The GridView that the space is defined on
 * \tparam R   The range field type of the local basis
 */
template<class GV, class R = double>
using LobattoAdaptiveBasis = DefaultGlobalBasis<LobattoAdaptivePreBasis<GV, R>>;

} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LOBATTOBASIS_HH
