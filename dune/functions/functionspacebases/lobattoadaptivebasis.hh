// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LOBATTOADAPTIVEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LOBATTOADAPTIVEBASIS_HH

#include <algorithm>
#include <map>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/localfunctions/lobatto.hh>
#include <dune/functions/functionspacebases/lobattobasis.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>


namespace Dune {
namespace Functions {

// For each GeometryType the association of polynomial degree to entity index
template<typename GridView>
class LobattoEntityOrders
{
  using IndexSet = typename GridView::IndexSet;
  using IndexType = typename IndexSet::IndexType;

public:
  /// Store a pointer to the indexSet of the given GridView
  LobattoEntityOrders (const GridView& gridView)
  {
    update(gridView);
  }

  /// On grid changes update the GridView and indexSet. Resets all stored polynomial degrees
  void update (const GridView& gridView)
  {
    indexSet_ = &gridView.indexSet();
    clear();
  }

  /// Resets all stored polynomial degrees
  void clear ()
  {
    orders_.clear();
    maxOrders_.clear();
    maxOrder_ = 1;
  }

  /// Set the polynomial degree of the given `entity` to `p`
  template<typename Entity>
  void set (const Entity& entity, std::uint8_t p)
  {
    set(entity.type(), indexSet_->index(entity), p);
  }

  /// Set the polynomial degree of the given entity with GeometryType `t` and index `idx` to `p`
  void set (GeometryType t, IndexType idx, std::uint8_t p)
  {
    if (!contains(t))
      orders_[t].resize(indexSet_->size(t), 1u);

    orders_[t][idx] = p;
    maxOrders_[t] = std::max(maxOrders_[t], p);
    maxOrder_ = std::max(maxOrder_, p);
  }

  /// Set the polynomial degree of all entities with GeometryType `t` to `p`
  void setAll (GeometryType t, std::uint8_t p)
  {
    orders_[t].clear();
    orders_[t].resize(indexSet_->size(t), p);
    maxOrders_[t] = std::max(maxOrders_[t], p);
    maxOrder_ = std::max(maxOrder_, p);
  }

  /// Return the polynomial degree for the given `entity`
  template<typename Entity>
  std::uint8_t get (const Entity& entity) const
  {
    return get(entity.type(), indexSet_->index(entity));
  }

  /// Return the polynomial degree for the given entity with GeometryType `t` and index `idx`
  /**
   * The polynomial degree is always >= 1. Thus, if no other polynomial degree is stored
   * for the given GeometryType and index, return the value `1`.
   **/
  std::uint8_t get (GeometryType t, IndexType idx) const
  {
    return contains(t) ? orders_.at(t).at(idx) : 1u;
  }

  /// Return the maximal polynomial order over all entities
  std::uint8_t max () const
  {
    return maxOrder_;
  }

  /// Return the maximal polynomial order over all entities of GeometryType `t`
  std::uint8_t max (GeometryType t) const
  {
    return maxOrders_.at(t);
  }

  /// Check whether a polynomial degree is stored for GeometryType `t`
  bool contains (GeometryType t) const
  {
    return orders_.count(t) > 0;
  }

  /// Begin-iterator over all orders
  auto begin () const { return orders_.begin(); }

  /// End-iterator over all orders
  auto end () const { return orders_.end(); }

private:
  const IndexSet* indexSet_ = nullptr;
  std::map<GeometryType, std::vector<std::uint8_t>> orders_;
  std::map<GeometryType, std::uint8_t> maxOrders_;
  std::uint8_t maxOrder_ = 1;
};


/**
 * \brief A pre-basis for an adaptive Lobatto basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * An adaptive Lobatto basis is a set of basis functions of local Lobatto shape functions
 * that allows to have different polynomial orders on each element. It even allows to have
 * different orders on sub-entities and thus opens the possibility for p-adaptivity.
 *
 * Note, for a stable discretization you need to fulfill a minimum-rule that essentially means that
 * the polynomial degree on sub-entities shared by some entities is a most the minimum of the
 * polynomial degrees of the entities.
 *
 * The implementation is based on the local `Lobatto[Cube]LocalFiniteElement` implementing the
 * functions from
 *
 *   "Higher-Order Finite Element Methods", P. Soling, K, Segeth, I. Dolezel,
 *   2004, Chapman & Hall/CRC
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam MI  Type to be used for multi-indices
 * \tparam R   Range type used for shape function values
 */
template<typename GV, typename MI, typename R>
class LobattoAdaptivePreBasis
{
  static const int dim = GV::dimension;

public:
  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  using Node = LobattoNode<GV, R, LobattoOrders<dim>>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  //! Constructor for a given grid view object and order vectors
  /**
   * \param gv      GridView to associate the basis to
   * \param orders  For each GeometryType a vector associating each entity a polynomial degree
   **/
  LobattoAdaptivePreBasis (const GridView& gv, const LobattoEntityOrders<GV>* orders)
    : gridView_(gv)
    , orders_(orders)
  {}

  // converting constructor
  template<typename MI_, typename R_>
  LobattoAdaptivePreBasis (const LobattoAdaptivePreBasis<GV,MI_,R_>& other)
    : gridView_(other.gridView_)
    , orders_(other.orders_)
  {}

  //! Construct an order container on the element, by collecting polynomial orders on
  //! all its subentities from the global vectors
  template<typename Entity>
  LobattoOrders<dim> getOrder (const Entity& e) const
  {
    assert(!!orders_);
    LobattoOrders<dim> order{e.type(),1};

    auto const& indexSet = gridView_.indexSet();
    auto refElem = referenceElement(e);
    for (int c = 0; c < dim; ++c) {
      for (int i = 0; i < refElem.size(c); ++i) {
        const auto type = refElem.type(i,c);
        if (orders_->contains(type)) {
          unsigned int p =  orders_->get(type, indexSet.subIndex(e,i,c));
          for (int k = 0; k < dim-c; ++k)
            order.set(i,c,k, p); // TODO: distinguish directions
        }
      }
    }

    return order;
  }

  //! Initialize the global indices
  void initializeIndices ()
  {
    std::map<GeometryType, std::map<std::uint8_t, size_type>> sizes; // [GeometryType][p]

    // count the total number of DOFs of associated to entity of codim c and polynomial degree p
    for (auto const& o : *orders_)
      for (std::uint8_t p : o.second)
      sizes[o.first][p] += LobattoGeometry::size(o.first,p,p,p);

    // transform sizes into entity offsets
    size_type offset = gridView_.size(dim); // number of vertices
    for (auto& st : sizes) {
      for (auto& stp : st.second) {
        auto size = stp.second;
        stp.second = offset;
        offset += size;
      }
    }

    // copmute global offsets
    for (auto const& o : *orders_) {
      GeometryType t = o.first;
      offsets_[t].resize(gridView_.size(t));
      for (std::size_t i = 0; i < std::size_t(gridView_.size(t)); ++i) {
        std::uint8_t p = orders_->get(t,i);
        offsets_[t][i] = sizes[t][p];
        sizes[t][p] += LobattoGeometry::size(t,p,p,p);
      }
    }

    // compute upper bound for node size
    maxNodeSize_ = 0;
    for (auto const& t : gridView_.indexSet().types(0)) {
      LobattoOrders<dim> maxOrders{t,1};
      auto refElem = referenceElement<double,dim>(t);
      for (int c = 0; c < dim; ++c) {
        for (int i = 0; i < refElem.size(c); ++i) {
          const auto type = refElem.type(i,c);
          if (orders_->contains(type)) {
            unsigned int p =  orders_->max(type);
            for (int k = 0; k < dim-c; ++k)
              maxOrders.set(i,c,k, p); // TODO: distinguish directions
          }
        }
      }
      maxNodeSize_ = std::max(maxNodeSize_, size_type(maxOrders.size()));
    }

    size_ = offset;
    ready_ = true;
  }

  void debug () const
  {
    std::cout << "offsets = {" << std::endl;
    for (auto const& st : offsets_) {
      std::cout << "  type " << st.first << ": " << std::endl;
      for (auto const& idx : st.second)
        std::cout << "    " << idx << std::endl;
    }
    std::cout << "}" << std::endl;

    std::cout << "size = " << size_ << std::endl;
    std::cout << "maxNodeSize = " << maxNodeSize_ << std::endl;
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView () const
  {
    return gridView_;
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update (const GridView& gv)
  {
    gridView_ = gv;
  }

  //! Create tree node
  Node makeNode () const
  {
    assert(ready_);
    return Node{[&](auto const& e) { return this->getOrder(e); }, gridView_.indexSet()};
  }

  //! Same as size(prefix) with empty prefix
  size_type size () const
  {
    assert(ready_);
    return size_;
  }

  //! Return number of possible values for next position in multi index
  size_type size (const SizePrefix prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ? size() : 0;
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension () const
  {
    return size();
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize () const
  {
    assert(ready_);
    return maxNodeSize_;
  }

  //! Collect the local indices on a leaf-basis node.
  template<typename It>
  It indices (const Node& node, It it) const
  {
    assert(ready_);
    const auto& gridIndexSet = gridView().indexSet();
    const auto& localFE = node.finiteElement();
    const auto& localCoefficients = localFE.localCoefficients();

    const auto refElem = referenceElement(node.element());

    for (size_type i = 0, end = localFE.size() ; i < end ; ++it, ++i)
    {
      Dune::LocalKey localKey = localCoefficients.localKey(i);
      size_type idx = gridIndexSet.subIndex(node.element(),localKey.subEntity(),localKey.codim());

      const auto type = refElem.type(localKey.subEntity(),localKey.codim());
      *it = {{ offset(type, idx) + localKey.index() }};
    }

    return it;
  }

private:
  size_type offset (GeometryType t, size_type idx) const
  {
    return t.isVertex() ? idx : offsets_.at(t)[idx];
  }

public:
  const LobattoEntityOrders<GV>* orders () const
  {
    return orders_;
  }

protected:
  GridView gridView_;
  const LobattoEntityOrders<GV>* orders_;

  std::map<GeometryType, std::vector<size_type>> offsets_;
  size_type size_ = 0;
  size_type maxNodeSize_ = 0;
  bool ready_ = false;
};


namespace BasisFactory {
namespace Imp {

template<typename R, typename GV>
class LobattoAdaptivePreBasisFactory
{
public:
  static const std::size_t requiredMultiIndexSize = 1;

  LobattoAdaptivePreBasisFactory (const LobattoEntityOrders<GV>* orders)
    : orders_(orders)
  {}

  template<typename MultiIndex>
  auto makePreBasis (const GV& gridView) const
  {
    return LobattoAdaptivePreBasis<GV, MultiIndex, R>(gridView, orders_);
  }

  const LobattoEntityOrders<GV>* orders_;
};

} // end namespace BasisFactory::Imp


/**
 * \brief Create a pre-basis factory that can create an adaptive Lobatto pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R       The range type of the local basis
 * \param  orders  Vector or polynomial orders for each codimension
 */
template<typename R=double, typename GV>
auto lobatto (const LobattoEntityOrders<GV>& orders)
{
  return Imp::LobattoAdaptivePreBasisFactory<R,GV>{&orders};
}

template<typename R=double, typename GV>
void lobatto (const LobattoEntityOrders<GV>&& orders) = delete;

} // end namespace BasisFactory


/** \brief Nodal basis of a scalar adaptive Lobatto finite element space
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * All arguments passed to the constructor will be forwarded to the constructor
 * of LobattoAdaptivePreBasis.
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam R The range type of the local basis
 */
template<typename GV, typename R=double>
using LobattoAdaptiveBasis
  = DefaultGlobalBasis<LobattoAdaptivePreBasis<GV, FlatMultiIndex<std::size_t>, R> >;

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LOBATTOADAPTIVEBASIS_HH
