// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LOBATTOBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LOBATTOBASIS_HH

#include <algorithm>
#include <map>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/localfunctions/lobatto.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>


namespace Dune {
namespace Functions {

// forward declaration
template<typename GV, typename R, typename Orders>
class LobattoNode;

/**
 * \brief A pre-basis for a Lobatto basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * A Lobatto basis is a set of continuous basis functions of local Lobatto shape functions
 * of arbitrary high order.
 *
 * The implementation is based on the local `Lobatto[Cube]LocalFiniteElement` implementing the
 * functions from
 *
 *   "Higher-Order Finite Element Methods", P. Soling, K, Segeth, I. Dolezel,
 *   2004, Chapman & Hall/CRC
 *
 * See \ref LobattoAdaptivePreBasis for an implementation with variable local polynomial degree.
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam MI  Type to be used for multi-indices
 * \tparam R   Range type used for shape function values
 */
template<typename GV, typename MI, typename R, typename Orders>
class LobattoPreBasis
{
  static const int dim = GV::dimension;

public:
  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Template mapping root tree path to type of created tree node
  using Node = LobattoNode<GV, R, Orders>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = MI;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<size_type, 1>;

  //! Constructor for a given grid view object and run-time order
  explicit LobattoPreBasis (const GridView& gv, unsigned int p = 1)
    : gridView_(gv)
  {
    for (auto type : gv.indexSet().types(0))
      orders_[type] = Orders{type, p};
  }

  //! Constructor for a given grid view object and run-time order
  LobattoPreBasis (const GridView& gv, const Orders& orders)
    : gridView_(gv)
  {
    for (auto type : gv.indexSet().types(0))
      orders_[type] = Orders{type, orders};
  }

  //! Initialize the global indices
  void initializeIndices ()
  {
    // number of DOFs per entity
    std::map<GeometryType, size_type> entityDofs;

    auto const& indexSet = gridView_.indexSet();
    for (auto type : indexSet.types(0)) {
      auto refElem = referenceElement<double,dim>(type);
      for (int c = 0; c <= dim; ++c) {
        for (int i = 0; i < refElem.size(c); ++i)
          entityDofs[refElem.type(i,c)] = orders_[type].size(i,c);
      }
    }

    size_type sum = 0;
    for (auto& ts : entityDofs) {
      offsets_[ts.first] = sum;
      sum += ts.second * (size_type)gridView_.size(ts.first);
    }

    size_ = sum;
    maxNodeSize_ = 0;
    for (auto t : indexSet.types(0))
      maxNodeSize_ = std::max(maxNodeSize_, size_type(orders_[t].size()));

    ready_ = true;
  }

  void debug () const
  {
    std::cout << "orders = " << orders_ << std::endl;

    std::cout << "offsets = {" << std::endl;
    for (auto&& ts : offsets_)
      std::cout << "  " << ts.first << " => " << ts.second << std::endl;
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
    return Node{[&](auto const& e) { return orders_.at(e.type()); }, gridView_.indexSet()};
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
    auto refElem = referenceElement(node.element());
    auto const& orders = orders_.at(refElem.type());

    for (size_type i = 0, end = localFE.size() ; i < end ; ++it, ++i)
    {
      Dune::LocalKey localKey = localCoefficients.localKey(i);
      size_type idx = gridIndexSet.subIndex(node.element(),localKey.subEntity(),localKey.codim());
      size_type entityDofs = orders.size(localKey.subEntity(), localKey.codim());
      GeometryType t = refElem.type(localKey.subEntity(), localKey.codim());

      *it = {{ offsets_.at(t) + entityDofs * idx + localKey.index() }};
    }
    return it;
  }

protected:
  GridView gridView_;
  std::map<GeometryType, Orders> orders_;
  std::map<GeometryType, size_type> offsets_;
  size_type size_;
  size_type maxNodeSize_;

  bool ready_ = false;
};



template<typename GV, typename R, typename Orders>
class LobattoNode
  : public LeafBasisNode
{
  static constexpr int dim = GV::dimension;
  using IndexSet = typename GV::IndexSet;

public:
  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;

  // TODO: use a GeometryType factory
  using FiniteElement = LobattoCubeLocalFiniteElement<typename GV::ctype, R, dim, Orders>;

  //! Default construction
  LobattoNode () = default;

  //! Constructor gets the vector or order information that is accessed in the bind method to
  //! build the corresponding local finite-element
  LobattoNode (std::function<Orders(const Element&)> orders, const IndexSet& indexSet)
    : orders_(std::move(orders))
    , indexSet_(&indexSet)
  {}

  //! Return current element, throw if unbound
  const Element& element () const
  {
    return *element_;
  }

  //! Return the LocalFiniteElement for the element we are bound to
  /**
   * The LocalFiniteElement implements the corresponding interfaces of dune-localfunctions
   */
  const FiniteElement& finiteElement () const
  {
    return *finiteElement_;
  }

  //! Bind to element.
  void bind (const Element& e)
  {
    element_ = &e;
    assert(e.type().isCube()); // currently only implemented for cubes.
    finiteElement_.emplace(orders_(e), Orientation<dim>{*element_, *indexSet_});
    this->setSize(finiteElement_->size());
  }

protected:
  std::function<Orders(const Element&)> orders_;
  const IndexSet* indexSet_ = nullptr;

  const Element* element_ = nullptr;
  std::optional<FiniteElement> finiteElement_ = std::nullopt;
};



namespace BasisFactory {
namespace Imp {

template <typename R, typename Orders>
class LobattoPreBasisFactory
{
public:
  static const std::size_t requiredMultiIndexSize = 1;

  LobattoPreBasisFactory (const Orders& orders)
    : orders_(orders)
  {}

  template <typename MultiIndex, typename GridView>
  auto makePreBasis (const GridView& gridView) const
  {
    if constexpr(!std::is_integral_v<Orders>)
      return LobattoPreBasis<GridView, MultiIndex, R, Orders>(gridView, orders_);
    else {
      using O = LobattoHomogeneousOrders<GridView::dimension>;
      return LobattoPreBasis<GridView, MultiIndex, R, O>(gridView, orders_);
    }
  }

  Orders orders_;
};

} // end namespace BasisFactory::Imp


/**
 * \brief Create a pre-basis factory that can create a  Lobatto pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R   The range type of the local basis
 * \param  p   Polynomial order on the elements
 */
template<typename R=double>
auto lobatto (unsigned int p = 1)
{
  return Imp::LobattoPreBasisFactory<R,unsigned int>{p};
}

/**
 * \brief Create a pre-basis factory that can create a  Lobatto pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R        The range type of the local basis
 * \param  orders   Polynomial orders on the elements
 */
template<typename R=double, int dim>
auto lobatto (const LobattoOrders<dim>& orders)
{
  return Imp::LobattoPreBasisFactory<R,LobattoOrders<dim>>{orders};
}

template<typename R=double, int dim>
auto lobatto (const LobattoHomogeneousOrders<dim>& orders)
{
  return Imp::LobattoPreBasisFactory<R,LobattoHomogeneousOrders<dim>>{orders};
}

} // end namespace BasisFactory


/** \brief Nodal basis of a scalar Lobatto finite element space
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * All arguments passed to the constructor will be forwarded to the constructor
 * of LobattoPreBasis.
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam R The range type of the local basis
 * \tparam Orders Type encoding the polynomial orders of the local shape functions
 */
template<typename GV, typename R=double, typename Orders=LobattoOrders<GV::dimension>>
using LobattoBasis = DefaultGlobalBasis<LobattoPreBasis<GV, FlatMultiIndex<std::size_t>, R, Orders> >;

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LOBATTOBASIS_HH