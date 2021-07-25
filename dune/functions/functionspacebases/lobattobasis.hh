// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LOBATTOBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LOBATTOBASIS_HH

#include <algorithm>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/localfunctions/lobatto.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>


namespace Dune {
namespace Functions {

// *****************************************************************************
// This is the reusable part of the LobattoBasis. It contains
//
//   LobattoPreBasis
//   LobattoNode
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

template<typename GV, typename R, typename Orders>
class LobattoNode;

template<typename GV, typename MI, typename R, typename Orders>
class LobattoPreBasis;


/**
 * \brief A pre-basis for a Lobatto basis
 *
 * \ingroup FunctionSpaceBasesImplementations
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
    : LobattoPreBasis{gv, Orders(std::uint8_t(p))}
  {}

  LobattoPreBasis (const GridView& gv, const Orders& orders)
    : gridView_(gv)
    , orders_(gv.size(0), orders)
  {}

  //! Initialize the global indices
  void initializeIndices ()
  {
    // number of DOFs per entity
    std::array<std::vector<size_type>, dim> entityDofs;

    // resize the vectors to number of entities and set number of DOFs to zero
    for (int c = 0; c < dim; ++c)
      entityDofs[c].resize(gridView_.size(c), 0);

    // traverse all elements and extract the local number of DOFs
    maxNodeSize_ = 0;
    auto const& indexSet = gridView_.indexSet();
    for (auto const& e : elements(gridView_)) {
      auto const& orders = orders_[indexSet.index(e)];
      maxNodeSize_ = std::max(maxNodeSize_, size_type(orders.size(e.type())));

      auto refElem = referenceElement(e);
      for (int c = 0; c < dim; ++c) {
        for (int i = 0; i < refElem.size(c); ++i) {
          entityDofs[c][indexSet.subIndex(e,i,c)] = orders.size(i,c);
        }
      }
    }

    // create partial sums of the DOF counts to obtain the offsets
    size_ = gridView_.size(dim);
    for (int d = 1; d <= dim; ++d) {
      int c = dim-d;
      offsets_[c].resize(gridView_.size(c));

      std::exclusive_scan(entityDofs[c].begin(), entityDofs[c].end(), offsets_[c].begin(), size_);
      size_ = offsets_[c].back() + entityDofs[c].back();
    }

    ready_ = true;
  }

  void debug() const
  {
    std::cout << "orders = {" << std::endl;
    for (std::size_t i = 0; i < orders_.size(); ++i)
      std::cout << "  " << i << ": " << orders_[i] << std::endl;
    std::cout << "}" << std::endl;

    std::cout << "offsets = {" << std::endl;
    for (int c = 0; c < dim; ++c) {
      std::cout << "  codim " << c << ":" << std::endl;
      for (std::size_t i = 0; i < offsets_[c].size(); ++i)
        std::cout << "    " << i << ": " << offsets_[c][i] << std::endl;
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

    // TODO: find a better way update the orders vector
    if (orders_.size() != std::size_t(gridView_.size(0))) {
      orders_.resize(gridView_.size(0));
      orders_.assign(orders_.size(), LobattoOrders<dim>{1});
    }
  }

  //! Create tree node
  Node makeNode () const
  {
    assert(ready_);
    return Node{orders_, gridView_.indexSet()};
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

    for (size_type i = 0, end = localFE.size() ; i < end ; ++it, ++i)
    {
      Dune::LocalKey localKey = localCoefficients.localKey(i);
      size_type idx = gridIndexSet.subIndex(node.element(),localKey.subEntity(),localKey.codim());

      /// TODO: for higher-orders we might need to transform the indices or basis functions
      /// depending on the orientation of the subEntity.
      *it = {{ offset(localKey.codim(), idx) + localKey.index() }};
    }
    return it;
  }

  /// Set new set of polynomial orders for the given entity.
  /// NOTE: Must call initializeIndices() afterwards.
  template<typename Entity>
  void setOrders (const Entity& entity, const Orders& orders)
  {
    orders_[gridView().indexSet().index(entity)] = orders;
    ready_ = false;
  }

private:
  size_type offset (unsigned int codim, size_type idx) const
  {
    return codim < dim ? offsets_[codim][idx] : idx;
  }

protected:
  GridView gridView_;
  std::vector<Orders> orders_;
  std::array<std::vector<size_type>, dim> offsets_;
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
  LobattoNode (const std::vector<Orders>& orders, const IndexSet& indexSet)
    : orders_(&orders)
    , indexSet_(&indexSet)
    , finiteElement_(std::nullopt)
    , element_(nullptr)
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
    finiteElement_.emplace((*orders_)[indexSet_->index(e)]);
    this->setSize(finiteElement_->size());
  }

protected:
  const std::vector<Orders>* orders_ = nullptr;
  const IndexSet* indexSet_ = nullptr;

  std::optional<FiniteElement> finiteElement_ = std::nullopt;
  const Element* element_ = nullptr;
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
template<typename R=double, unsigned int dim>
auto lobatto (const LobattoOrders<dim>& orders)
{
  return Imp::LobattoPreBasisFactory<R,LobattoOrders<dim>>{orders};
}

template<typename R=double, unsigned int dim>
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
