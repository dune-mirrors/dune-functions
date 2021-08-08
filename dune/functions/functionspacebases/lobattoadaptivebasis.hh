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

/**
 * \brief A pre-basis for an adaptive Lobatto basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam MI  Type to be used for multi-indices
 * \tparam R   Range type used for shape function values
 */
template<typename GV, typename MI, typename R, typename Orders>
class LobattoAdaptivePreBasis
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
  explicit LobattoAdaptivePreBasis (const GridView& gv, unsigned int p = 1)
    : gridView_(gv)
    , orders_(gv.size(0))
  {
    auto const& indexSet = gv.indexSet();
    for (auto const& e : elements(gv))
      orders_[indexSet.index(e)] = Orders{e.type(), p};
  }

  LobattoAdaptivePreBasis (const GridView& gv, const Orders& orders)
    : gridView_(gv)
    , orders_(gv.size(0))
  {
    auto const& indexSet = gv.indexSet();
    for (auto const& e : elements(gv))
      orders_[indexSet.index(e)] = Orders{e.type(), orders};
  }

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
      maxNodeSize_ = std::max(maxNodeSize_, size_type(orders.size()));

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

      auto it = offsets_[c].begin();
      for (auto const& d : entityDofs[c]) {
        *it++ = size_;
        size_ += d;
      }
    }

    ready_ = true;
  }

  void debug () const
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
      auto const& indexSet = gv.indexSet();
      for (auto const& e : elements(gv))
        orders_[indexSet.index(e)] = Orders{e.type(), 1};
    }
  }

  //! Create tree node
  Node makeNode () const
  {
    assert(ready_);
    return Node{[&](auto const& e) { return orders_[gridView_.indexSet().index(e)]; },
                gridView_.indexSet()};
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

      *it = {{ offset(localKey.codim(), idx) + localKey.index() }};
    }
    return it;
  }

  //! Set new set of polynomial orders for the given entity.
  // NOTE: Must call initializeIndices() afterwards.
  template<typename Entity>
  void setOrders (const Entity& entity, const Orders& orders)
  {
    const auto& indexSet = gridView().indexSet();
    orders_[indexSet.index(entity)] = orders;
    ready_ = false;
  }

  //! Enforce that the polynomial degree on entities shared by two elements is the minimum of the
  //! polynomial degrees of the element interiors.
  void enforceMinimumRule ()
  {
    const int dim = GridView::dimension;
    if constexpr(dim > 1) {
      const auto& indexSet = gridView().indexSet();
      for (auto const& e : elements(gridView())) {
        auto& inside = orders_[indexSet.index(e)];
        for (auto const& is : intersections(gridView(), e))
        {
          if (is.neighbor()) {
            auto& outside = orders_[indexSet.index(is.outside())];
            std::uint8_t min_p = std::min(inside(0,0,0), outside(0,0,0));

            for (int k = 0; k < dim-1; ++k) {
              outside.set(is.indexInOutside(),1,k,
                std::min(min_p, outside(is.indexInOutside(),1,k)));
              inside.set(is.indexInInside(),1,k,
                std::min(min_p, inside(is.indexInInside(),1,k)));
            }
            // TODO: set polynomial degree also for sub-entities of intersections
          }
        }
      }
      ready_ = false;
    }
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


namespace BasisFactory {
namespace Imp {

template <typename R, typename Orders>
class LobattoAdaptivePreBasisFactory
{
public:
  static const std::size_t requiredMultiIndexSize = 1;

  LobattoAdaptivePreBasisFactory (const Orders& orders)
    : orders_(orders)
  {}

  template <typename MultiIndex, typename GridView>
  auto makePreBasis (const GridView& gridView) const
  {
    if constexpr(!std::is_integral_v<Orders>)
      return LobattoAdaptivePreBasis<GridView, MultiIndex, R, Orders>(gridView, orders_);
    else {
      using O = LobattoHomogeneousOrders<GridView::dimension>;
      return LobattoAdaptivePreBasis<GridView, MultiIndex, R, O>(gridView, orders_);
    }
  }

  Orders orders_;
};

} // end namespace BasisFactory::Imp


/**
 * \brief Create a pre-basis factory that can create an adaptive Lobatto pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R   The range type of the local basis
 * \param  p   Polynomial order on the elements
 */
template<typename R=double>
auto lobattoAdaptive (unsigned int p = 1)
{
  return Imp::LobattoAdaptivePreBasisFactory<R,unsigned int>{p};
}

/**
 * \brief Create a pre-basis factory that can create an adaptive Lobatto pre-basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam R        The range type of the local basis
 * \param  orders   Polynomial orders on the elements
 */
template<typename R=double, int dim>
auto lobattoAdaptive (const LobattoOrders<dim>& orders)
{
  return Imp::LobattoAdaptivePreBasisFactory<R,LobattoOrders<dim>>{orders};
}

template<typename R=double, int dim>
auto lobattoAdaptive (const LobattoHomogeneousOrders<dim>& orders)
{
  return Imp::LobattoAdaptivePreBasisFactory<R,LobattoHomogeneousOrders<dim>>{orders};
}

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
 * \tparam Orders Type encoding the polynomial orders of the local shape functions
 */
template<typename GV, typename R=double, typename Orders=LobattoOrders<GV::dimension>>
using LobattoAdaptiveBasis = DefaultGlobalBasis<LobattoAdaptivePreBasis<GV, FlatMultiIndex<std::size_t>, R, Orders> >;

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LOBATTOADAPTIVEBASIS_HH
