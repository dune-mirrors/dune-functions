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
   * \param orders  For each codim a vector associating each entity a polynomial degree
   **/
  LobattoAdaptivePreBasis (const GridView& gv,
                           std::array<std::vector<std::uint8_t>, dim> const* orders)
    : gridView_(gv)
    , orders_(orders)
  {}

  //! Construct an order container on the element, by collecting polynomial orders on
  //! all its subentities from the global vectors
  template <class Entity>
  LobattoOrders<dim> getOrder (const Entity& e) const
  {
    LobattoOrders<dim> order{e.type(),1};

    auto const& indexSet = gridView_.indexSet();
    auto refElem = referenceElement(e);
    for (int c = 0; c < dim; ++c) {
      for (int i = 0; i < refElem.size(c); ++i) {
        unsigned int p = (*orders_)[c].empty() ? 0u : (*orders_)[c][indexSet.subIndex(e,i,c)];
        for (int k = 0; k < dim-c; ++k)
          order.set(i,c,k, p); // TODO: distinguish directions
      }
    }

    return order;
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
      LobattoOrders<dim> order = getOrder(e);
      maxNodeSize_ = std::max(maxNodeSize_, size_type(order.size()));

      auto refElem = referenceElement(e);
      for (int c = 0; c < dim; ++c)
        for (int i = 0; i < refElem.size(c); ++i)
          entityDofs[c][indexSet.subIndex(e,i,c)] = order.size(i,c);
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

    for (size_type i = 0, end = localFE.size() ; i < end ; ++it, ++i)
    {
      Dune::LocalKey localKey = localCoefficients.localKey(i);
      size_type idx = gridIndexSet.subIndex(node.element(),localKey.subEntity(),localKey.codim());

      *it = {{ offset(localKey.codim(), idx) + localKey.index() }};
    }
    return it;
  }


private:
  size_type offset (unsigned int codim, size_type idx) const
  {
    return codim < dim ? offsets_[codim][idx] : idx;
  }

protected:
  GridView gridView_;
  const std::array<std::vector<std::uint8_t>, dim>* orders_;
  std::array<std::vector<size_type>, dim> offsets_;
  size_type size_;
  size_type maxNodeSize_;

  bool ready_ = false;
};


namespace BasisFactory {
namespace Imp {

template <typename R, std::size_t dim>
class LobattoAdaptivePreBasisFactory
{
public:
  static const std::size_t requiredMultiIndexSize = 1;

  LobattoAdaptivePreBasisFactory (const std::array<std::vector<std::uint8_t>, dim>* orders)
    : orders_(orders)
  {}

  template <typename MultiIndex, typename GridView>
  auto makePreBasis (const GridView& gridView) const
  {
    return LobattoAdaptivePreBasis<GridView, MultiIndex, R>(gridView, orders_);
  }

  const std::array<std::vector<std::uint8_t>, dim>* orders_;
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
template<typename R=double, std::size_t dim>
auto lobatto (const std::array<std::vector<std::uint8_t>, dim>& orders)
{
  return Imp::LobattoAdaptivePreBasisFactory<R,dim>{&orders};
}

template<typename R=double, std::size_t dim>
void lobatto (const std::array<std::vector<std::uint8_t>, dim>&& orders) = delete;

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
using LobattoAdaptiveBasis = DefaultGlobalBasis<LobattoAdaptivePreBasis<GV, FlatMultiIndex<std::size_t>, R> >;

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LOBATTOADAPTIVEBASIS_HH
