// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASIXBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASIXBASIS_HH

#include <array>
#include <cassert>
#include <iterator>
#include <span>
#include <type_traits>

#include <dune/localfunctions/basix/globalbasix.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>

#include <dune/geometry/type.hh>

namespace Dune::Functions {

/**
  * \brief A pre-basis based on the basix framework
  *
  * \ingroup FunctionSpaceBasesImplementations
  *
  * \tparam GV        The grid view that the FE basis is defined on.
  * \tparam rangeClass  The class of the basis function range, i.e., RangeClass:scalar, vector or matrix.
  */
template<class GV, RangeClass rangeClass>
class BasixPreBasis
    : public LeafPreBasisMapperMixin<GV>
{
  using Base = LFEPreBasisMixin<GV>;
  using Geometry = typename GV::template Codim<0>::Entity::Geometry;
  using FiniteElement = BasixFiniteElement<Geometry,rangeClass,typename GV::ctype>;
  using Basix = typename FiniteElement::Basix;

  static auto makeLayout (const Basix& b)
  {
    auto& entity_dofs = b.entity_dofs();
    assert(b.cell_type() != basix::element::cell::type::prism);
    assert(b.cell_type() != basix::element::cell::type::pyramid);
    // TODO: need to count the DOFs per GeometryType and not per dimension
    std::array<std::size_t, GV::dimension+1> sizes{};
    for (std::size_t d = 0; d < entity_dofs.size(); ++d) {
      sizes[d] = entity_dofs[d][0].size();
    }
    return [sizes](GeometryType gt, int) -> std::size_t { return sizes[gt.dim()]; };
  }

  class Node;

public:
  using size_type = typename Base::size_type;

public:
  BasixPreBasis (const GV& gridView, const Basix& basix)
    : Base(gridView, makeLayout(basix))
    , fe_{basix}
  {
    // TODO: It probably only works for single GeometryType grids
  }

  //! Create tree node
  Node makeNode () const
  {
    return Node(fe_);
  }

  //! Fill cache with global indices of DOFs associated to the given bound node
  template<class Node, class It>
  It indices (const Node& node, It it) const
  {
    static_assert(std::is_same_v<typename std::iterator_traits<It>::iterator_category, std::random_access_iterator_tag>);

    auto first = it;
    auto last = Base::indices(node, it);
    assert(last - first == node.size());

    if (fe_.basix().dof_transformations_are_permutations())
    {
      // TODO: encodes bitwise (from right to left) whether the i'th edge is flipped in the real element.
      // cold be extracted from the nedelec basis, for example.
      std::uint32_t cell_info = 0;
      fe_.basix().permute(std::span{&*first, node.size()}, cell_info);
    }

    return last;
  }

protected:
  FiniteElement fe_;
};


template <class GV, int dimRange>
class BasixPreBasis<GV,dimRange>::Node
  : public LeafBasisNode
{
  static constexpr int dim = GV::dimension;

public:
  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = typename BasixPreBasis<GV,dimRange>::FiniteElement;

  //! Constructor; stores a pointer to the passed local finite-element `fe`.
  explicit Node (const FiniteElement& fe)
    : fe_{&fe}
    , element_{nullptr}
  {}

  //! Return current element; might raise an error if unbound
  const Element& element () const
  {
    assert(!!element_);
    return *element_;
  }

  /**
   * \brief Return the LocalFiniteElement for the element we are bound to; might raise an error if unbound.
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement () const
  {
    assert(!!fe_);
    return *fe_;
  }

  //! Bind to element. Stores a pointer to the passed element reference.
  void bind (const Element& e)
  {
    element_ = &e;
    fe_->bind(element_->geometry(), 0);
    this->setSize(fe_->size());
  }

protected:
  const FiniteElement* fe_;
  const Element* element_;
};


namespace BasisFactory {

/**
  * \brief A factory that can create a Basix pre-basis
  *
  * \ingroup FunctionSpaceBasesImplementations
  *
  * \tparam Basix    The type of the basix finite-element
  * \tparam dimRange Dimension of the range of the basis functions
  */
template<class Basix, int dimRange = 1>
auto basix (Basix basix)
{
  return [basix=std::move(basix)](const auto& gridView) {
    return BasixPreBasis<std::decay_t<decltype(gridView)>, dimRange>(gridView, std::move(basix));
  };
}

} // end namespace BasisFactory
} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASIXBASIS_HH
