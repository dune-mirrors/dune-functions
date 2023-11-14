// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LFEPREBASISMIXIN_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LFEPREBASISMIXIN_HH

#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/nodes.hh>

namespace Dune::Functions {

template<typename GV, class LFE>
class LFENode;

template<typename GV, class LFE>
class LFEPreBasisMixin;


/**
 * \brief A pre-basis mixin class parametrized with a local finite-element and a DOF layout.
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * This mixin class allows for simple construction of leaf pre-bases that are based on
 * a local finite-element and a DOF layout only. Examples are the refined Lagrange pre-bases,
 * or a hierarchical Lagrange pre-basis. Note that the layout is currently not capable of
 * describing a reordering of local DOFs if there are multiple assigned to a grid entity.
 * Thus higher-order continuous finite-elements are currently not possible to describe by
 * this mixin class.
 *
 * \b Example
 * \code{.cpp}
   template <class GV, class R = double>
   class RefinedP0PreBasis :
      public LFEPreBasisMixin<GV, RefinedP0LocalFiniteElement<typename GV::ctype,R,GV::dimension>>
   {
     using LFE = RefinedP0LocalFiniteElement<typename GV::ctype,R,GV::dimension>;
     using Base = LFEPreBasisMixin<GV, LFE>;
     static const int dim = GV::dimension;
   public:
     RefinedP0PreBasis (const GV& gv) :
       Base(gv, [](GeometryType gt, int) { return (gt.dim()==dim) ? (1<<dim) : 0; })
     {}
   };
 * \endcode
 *
 * \tparam GV   The grid view that the FE basis is defined on
 * \tparam LFE  The local finite-element type
 */
template <class GV, class LFE>
class LFEPreBasisMixin :
  public LeafPreBasisMapperMixIn< GV >
{
  using Base = LeafPreBasisMapperMixIn< GV >;

public:

  //! The grid view that the FE basis is defined on
  using GridView = GV;

  //! Type of the rtree node
  using Node = LFENode<GV, LFE>;

  /**
   * \brief Constructor for a given grid view object and layout.
   *
   * Requires that the local-finite element is default constructible.
   */
  template <class LFE_ = LFE,
    std::enable_if_t<std::is_default_constructible_v<LFE_>, int> = 0>
  LFEPreBasisMixin (const GridView& gv, MCMGLayout layout)
    : Base{gv, layout}
    , lfe_{}
  {}

  /**
   * \brief Constructor for a given grid view object, local finite-element and layout.
   *
   * Requires that the local-finite element is copyable or movable.
   */
  template <class LFE_>
  LFEPreBasisMixin (const GridView& gv, LFE_&& lfe, MCMGLayout layout)
    : Base{gv, layout}
    , lfe_{std::forward<LFE_>(lfe)}
  {}

  //! Create tree node
  Node makeNode () const
  {
    return Node{lfe_};
  }

private:
  LFE lfe_;
};

template <class GV, class LFE>
LFEPreBasisMixin(const GV&, const LFE&, MCMGLayout)
  -> LFEPreBasisMixin<GV,LFE>;




template <class GV, class LFE>
class LFENode
  : public LeafBasisNode
{
  static constexpr int dim = GV::dimension;

public:
  using size_type = std::size_t;
  using Element = typename GV::template Codim<0>::Entity;
  using FiniteElement = LFE;

  explicit LFENode (const LFE& lfe)
    : lfe_{&lfe}
    , element_{nullptr}
  {}

  //! Return current element, throw if unbound
  const Element& element () const
  {
    return *element_;
  }

  /** \brief Return the LocalFiniteElement for the element we are bound to
   *
   * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions module
   */
  const FiniteElement& finiteElement () const
  {
    return *lfe_;
  }

  //! Bind to element.
  void bind (const Element& e)
  {
    element_ = &e;
    this->setSize(lfe_->size());
  }

protected:
  const FiniteElement* lfe_;
  const Element* element_;
};


} // end namespace Dune::Functions


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LFEPREBASISMIXIN_HH
