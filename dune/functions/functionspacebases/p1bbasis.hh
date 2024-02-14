// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIMPLEXP1B_BASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIMPLEXP1B_BASIS_HH

#include <cassert>
#include <type_traits>

#include <dune/functions/functionspacebases/lfeprebasismixin.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/bubble/simplexp1b.hh>

namespace Dune::Functions {

/**
  * \brief Piecewise linear Lagrange functions enriched with element bubble functions.
  *
  * The set of basis functions contains the classical linear Lagrange basis functions
  * and a single "bubble" function per element that vanishes on all facacets.
  *
  * A classical example where this kind of basis is used in the discretization
  * of the Stokes equation with the stable mixed-element called MINI element,
  * see
  *
  *   Arnold, D.N., Brezzi, F. and Fortin, M. A stable finite element for the
  *   Stokes equations. Calcolo 21, 337-344 (1984). doi: 10.1007/BF02576171
  *
  * The velocity field is discretized with continuous piecewise linear
  * functions enriched by a bubble function.
  *
  * \note The implementation here is restricted to simplex elements.
  *
  * \tparam GV  The type of the GridView the basis is defined on
  * \tparam R   The range field type of the local basis functions [double]
  **/
template <class GV, class R = double>
class P1BPreBasis :
  public LFEPreBasisMixin<GV, SimplexP1BLocalFiniteElement<typename GV::ctype,R,GV::dimension>>
{
  static constexpr int dim = GV::dimension;

  using LFE = SimplexP1BLocalFiniteElement<typename GV::ctype,R,dim>;
  using Base = LFEPreBasisMixin<GV, LFE>;

public:
  explicit P1BPreBasis (const GV& gv) :
    Base(gv, [](GeometryType gt, int) { return (gt.dim()==0) ? 1 : (gt.dim()==dim) ? 1 : 0; })
  {
    assert(gv.indexSet().types(0).size() == 1);
    assert(gv.indexSet().types(0).front() == GeometryTypes::simplex(dim));
  }
};

namespace BasisFactory {

/**
 * \ingroup FunctionSpaceBasesImplementations
 * \brief Create a pre-basis factory that can create a SimplexP1B pre-basis
 *
 * \tparam R  The range field type of the local basis functions [double]
 */
template <class R = double>
auto p1b ()
{
  return [](const auto& gridView) {
    return P1BPreBasis<std::decay_t<decltype(gridView)>, R>(gridView);
  };
}

} // end namespace BasisFactory
} // end namespace Dune::Functions


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_SIMPLEXP1B_BASIS_HH
