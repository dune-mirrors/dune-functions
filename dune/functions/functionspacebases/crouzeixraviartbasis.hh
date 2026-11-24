// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CROUZEIXRAVIARTBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CROUZEIXRAVIARTBASIS_HH

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/rannacherturekbasis.hh>


namespace Dune::Functions
{

namespace BasisFactory {

/**
 * \brief Create a pre-basis factory that can create a Crouzeix-Raviart pre-basis
 *
 * The Crouzeix-Raviart basis is only defined for simplex elements.
 * A Rannacher-Turek basis is used if the grid contains quadrilateral elements.
 *
 * \ingroup FunctionSpaceBasesImplementations
 */
template<class Dummy=void>
auto crouzeixRaviart()
{
  return [](const auto& gridView) {
    // The RannacherTurekPreBasis implements Crouzeix-Raviart FE if the grid is simplicial.
    return RannacherTurekPreBasis<std::decay_t<decltype(gridView)>>(gridView);
  };
}

} // end namespace BasisFactory




/** \brief Crouzeix-Raviart basis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * The Crouzeix-Raviart basis is only defined for simplex elements.
 * A Rannacher-Turek basis is used if the grid contains quadrilateral elements.
 * The implementation is therefore identical to the RannacherTurekBasis
 * implementation.
 *
 * \tparam GV The GridView that the space is defined on
 */
template<typename GV>
using CrouzeixRaviartBasis = DefaultGlobalBasis<RannacherTurekPreBasis<GV> >;

} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CROUZEIXRAVIARTBASIS_HH
