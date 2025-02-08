// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICNODETORANGEMAP_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICNODETORANGEMAP_HH


#include <utility>
#include <type_traits>

#include <dune/common/concept.hh>

#include <dune/functions/functionspacebases/concepts.hh>
#include <dune/functions/common/indexaccess.hh>

namespace Dune {
namespace Functions {



/**
 * \brief A simple node to range map using the nested tree indices
 *
 * This map directly usses the tree path entries of the given
 * node to access the nested container.
 *
 * If the container does not provide any operator[] access,
 * it is simply forwarded for all nodes.
 */
struct HierarchicNodeToRangeMap
{
  template<class Node, class TreePath, class Range>
  decltype(auto) operator()(const Node&, [[maybe_unused]] const TreePath& treePath, Range&& y) const
  {
    if constexpr(Concept::HasIndexAccess<Range, Dune::index_constant<0>>)
      return resolveStaticMultiIndex(y, treePath);
    else
      return std::forward<Range>(y);
  }
};



} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HIERARCHICNODETORANGEMAP_HH
