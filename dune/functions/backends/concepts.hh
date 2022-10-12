// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_BACKEND_CONCEPTS_HH
#define DUNE_FUNCTIONS_BACKEND_CONCEPTS_HH


#include <utility>

#include <dune/common/concept.hh>

namespace Dune {
namespace Functions {
namespace Concept {

using namespace Dune::Concept;


// Concept for a ConstVectorBackend
template<class V, class GlobalBasis>
concept ConstVectorBackend = requires(const V& v, typename GlobalBasis::MultiIndex mi) {
  v[mi];
};

// Concept for a VectorBackend
template<class V, class GlobalBasis>
concept VectorBackend = ConstVectorBackend<V,GlobalBasis>
&& requires(V& v, typename GlobalBasis::MultiIndex mi, const GlobalBasis& sizeProvider) {
  v.resize(sizeProvider);
  v[mi] = v[mi];
};

} // namespace Dune::Functions::Concept
} // namespace Dune::Functions
} // namespace Dune


#endif // DUNE_FUNCTIONS_BACKEND_CONCEPTS_HH
