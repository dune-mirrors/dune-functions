// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_PRECOMPUTEIDENTITY_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_PRECOMPUTEIDENTITY_HH

#include <cstddef>

namespace Dune::Functions {

/**
 * \brief Opaque identity of data affecting reference-basis precomputation.
 *
 * The object address distinguishes simultaneously existing reference finite
 * elements. The optional generation distinguishes successive contents stored
 * at the same address, for example after rebuilding an hp finite-element map.
 */
struct LocalBasisPrecomputeIdentity
{
  void const* object = nullptr;
  std::size_t generation = 0;

  bool operator==(LocalBasisPrecomputeIdentity const&) const = default;
};

} // namespace Dune::Functions

#endif
