// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_DERIVATIVE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_DERIVATIVE_HH

#include <array>

namespace Dune::Functions::Derivatives {

/**
 * \brief Tag types for selecting values or derivatives of local shape functions.
 *
 * These tags form the experimental dispatch vocabulary used by transformed
 * local finite elements.  A tag object is passed to evaluate(), precompute(),
 * or finalize() to select the quantity to compute without overloading on the
 * output container type.
 *
 * The tags are intentionally lightweight value types.  Tags without data denote
 * a whole derivative quantity, while Partial stores the selected coordinate
 * direction.
 */

//! Select evaluation of shape-function values.
struct Value { bool operator==(Value const&) const = default; };

//! Select evaluation of the Jacobian/first derivative of shape functions.
struct Jacobian { bool operator==(Jacobian const&) const = default; };

//! Select evaluation of gradients.
struct Gradient { bool operator==(Gradient const&) const = default; };

//! Select evaluation of divergences.
struct Divergence { bool operator==(Divergence const&) const = default; };

//! Select evaluation of curls.
struct Curl { bool operator==(Curl const&) const = default; };

//! Select evaluation of Hessians.
struct Hessian { bool operator==(Hessian const&) const = default; };

//! Select evaluation of Laplacians.
struct Laplacian { bool operator==(Laplacian const&) const = default; };

//! Select evaluation of a partial derivative in coordinate direction i.
template <std::size_t dim>
struct Partial {
  std::array<unsigned int,dim> orders;
  bool operator==(Partial const&) const = default;
};

} // end namespace Dune::Functions::Derivatives

#endif
