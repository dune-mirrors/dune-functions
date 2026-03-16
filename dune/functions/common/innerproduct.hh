// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_COMMON_INNERPRODUCT_HH
#define DUNE_FUNCTIONS_COMMON_INNERPRODUCT_HH

#include <complex>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/promotiontraits.hh>

namespace Dune::Functions {
namespace Impl {

template <class T>
  requires requires(T z) { conj(z); }
auto conj(T z) { return conj(z); }

template <class T>
auto conj(T x) { return x; }

} // end namespace Impl

//! Returns the inner product a^T * b
template <class T0, int n, class T1>
auto innerProduct(FieldVector<T0,n> const& a, FieldVector<T1,n> const& b)
{
  using T = typename Dune::PromotionTraits<T0,T1>::PromotedType;

  T c = 0;
  for (int i = 0; i < n; ++i)
    c += Impl::conj(a[i]) * b[i];
  return c;
}

//! Returns the matrix inner product A^T : B
template <class T0, int n0, int n1, class T1>
auto innerProduct(FieldMatrix<T0,n0,n1> const& A, FieldMatrix<T1,n0,n1> const& B)
{
  using T = typename Dune::PromotionTraits<T0,T1>::PromotedType;

  T c = 0;
  for (int i0 = 0; i0 < n0; ++i0)
    for (int i1 = 0; i1 < n1; ++i1)
      c += Impl::conj(A[i0][i1]) * B[i0][i1];
  return c;
}

} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_COMMON_INNERPRODUCT_HH
