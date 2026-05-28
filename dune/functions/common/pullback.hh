// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_COMMON_PULLBACK_HH
#define DUNE_FUNCTIONS_COMMON_PULLBACK_HH

#include <dune/common/diagonalmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/promotiontraits.hh>
#include <dune/common/rangeutilities.hh>

namespace Dune::Functions::Impl {

//! Returns the product b^T * A * b, which results in a scalar
template <class S, int n, class T>
auto pullback(FieldMatrix<S,n,n> const& A, FieldVector<T,n> const& b)
{
  using U = typename Dune::PromotionTraits<S,T>::PromotedType;

  U c = 0;
  for (int i1 = 0; i1 < n; ++i1)
    for (int i2 = 0; i2 < n; ++i2)
      c += b[i1] * A[i1][i2] * b[i2];
  return c;
}

//! Returns the product B^T * A * B, which results in a matrix of dimension m x m
template <class S, int n, class T, int m>
auto pullback(FieldMatrix<S,n,n> const& A, FieldMatrix<T,n,m> const& B)
{
  using U = typename Dune::PromotionTraits<S,T>::PromotedType;

  FieldMatrix<U,m,m> C{};
  for (int i0 = 0; i0 < m; ++i0)
    for (int i1 = 0; i1 < n; ++i1)
      for (int i2 = 0; i2 < n; ++i2)
        for (int i3 = 0; i3 < m; ++i3)
          C[i0][i3] += B[i1][i0] * A[i1][i2] * B[i2][i3];
  return C;
}

//! Returns the product B^T * A * B, which results in a matrix of dimension m x m
template <class MatrixA, class MatrixB>
auto pullback(MatrixA const& A, MatrixB const& B)
{
  using U = typename Dune::PromotionTraits<typename MatrixA::value_type, typename MatrixB::value_type>::PromotedType;

  constexpr int m = MatrixB::M();
  FieldMatrix<U,m,m> C{};
  for (auto i1 : Dune::range(A.N()))
    for (auto [a12,i2] : Dune::sparseRange(A[i1]))
      for (auto [b23,i3] : Dune::sparseRange(B[i2]))
        for (auto [b10,i0] : Dune::sparseRange(B[i1]))
          C[i0][i3] += transposed(b10) * a12 * b23;
  return C;
}

} // end namespace Dune::Functions::Impl

#endif // DUNE_FUNCTIONS_COMMON_PULLBACK_HH
