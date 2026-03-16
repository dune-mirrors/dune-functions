// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_COMMON_MULTIDOT_HH
#define DUNE_FUNCTIONS_COMMON_MULTIDOT_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/promotiontraits.hh>

namespace Dune::Functions {

//! Returns the product b0^T * A * b1, which results in a scalar
template <class T1, int n1, int n2, class T0, class T2>
auto multiDot(FieldMatrix<T1,n1,n2> const& A, FieldVector<T0,n1> const& b0, FieldVector<T2,n2> const& b1)
{
  using T = typename Dune::PromotionTraits<typename Dune::PromotionTraits<T0,T1>::PromotedType, T2>::PromotedType;

  T c = 0;
  for (int i1 = 0; i1 < n1; ++i1)
    for (int i2 = 0; i2 < n2; ++i2)
      c += b0[i1] * A[i1][i2] * b1[i2];
  return c;
}

//! Returns the product B0^T * A * B1, which results in a matrix of dimension n0 x n3
template <class T1, int n1, int n2, class T0, int n0, class T2, int n3>
auto multiDot(FieldMatrix<T1,n1,n2> const& A, FieldMatrix<T0,n0,n1> const& B0, FieldMatrix<T2,n2,n3> const& B1)
{
  using T = typename Dune::PromotionTraits<typename Dune::PromotionTraits<T0,T1>::PromotedType, T2>::PromotedType;

  FieldMatrix<T,n0,n3> C{};
  for (int i0 = 0; i0 < n0; ++i0)
    for (int i1 = 0; i1 < n1; ++i1)
      for (int i2 = 0; i2 < n2; ++i2)
        for (int i3 = 0; i3 < n3; ++i3)
          C[i0][i3] += B0[i0][i1] * A[i1][i2] * B1[i2][i3];
  return C;
}

} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_COMMON_MULTIDOT_HH
