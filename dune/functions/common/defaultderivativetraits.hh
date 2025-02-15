// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_COMMON_DEFAULT_DERIVATIVE_TRAITS_HH
#define DUNE_FUNCTIONS_COMMON_DEFAULT_DERIVATIVE_TRAITS_HH

#include <type_traits>
#include <utility>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/promotiontraits.hh>
#include <dune/common/typetraits.hh>

namespace Dune {
namespace Functions {



/**
 * \brief Dummy range class to be used if no proper type is available
 *
 * \ingroup FunctionUtility
 */
class InvalidRange
{};


/**
 * \brief Default implementation for derivative traits
 *
 * \ingroup FunctionUtility
 *
 * This class provides sensible defaults for the range
 * of derivatives of functions with some common \p Domain
 * and \p Range types.
 */
template<class Signature>
struct DefaultDerivativeTraits
{
  //! Range of derivative for function with given signature
  typedef InvalidRange Range;
};


/**
 * \brief Default implementation for derivative traits
 *
 * \ingroup FunctionUtility
 *
 * Specialization for Signature = double(double)
 */
template<typename K1, typename K2>
  requires (Dune::IsNumber<K1>::value && Dune::IsNumber<K2>::value)
struct DefaultDerivativeTraits< K1(K2) >
{
  //! \copydoc DefaultDerivativeTraits::Range
  typedef typename Dune::PromotionTraits<K1, K2>::PromotedType K;
  typedef K Range;
};

/**
 * \brief Default implementation for derivative traits
 *
 * \ingroup FunctionUtility
 *
 * \tparam K Scalar range type
 *
 * Specialization for Signature = K(FieldVector<K,n>)
 */
template<typename K1, typename K2, int n>
  requires (Dune::IsNumber<K1>::value && Dune::IsNumber<K2>::value)
struct DefaultDerivativeTraits<K1(FieldVector<K2,n>)>
{
  //! \copydoc DefaultDerivativeTraits::Range
  typedef typename Dune::PromotionTraits<K1, K2>::PromotedType K;
  typedef FieldVector<K,n> Range;
};

/**
 * \brief Default implementation for derivative traits
 *
 * \ingroup FunctionUtility
 *
 * \tparam K Scalar range type
 *
 * Specialization for Signature = FieldVector<K,m>(FieldVector<K,n>)
 */
template<typename K1, typename K2, int n, int m>
  requires (Dune::IsNumber<K1>::value && Dune::IsNumber<K2>::value)
struct DefaultDerivativeTraits<FieldVector<K1,m>(FieldVector<K2,n>)>
{
  //! \copydoc DefaultDerivativeTraits::Range
  typedef typename Dune::PromotionTraits<K1, K2>::PromotedType K;
  typedef FieldMatrix<K,m,n> Range;
};

/**
 * \brief Default implementation for derivative traits
 *
 * \ingroup FunctionUtility
 *
 * \tparam K Scalar range type
 *
 * Specialization for Signature = FieldMatrix<K,1,m>(FieldVector<K,n>)
 */
template<typename K1, typename K2, int n, int m>
  requires (Dune::IsNumber<K1>::value && Dune::IsNumber<K2>::value)
struct DefaultDerivativeTraits<FieldMatrix<K1,1,m>(FieldVector<K2,n>)>
{
  //! \copydoc DefaultDerivativeTraits::Range
  typedef typename Dune::PromotionTraits<K1, K2>::PromotedType K;
  typedef FieldMatrix<K,m,n> Range;
};


}} // namespace Dune::Functions


#endif // DUNE_FUNCTIONS_COMMON_DEFAULT_DERIVATIVE_TRAITS_HH
