// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_FUNCTIONS_GENERICLOCALFINITEELEMENTS_POLYNOMIALBASISCOEFFICIENTS_HH
#define DUNE_FUNCTIONS_GENERICLOCALFINITEELEMENTS_POLYNOMIALBASISCOEFFICIENTS_HH

#include <array>
#include <cmath>
#include <string>

#include <dune/common/fmatrix.hh>
#include <dune/localfunctions/utility/field.hh>

namespace Dune
{
namespace Impl
{
/** Small Helper class to fill a FieldMatrix into a SparseCoefficientMatrix. */
template<class F, unsigned int Rows, unsigned int Cols = Rows>
struct PolynomialBasisMatrix {
    using Field = F;

    PolynomialBasisMatrix() = delete;

    PolynomialBasisMatrix(Dune::FieldMatrix<Field, Rows, Cols> &&mat) : mat_(mat) {}

    static unsigned int constexpr cols() { return Cols; };
    static unsigned int constexpr rows() { return Rows; };

    template<class Vector>
    void row(unsigned int const row, Vector &vec) const
    {
      unsigned int const N = cols();
      assert(vec.size() == N);

      for (unsigned int i = 0; i < N; ++i) {
        field_cast(mat_[row][i], vec[i]);
      }
    }

    Dune::FieldMatrix<Field, Rows, Cols> mat_;
};
template<class F, unsigned int dim>
using HermiteMatrix = PolynomialBasisMatrix<F, (dim == 1)   ? 4
                                               : (dim == 2) ? 10
                                                            : 20>;

template<class F>
using ArgyrisMatrix = PolynomialBasisMatrix<F, 21>;

template<class F>
using MorleyMatrix = PolynomialBasisMatrix<F, 6>;

namespace PolynomialBasisCoefficients
{
/**
 * @brief Get the Hermite Coefficients Matrix
 *
 * @tparam F Field type
 * @tparam dim dimesion of domain of Reference triangle
 * @return HermiteVecMatrix<F,dim> where size of the underlying matrix depends on dim
 */
template<class F, int dim, bool reduced>
constexpr auto getHermiteCoefficients()
{
  static_assert(dim > 0 and dim < 4 and not(reduced and dim != 2));
  if constexpr (dim == 1) {
    return HermiteMatrix<F, 1>({{1, 0, -3, 2}, {0, 1, -2, 1}, {0, 0, 3, -2}, {0, 0, -1, 1}});
  } else if constexpr (dim == 2) {
    if constexpr (reduced) {
      auto w = std::array<F, 9>{1. / 3,  1. / 18, 1. / 18, 1. / 3, -1. / 9,
                                1. / 18, 1. / 3,  1. / 18, -1. / 9};
      return PolynomialBasisMatrix<F, 9, 10>({
          {1, 0, 0, -3, -13 + w[0] * 27, -3, 2, 13 - w[0] * 27, 13 - w[0] * 27, 2},
          {0, 1, 0, -2, -3 + w[1] * 27, 0, 1, 3 - w[1] * 27, 2 - w[1] * 27, 0},
          {0, 0, 1, 0, -3 + w[2] * 27, -2, 0, 2 - w[2] * 27, 3 - w[2] * 27, 1},
          {0, 0, 0, 3, -7 + w[3] * 27, 0, -2, 7 - w[3] * 27, 7 - w[3] * 27, 0},
          {0, 0, 0, -1, 2 + w[4] * 27, 0, 1, -2 - w[4] * 27, -2 - w[4] * 27, 0},
          {0, 0, 0, 0, -1 + w[5] * 27, 0, 0, 2 - w[5] * 27, 1 - w[5] * 27, 0},
          {0, 0, 0, 0, -7 + w[6] * 27, 3, 0, 7 - w[6] * 27, 7 - w[6] * 27, -2},
          {0, 0, 0, 0, -1 + w[7] * 27, 0, 0, 1 - w[7] * 27, 2 - w[7] * 27, 0},
          {0, 0, 0, 0, 2 + w[8] * 27, -1, 0, -2 - w[8] * 27, -2 - w[8] * 27, 1},
      });
    } else
      return HermiteMatrix<F, 2>({{1, 0, 0, -3, -13, -3, 2, 13, 13, 2},
                                  {0, 1, 0, -2, -3, 0, 1, 3, 2, 0},
                                  {0, 0, 1, 0, -3, -2, 0, 2, 3, 1}, // l_2
                                  {0, 0, 0, 3, -7, 0, -2, 7, 7, 0},
                                  {0, 0, 0, -1, 2, 0, 1, -2, -2, 0},
                                  {0, 0, 0, 0, -1, 0, 0, 2, 1, 0},
                                  {0, 0, 0, 0, -7, 3, 0, 7, 7, -2}, // l_6
                                  {0, 0, 0, 0, -1, 0, 0, 1, 2, 0},
                                  {0, 0, 0, 0, 2, -1, 0, -2, -2, 1},
                                  {0, 0, 0, 0, 27, 0, 0, -27, -27, 0}}); // l_9, inner dof
  } else if constexpr (dim == 3) {
    return HermiteMatrix<F, 3>({{1, 0,  0,  0, -3, -13, -3, -13, -13, -3, // deg 0 to 2
                                 2, 13, 13, 2, 13, 33,  13, 13,  13,  2}, // deg 3
                                {0, 1, 0, 0, -2, -3, 0, -3, 0, 0, 1, 3, 2, 0, 3, 4, 0, 2, 0, 0},
                                {0, 0, 1, 0, 0, -3, -2, 0, -3, 0, 0, 2, 3, 1, 0, 4, 3, 0, 2, 0},
                                {0, 0, 0, 1, 0, 0, 0, -3, -3, -2, 0, 0, 0, 0, 2, 4, 2, 3, 3, 1},
                                {0,  0, 0, 0, 3, -7, 0, -7, 0, 0, // l_4
                                 -2, 7, 7, 0, 7, 7,  0, 7,  0, 0},
                                {0, 0, 0, 0, -1, 2, 0, 2, 0, 0, 1, -2, -2, 0, -2, -2, 0, -2, 0, 0},
                                {0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0},
                                {0, 0, 0, 0,  0, -7, 3, 0, -7, 0, // l_8
                                 0, 7, 7, -2, 0, 7,  7, 0, 7,  0},
                                {0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 2, -1, 0, 2, 0, 0, -2, -2, 1, 0, -2, -2, 0, -2, 0},
                                {0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0},
                                {0, 0, 0, 0, 0, 0, 0, -7, -7, 3, // l_12
                                 0, 0, 0, 0, 7, 7, 7, 7,  7,  -2},
                                {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0},
                                {0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0},
                                {0, 0, 0, 0, 0, 0, 0, 2, 2, -1, 0, 0, 0, 0, -2, -2, -2, -2, -2, 1},
                                // l_16, from here on inner dofs
                                {0, 0,   0,   0, 0, 27,  0, 0, 0, 0, // bottom
                                 0, -27, -27, 0, 0, -27, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0,   0,   0, 27,  0, 0, // front
                                 0, 0, 0, 0, -27, -27, 0, -27, 0, 0},
                                {0, 0, 0, 0, 0, 0,   0,   0, 27,  0, // left
                                 0, 0, 0, 0, 0, -27, -27, 0, -27, 0},
                                {0, 0, 0, 0, 0, 0,  0, 0, 0, 0, // right
                                 0, 0, 0, 0, 0, 27, 0, 0, 0, 0}});
  }
}

/**
 * @brief Get the Argyris Coefficient Matrix
 * values taken from https://defelement.com/elements/argyris.html.
 * @tparam F Field Type
 * @return ArgyrisVecMatrix<F>
 */
template<class F>
ArgyrisMatrix<F> getArgyrisCoefficients()
{
  F sqrt2 = -8. * std::sqrt(2.);
  return ArgyrisMatrix<F>({

      // vertex functionals
      // l_0
      {/*0th order*/ 1,
       /*1th order*/ 0,   0,
       /*2th order*/ 0,   0, 0,
       /*3th order*/ -10, 0, 0,   -10,
       /*4th order*/ 15,  0, -30, 0,   15,
       /*5th order*/ -6,  0, 30,  30,  0,  -6},
      // l_1
      {0, 1, 0, 0, 0, 0, -6, 0, -11, 0, 8, 0, 10, 18, 0, -3, 0, 1, -10, -8, 0},
      // l_2, l_1 mirrored
      {0, 0, 1, 0, 0, 0, 0, -11, 0, -6, 0, 18, 10, 0, 8, 0, -8, -10, 1, 0, -3},
      // l_3
      {0, 0, 0, 0.5, 0, 0, -1.5, 0, 0, 0, 1.5, 0, -1.5, 0, 0, -0.5, 0, 1.5, 1, 0, 0},
      // l_4
      {0, 0, 0, 0, 1, 0, 0, -4, -4, 0, 0, 5, 10, 5, 0, 0, -2, -6, -6, -2, 0},
      // l_5, l_3 mirrored
      {0, 0, 0, 0, 0, 0.5, 0, 0, 0, -1.5, 0, 0, -1.5, 0, 1.5, 0, 0, 1, 1.5, 0, -0.5},
      // l_6
      {0, 0, 0, 0, 0, 0, 10, 0, 0, 0, -15, 0, 15, 0, 0, 6, 0, -15, -15, 0, 0},
      // l_7
      {0, 0, 0, 0, 0, 0, -4, 0, 0, 0, 7, 0, -3.5, 0, 0, -3, 0, 3.5, 3.5, 0, 0},
      // l_8
      {0, 0, 0, 0, 0, 0, 0, -5, 0, 0, 0, 14, 18.5, 0, 0, 0, -8, -18.5, -13.5, 0, 0},
      // l_9
      {0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, -1, 0, 0.25, 0, 0, 0.5, 0, -0.25, -0.25, 0, 0},
      // l_10
      {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -3, -3.5, 0, 0, 0, 2, 3.5, 2.5, 0, 0},
      // l_11
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.25, 0, 0, 0, 0, -0.75, -1.25, 0, 0},
      // l_12 mirrors l_6
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 15, 0, -15, 0, 0, -15, -15, 0, 6},
      // l_13 mirrors l_8
      {0, 0, 0, 0, 0, 0, 0, 0, -5, 0, 0, 0, 18.5, 14, 0, 0, 0, -13.5, -18.5, -8, 0},
      // l_14 mirrors l_7
      {0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 0, 0, -3.5, 0, 7, 0, 0, 3.5, 3.5, 0, -3},
      // l_15 mirrors l_11
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.25, 0, 0, 0, 0, -1.25, -0.75, 0, 0},
      // l_16 mirrors l_10
      {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -3.5, -3, 0, 0, 0, 2.5, 3.5, 2, 0},
      // l_17 mirrors l_9
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0.25, 0, -1, 0, 0, -0.25, -0.25, 0, 0.5},
      // edge functionals
      // l_18
      {0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, -32, -32, 0, 0, 0, 16, 32, 16, 0, 0},
      // l_19
      {0, 0, 0, 0, 0, 0, 0, 0, -16, 0, 0, 0, 32, 32, 0, 0, 0, -16, -32, -16, 0},
      // l_20
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1. * sqrt2, 0, 0, 0, 0, sqrt2, sqrt2, 0, 0}});
}

/**
 * @brief Get the Morley Coefficient Matrix
 * values taken from https://defelement.com/elements/argyris.html.
 * @tparam F Field Type
 * @return MorleyVecMatrix<F>
 */
template<class F>
MorleyMatrix<F> getMorleyCoefficients()
{
  F sqrt2 = 0.5 * std::sqrt(2.);

  return MorleyMatrix<F>({{1, -1, -1, 0, 2, 0},
                          {0, 0.5, 0.5, 0.5, -1, -0.5},
                          {0, 0.5, 0.5, -0.5, -1, 0.5},
                          {0, 0, 1, 0, 0, -1},
                          {0, -1., 0, 1., 0, 0},
                          {0, sqrt2, sqrt2, -sqrt2, -2. * sqrt2, -sqrt2}});
}
} // namespace PolynomialBasisCoefficients
} // namespace Impl
} // namespace Dune
#endif
