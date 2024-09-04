// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_ARNOLDWINTHER_INC_HH
#define DUNE_LOCALFUNCTIONS_ARNOLDWINTHER_INC_HH
namespace Dune {
namespace Impl {
// generated with sympy from symfem library
template <class D, class R>
void ArnoldWintherLocalBasis<D, R>::evaluateFunction(
    const typename Traits::DomainType &in, std::vector<typename Traits::RangeType> &out) const {
  out.resize(size());
  auto iter = out.begin(); // generated with sympy from symfem library
  auto const &x = in[0], y = in[1];
  double diag;
  // 1th basis function
  diag = x * ((15.0 / 2.0) * x * y - 3 * y);

  *(iter++) = {{x * x * (3.0 / 2.0 - 5.0 / 2.0 * x) + y * (y * (18 - 10 * y) - 9) + 1, diag},
               {diag, x * y * ((15.0 / 2.0) * y - 6) + y * (3.0 / 2.0 - 3.0 / 2.0 * y)}};

  // 2th basis function
  diag =
      x * (x * (-10 * x + 90 * y + 18) + y * (45 * y - 42) - 9) + y * (y * (18 - 10 * y) - 9) + 1;

  *(iter++) = {{x * (x * (-30 * x - 45 * y + 21) + y * (30 * y - 36) + 9), diag},
               {diag, x * (30 * x * y + y * (90 * y - 108)) + y * (y * (15 * y - 51) + 36)}};

  // 3th basis function
  diag = x * ((45.0 / 2.0) * x * y + y * (15 * y - 15));

  *(iter++) = {{x * x * (-15.0 / 2.0 * x - 15 * y + 15.0 / 2.0), diag},
               {diag, x * (x * (18 - 10 * x) + y * ((45.0 / 2.0) * y - 18) - 9) +
                          y * (y * (5 * y - 27.0 / 2.0) + 15.0 / 2.0) + 1}};

  // 4th basis function
  diag = x * (-15.0 / 2.0 * x * y + 3 * y);

  *(iter++) = {{x * x * ((5.0 / 2.0) * x - 3.0 / 2.0), diag},
               {diag, x * y * (6 - 15.0 / 2.0 * y) + y * ((3.0 / 2.0) * y - 3.0 / 2.0)}};

  // 5th basis function
  diag = x * (x * (10 * x - 135.0 / 2.0 * y - 12) + y * (45 - 45 * y) + 3);

  *(iter++) = {{x * x * ((45.0 / 2.0) * x + 45 * y - 45.0 / 2.0), diag},
               {diag, x * (-30 * x * y + y * (78 - 135.0 / 2.0 * y)) +
                          y * (y * (81.0 / 2.0 - 15 * y) - 51.0 / 2.0)}};

  // 6th basis function
  diag = 0;

  *(iter++) = {{0, diag}, {diag, x * (x * (10 * x - 12) + 3)}};

  // 7th basis function
  diag = 0;

  *(iter++) = {{y * (y * (10 * y - 12) + 3), diag}, {diag, 0}};

  // 8th basis function
  diag = x * (-45.0 / 2.0 * x * y + 9 * y) + y * (y * (10 * y - 12) + 3);

  *(iter++) = {{x * (x * ((15.0 / 2.0) * x - 9.0 / 2.0) + y * (24 - 30 * y) - 3), diag},
               {diag, x * y * (18 - 45.0 / 2.0 * y) + y * ((9.0 / 2.0) * y - 9.0 / 2.0)}};

  // 9th basis function
  diag = x * (-45.0 / 2.0 * x * y + y * (15 - 15 * y));

  *(iter++) = {{x * x * ((15.0 / 2.0) * x + 15 * y - 15.0 / 2.0), diag},
               {diag, x * y * (18 - 45.0 / 2.0 * y) + y * (y * (27.0 / 2.0 - 5 * y) - 15.0 / 2.0)}};

  // 10th basis function
  diag = x * (-135 * x * y + y * (90 - 90 * y));

  *(iter++) = {
      {x * x * (45 * x + 90 * y - 45), diag},
      {diag, x * (x * (60 * x - 96) + y * (132 - 135 * y) + 36) + y * (y * (93 - 30 * y) - 63)}};

  // 11th basis function
  diag = x * (x * (-60 * x + 495 * y + 96) + y * (360 * y - 318) - 36);

  *(iter++) = {{x * (x * (-165 * x - 360 * y + 147) - 24 * y + 18), diag},
               {diag, x * (180 * x * y + y * (495 * y - 588)) + y * (y * (120 * y - 327) + 207)}};

  // 12th basis function
  diag = x * (270 * x * y + y * (180 * y - 180));

  *(iter++) = {{x * x * (-90 * x - 180 * y + 90), diag},
               {diag, x * (x * (180 - 120 * x) + y * (270 * y - 264) - 60) +
                          y * (y * (60 * y - 174) + 114)}};

  // 13th basis function
  diag = x * (x * (120 * x - 990 * y - 180) + y * (660 - 720 * y) + 60);

  *(iter++) = {{x * (x * (330 * x + 720 * y - 306) + 24 * y - 24), diag},
               {diag, x * (-360 * x * y + y * (1152 - 990 * y)) + y * (y * (642 - 240 * y) - 402)}};

  // 14th basis function
  diag = x * (-45 * x * y + 18 * y);

  *(iter++) = {{x * (x * (15 * x + 3) + 24 * y - 18) + y * (y * (60 * y - 96) + 36), diag},
               {diag, x * y * (36 - 45 * y) + y * (9 * y - 9)}};

  // 15th basis function
  diag = x * (-45 * x * y + y * (90 * y - 42)) + y * (y * (60 * y - 96) + 36);

  *(iter++) = {{x * (x * (15 * x - 90 * y + 21) + y * (192 - 180 * y) - 36), diag},
               {diag, x * y * (60 - 45 * y) + y * (y * (30 * y - 21) - 9)}};

  // 16th basis function
  diag = x * (90 * x * y - 36 * y);

  *(iter++) = {{x * (x * (6 - 30 * x) - 48 * y + 24) + y * (y * (180 - 120 * y) - 60), diag},
               {diag, x * y * (90 * y - 72) + y * (18 - 18 * y)}};

  // 17th basis function
  diag = x * (90 * x * y + y * (60 - 180 * y)) + y * (y * (180 - 120 * y) - 60);

  *(iter++) = {{x * (x * (-30 * x + 180 * y - 30) + y * (360 * y - 360) + 60), diag},
               {diag, x * y * (90 * y - 96) + y * (y * (54 - 60 * y) + 6)}};

  // 18th basis function
  diag = x * (-45 * M_SQRT2 * x * y + y * (-45 * M_SQRT2 * y + 36 * M_SQRT2));

  *(iter++) = {{x * (x * (15 * M_SQRT2 * x + 45 * M_SQRT2 * y - 12 * M_SQRT2) - 3 * M_SQRT2), diag},
               {diag, x * y * (-45 * M_SQRT2 * y + 48 * M_SQRT2) +
                          y * (y * (-15 * M_SQRT2 * y + 36 * M_SQRT2) - 21 * M_SQRT2)}};

  // 19th basis function
  diag = x * (-90 * M_SQRT2 * x * y + y * (-45 * M_SQRT2 * y + 54 * M_SQRT2));

  *(iter++) = {{x * (x * (30 * M_SQRT2 * x + 45 * M_SQRT2 * y - 33 * M_SQRT2) + 3 * M_SQRT2), diag},
               {diag, x * y * (-90 * M_SQRT2 * y + 84 * M_SQRT2) +
                          y * (y * (-15 * M_SQRT2 * y + 45 * M_SQRT2) - 30 * M_SQRT2)}};

  // 20th basis function
  diag = x * (90 * M_SQRT2 * x * y + y * (90 * M_SQRT2 * y - 72 * M_SQRT2));

  *(iter++) = {
      {x * (x * (-30 * M_SQRT2 * x - 90 * M_SQRT2 * y + 30 * M_SQRT2) + 12 * M_SQRT2 * y), diag},
      {diag, x * y * (90 * M_SQRT2 * y - 84 * M_SQRT2) +
                 y * (y * (30 * M_SQRT2 * y - 66 * M_SQRT2) + 36 * M_SQRT2)}};

  // 21th basis function
  diag = x * (180 * M_SQRT2 * x * y + y * (90 * M_SQRT2 * y - 108 * M_SQRT2));

  *(iter++) = {
      {x * (x * (-60 * M_SQRT2 * x - 90 * M_SQRT2 * y + 60 * M_SQRT2) - 12 * M_SQRT2 * y), diag},
      {diag, x * y * (180 * M_SQRT2 * y - 156 * M_SQRT2) +
                 y * (y * (30 * M_SQRT2 * y - 84 * M_SQRT2) + 54 * M_SQRT2)}};

  // 22th basis function
  diag = 0;

  *(iter++) = {{x * (-24 * x - 24 * y + 24), diag}, {diag, 0}};

  // 23th basis function
  diag = 24 * x * y;

  *(iter++) = {{x * (-24 * x - 48 * y + 24), diag}, {diag, -48 * x * y + y * (24 - 24 * y)}};

  // 24th basis function
  diag = 0;

  *(iter++) = {{0, diag}, {diag, -24 * x * y + y * (24 - 24 * y)}};
}
// generated with sympy from symfem library
template <class D, class R>
void ArnoldWintherLocalBasis<D, R>::evaluateJacobian(
    const typename Traits::DomainType &in, std::vector<typename Traits::JacobianType> &out) const {
  out.resize(size());
  auto iter = out.begin(); // generated with sympy from symfem library
  auto const &x = in[0], y = in[1];
  FieldVector<double, 2> diag;
  // 1th basis function
  diag = {15 * x * y - 3 * y, x * ((15.0 / 2.0) * x - 3)};

  *(iter++) = {{{x * (3 - 15.0 / 2.0 * x), y * (36 - 30 * y) - 9}, diag},
               {diag, {y * ((15.0 / 2.0) * y - 6), x * (15 * y - 6) - 3 * y + 3.0 / 2.0}}};

  // 2th basis function
  diag = {x * (-30 * x + 180 * y + 36) + y * (45 * y - 42) - 9,
          x * (90 * x + 90 * y - 42) + y * (36 - 30 * y) - 9};

  *(iter++) = {
      {{x * (-90 * x - 90 * y + 42) + y * (30 * y - 36) + 9, x * (-45 * x + 60 * y - 36)}, diag},
      {diag,
       {60 * x * y + y * (90 * y - 108), x * (30 * x + 180 * y - 108) + y * (45 * y - 102) + 36}}};

  // 3th basis function
  diag = {45 * x * y + y * (15 * y - 15), x * ((45.0 / 2.0) * x + 30 * y - 15)};

  *(iter++) = {{{x * (-45.0 / 2.0 * x - 30 * y + 15), -15 * x * x}, diag},
               {diag,
                {x * (36 - 30 * x) + y * ((45.0 / 2.0) * y - 18) - 9,
                 x * (45 * y - 18) + y * (15 * y - 27) + 15.0 / 2.0}}};

  // 4th basis function
  diag = {-15 * x * y + 3 * y, x * (3 - 15.0 / 2.0 * x)};

  *(iter++) = {{{x * ((15.0 / 2.0) * x - 3), 0}, diag},
               {diag, {y * (6 - 15.0 / 2.0 * y), x * (6 - 15 * y) + 3 * y - 3.0 / 2.0}}};

  // 5th basis function
  diag = {x * (30 * x - 135 * y - 24) + y * (45 - 45 * y) + 3,
          x * (-135.0 / 2.0 * x - 90 * y + 45)};

  *(iter++) = {{{x * ((135.0 / 2.0) * x + 90 * y - 45), 45 * x * x}, diag},
               {diag,
                {-60 * x * y + y * (78 - 135.0 / 2.0 * y),
                 x * (-30 * x - 135 * y + 78) + y * (81 - 45 * y) - 51.0 / 2.0}}};

  // 6th basis function
  diag = {0, 0};

  *(iter++) = {{{0, 0}, diag}, {diag, {x * (30 * x - 24) + 3, 0}}};

  // 7th basis function
  diag = {0, 0};

  *(iter++) = {{{0, y * (30 * y - 24) + 3}, diag}, {diag, {0, 0}}};

  // 8th basis function
  diag = {-45 * x * y + 9 * y, x * (9 - 45.0 / 2.0 * x) + y * (30 * y - 24) + 3};

  *(iter++) = {{{x * ((45.0 / 2.0) * x - 9) + y * (24 - 30 * y) - 3, x * (24 - 60 * y)}, diag},
               {diag, {y * (18 - 45.0 / 2.0 * y), x * (18 - 45 * y) + 9 * y - 9.0 / 2.0}}};

  // 9th basis function
  diag = {-45 * x * y + y * (15 - 15 * y), x * (-45.0 / 2.0 * x - 30 * y + 15)};

  *(iter++) = {
      {{x * ((45.0 / 2.0) * x + 30 * y - 15), 15 * x * x}, diag},
      {diag, {y * (18 - 45.0 / 2.0 * y), x * (18 - 45 * y) + y * (27 - 15 * y) - 15.0 / 2.0}}};

  // 10th basis function
  diag = {-270 * x * y + y * (90 - 90 * y), x * (-135 * x - 180 * y + 90)};

  *(iter++) = {{{x * (135 * x + 180 * y - 90), 90 * x * x}, diag},
               {diag,
                {x * (180 * x - 192) + y * (132 - 135 * y) + 36,
                 x * (132 - 270 * y) + y * (186 - 90 * y) - 63}}};

  // 11th basis function
  diag = {x * (-180 * x + 990 * y + 192) + y * (360 * y - 318) - 36, x * (495 * x + 720 * y - 318)};

  *(iter++) = {{{x * (-495 * x - 720 * y + 294) - 24 * y + 18, x * (-360 * x - 24)}, diag},
               {diag,
                {360 * x * y + y * (495 * y - 588),
                 x * (180 * x + 990 * y - 588) + y * (360 * y - 654) + 207}}};

  // 12th basis function
  diag = {540 * x * y + y * (180 * y - 180), x * (270 * x + 360 * y - 180)};

  *(iter++) = {{{x * (-270 * x - 360 * y + 180), -180 * x * x}, diag},
               {diag,
                {x * (360 - 360 * x) + y * (270 * y - 264) - 60,
                 x * (540 * y - 264) + y * (180 * y - 348) + 114}}};

  // 13th basis function
  diag = {x * (360 * x - 1980 * y - 360) + y * (660 - 720 * y) + 60,
          x * (-990 * x - 1440 * y + 660)};

  *(iter++) = {{{x * (990 * x + 1440 * y - 612) + 24 * y - 24, x * (720 * x + 24)}, diag},
               {diag,
                {-720 * x * y + y * (1152 - 990 * y),
                 x * (-360 * x - 1980 * y + 1152) + y * (1284 - 720 * y) - 402}}};

  // 14th basis function
  diag = {-90 * x * y + 18 * y, x * (18 - 45 * x)};

  *(iter++) = {{{x * (45 * x + 6) + 24 * y - 18, 24 * x + y * (180 * y - 192) + 36}, diag},
               {diag, {y * (36 - 45 * y), x * (36 - 90 * y) + 18 * y - 9}}};

  // 15th basis function
  diag = {-90 * x * y + y * (90 * y - 42), x * (-45 * x + 180 * y - 42) + y * (180 * y - 192) + 36};

  *(iter++) = {
      {{x * (45 * x - 180 * y + 42) + y * (192 - 180 * y) - 36, x * (-90 * x - 360 * y + 192)},
       diag},
      {diag, {y * (60 - 45 * y), x * (60 - 90 * y) + y * (90 * y - 42) - 9}}};

  // 16th basis function
  diag = {180 * x * y - 36 * y, x * (90 * x - 36)};

  *(iter++) = {{{x * (12 - 90 * x) - 48 * y + 24, -48 * x + y * (360 - 360 * y) - 60}, diag},
               {diag, {y * (90 * y - 72), x * (180 * y - 72) - 36 * y + 18}}};

  // 17th basis function
  diag = {180 * x * y + y * (60 - 180 * y), x * (90 * x - 360 * y + 60) + y * (360 - 360 * y) - 60};

  *(iter++) = {
      {{x * (-90 * x + 360 * y - 60) + y * (360 * y - 360) + 60, x * (180 * x + 720 * y - 360)},
       diag},
      {diag, {y * (90 * y - 96), x * (180 * y - 96) + y * (108 - 180 * y) + 6}}};

  // 18th basis function
  diag = {-90 * M_SQRT2 * x * y + y * (-45 * M_SQRT2 * y + 36 * M_SQRT2),
          x * (-45 * M_SQRT2 * x - 90 * M_SQRT2 * y + 36 * M_SQRT2)};

  *(iter++) = {{{x * (45 * M_SQRT2 * x + 90 * M_SQRT2 * y - 24 * M_SQRT2) - 3 * M_SQRT2,
                 45 * M_SQRT2 * x * x},
                diag},
               {diag,
                {y * (-45 * M_SQRT2 * y + 48 * M_SQRT2),
                 x * (-90 * M_SQRT2 * y + 48 * M_SQRT2) + y * (-45 * M_SQRT2 * y + 72 * M_SQRT2) -
                     21 * M_SQRT2}}};

  // 19th basis function
  diag = {-180 * M_SQRT2 * x * y + y * (-45 * M_SQRT2 * y + 54 * M_SQRT2),
          x * (-90 * M_SQRT2 * x - 90 * M_SQRT2 * y + 54 * M_SQRT2)};

  *(iter++) = {{{x * (90 * M_SQRT2 * x + 90 * M_SQRT2 * y - 66 * M_SQRT2) + 3 * M_SQRT2,
                 45 * M_SQRT2 * x * x},
                diag},
               {diag,
                {y * (-90 * M_SQRT2 * y + 84 * M_SQRT2),
                 x * (-180 * M_SQRT2 * y + 84 * M_SQRT2) + y * (-45 * M_SQRT2 * y + 90 * M_SQRT2) -
                     30 * M_SQRT2}}};

  // 20th basis function
  diag = {180 * M_SQRT2 * x * y + y * (90 * M_SQRT2 * y - 72 * M_SQRT2),
          x * (90 * M_SQRT2 * x + 180 * M_SQRT2 * y - 72 * M_SQRT2)};

  *(iter++) = {{{x * (-90 * M_SQRT2 * x - 180 * M_SQRT2 * y + 60 * M_SQRT2) + 12 * M_SQRT2 * y,
                 x * (-90 * M_SQRT2 * x + 12 * M_SQRT2)},
                diag},
               {diag,
                {y * (90 * M_SQRT2 * y - 84 * M_SQRT2), x * (180 * M_SQRT2 * y - 84 * M_SQRT2) +
                                                            y * (90 * M_SQRT2 * y - 132 * M_SQRT2) +
                                                            36 * M_SQRT2}}};

  // 21th basis function
  diag = {360 * M_SQRT2 * x * y + y * (90 * M_SQRT2 * y - 108 * M_SQRT2),
          x * (180 * M_SQRT2 * x + 180 * M_SQRT2 * y - 108 * M_SQRT2)};

  *(iter++) = {{{x * (-180 * M_SQRT2 * x - 180 * M_SQRT2 * y + 120 * M_SQRT2) - 12 * M_SQRT2 * y,
                 x * (-90 * M_SQRT2 * x - 12 * M_SQRT2)},
                diag},
               {diag,
                {y * (180 * M_SQRT2 * y - 156 * M_SQRT2),
                 x * (360 * M_SQRT2 * y - 156 * M_SQRT2) + y * (90 * M_SQRT2 * y - 168 * M_SQRT2) +
                     54 * M_SQRT2}}};

  // 22th basis function
  diag = {0, 0};

  *(iter++) = {{{-48 * x - 24 * y + 24, -24 * x}, diag}, {diag, {0, 0}}};

  // 23th basis function
  diag = {24 * y, 24 * x};

  *(iter++) = {{{-48 * x - 48 * y + 24, -48 * x}, diag}, {diag, {-48 * y, -48 * x - 48 * y + 24}}};

  // 24th basis function
  diag = {0, 0};

  *(iter++) = {{{0, 0}, diag}, {diag, {-24 * y, -24 * x - 48 * y + 24}}};
}
// generated with sympy from symfem library
template <class D, class R>
void ArnoldWintherLocalBasis<D, R>::partial(const std::array<unsigned int, dim> &order,
                                            const typename Traits::DomainType &in,
                                            std::vector<typename Traits::RangeType> &out) const {
  out.resize(size());
  auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
  if (totalOrder == 0) {
    evaluateFunction(in, out);
  } else if (totalOrder == 1) {
    if (order[0] == 1) {
      auto iter = out.begin();
      // generated with sympy from symfem library
      auto const &x = in[0], y = in[1];
      double diag;
      // 1th basis function
      diag = 15 * x * y - 3 * y;

      *(iter++) = {{x * (3 - 15.0 / 2.0 * x), diag}, {diag, y * ((15.0 / 2.0) * y - 6)}};

      // 2th basis function
      diag = x * (-30 * x + 180 * y + 36) + y * (45 * y - 42) - 9;

      *(iter++) = {{x * (-90 * x - 90 * y + 42) + y * (30 * y - 36) + 9, diag},
                   {diag, 60 * x * y + y * (90 * y - 108)}};

      // 3th basis function
      diag = 45 * x * y + y * (15 * y - 15);

      *(iter++) = {{x * (-45.0 / 2.0 * x - 30 * y + 15), diag},
                   {diag, x * (36 - 30 * x) + y * ((45.0 / 2.0) * y - 18) - 9}};

      // 4th basis function
      diag = -15 * x * y + 3 * y;

      *(iter++) = {{x * ((15.0 / 2.0) * x - 3), diag}, {diag, y * (6 - 15.0 / 2.0 * y)}};

      // 5th basis function
      diag = x * (30 * x - 135 * y - 24) + y * (45 - 45 * y) + 3;

      *(iter++) = {{x * ((135.0 / 2.0) * x + 90 * y - 45), diag},
                   {diag, -60 * x * y + y * (78 - 135.0 / 2.0 * y)}};

      // 6th basis function
      diag = 0;

      *(iter++) = {{0, diag}, {diag, x * (30 * x - 24) + 3}};

      // 7th basis function
      diag = 0;

      *(iter++) = {{0, diag}, {diag, 0}};

      // 8th basis function
      diag = -45 * x * y + 9 * y;

      *(iter++) = {{x * ((45.0 / 2.0) * x - 9) + y * (24 - 30 * y) - 3, diag},
                   {diag, y * (18 - 45.0 / 2.0 * y)}};

      // 9th basis function
      diag = -45 * x * y + y * (15 - 15 * y);

      *(iter++) = {{x * ((45.0 / 2.0) * x + 30 * y - 15), diag}, {diag, y * (18 - 45.0 / 2.0 * y)}};

      // 10th basis function
      diag = -270 * x * y + y * (90 - 90 * y);

      *(iter++) = {{x * (135 * x + 180 * y - 90), diag},
                   {diag, x * (180 * x - 192) + y * (132 - 135 * y) + 36}};

      // 11th basis function
      diag = x * (-180 * x + 990 * y + 192) + y * (360 * y - 318) - 36;

      *(iter++) = {{x * (-495 * x - 720 * y + 294) - 24 * y + 18, diag},
                   {diag, 360 * x * y + y * (495 * y - 588)}};

      // 12th basis function
      diag = 540 * x * y + y * (180 * y - 180);

      *(iter++) = {{x * (-270 * x - 360 * y + 180), diag},
                   {diag, x * (360 - 360 * x) + y * (270 * y - 264) - 60}};

      // 13th basis function
      diag = x * (360 * x - 1980 * y - 360) + y * (660 - 720 * y) + 60;

      *(iter++) = {{x * (990 * x + 1440 * y - 612) + 24 * y - 24, diag},
                   {diag, -720 * x * y + y * (1152 - 990 * y)}};

      // 14th basis function
      diag = -90 * x * y + 18 * y;

      *(iter++) = {{x * (45 * x + 6) + 24 * y - 18, diag}, {diag, y * (36 - 45 * y)}};

      // 15th basis function
      diag = -90 * x * y + y * (90 * y - 42);

      *(iter++) = {{x * (45 * x - 180 * y + 42) + y * (192 - 180 * y) - 36, diag},
                   {diag, y * (60 - 45 * y)}};

      // 16th basis function
      diag = 180 * x * y - 36 * y;

      *(iter++) = {{x * (12 - 90 * x) - 48 * y + 24, diag}, {diag, y * (90 * y - 72)}};

      // 17th basis function
      diag = 180 * x * y + y * (60 - 180 * y);

      *(iter++) = {{x * (-90 * x + 360 * y - 60) + y * (360 * y - 360) + 60, diag},
                   {diag, y * (90 * y - 96)}};

      // 18th basis function
      diag = -90 * M_SQRT2 * x * y + y * (-45 * M_SQRT2 * y + 36 * M_SQRT2);

      *(iter++) = {{x * (45 * M_SQRT2 * x + 90 * M_SQRT2 * y - 24 * M_SQRT2) - 3 * M_SQRT2, diag},
                   {diag, y * (-45 * M_SQRT2 * y + 48 * M_SQRT2)}};

      // 19th basis function
      diag = -180 * M_SQRT2 * x * y + y * (-45 * M_SQRT2 * y + 54 * M_SQRT2);

      *(iter++) = {{x * (90 * M_SQRT2 * x + 90 * M_SQRT2 * y - 66 * M_SQRT2) + 3 * M_SQRT2, diag},
                   {diag, y * (-90 * M_SQRT2 * y + 84 * M_SQRT2)}};

      // 20th basis function
      diag = 180 * M_SQRT2 * x * y + y * (90 * M_SQRT2 * y - 72 * M_SQRT2);

      *(iter++) = {
          {x * (-90 * M_SQRT2 * x - 180 * M_SQRT2 * y + 60 * M_SQRT2) + 12 * M_SQRT2 * y, diag},
          {diag, y * (90 * M_SQRT2 * y - 84 * M_SQRT2)}};

      // 21th basis function
      diag = 360 * M_SQRT2 * x * y + y * (90 * M_SQRT2 * y - 108 * M_SQRT2);

      *(iter++) = {
          {x * (-180 * M_SQRT2 * x - 180 * M_SQRT2 * y + 120 * M_SQRT2) - 12 * M_SQRT2 * y, diag},
          {diag, y * (180 * M_SQRT2 * y - 156 * M_SQRT2)}};

      // 22th basis function
      diag = 0;

      *(iter++) = {{-48 * x - 24 * y + 24, diag}, {diag, 0}};

      // 23th basis function
      diag = 24 * y;

      *(iter++) = {{-48 * x - 48 * y + 24, diag}, {diag, -48 * y}};

      // 24th basis function
      diag = 0;

      *(iter++) = {{0, diag}, {diag, -24 * y}};
    } else if (order[1] == 1) {
      auto iter = out.begin();
      // generated with sympy from symfem library
      auto const &x = in[0], y = in[1];
      double diag;
      // 1th basis function
      diag = x * ((15.0 / 2.0) * x - 3);

      *(iter++) = {{y * (36 - 30 * y) - 9, diag}, {diag, x * (15 * y - 6) - 3 * y + 3.0 / 2.0}};

      // 2th basis function
      diag = x * (90 * x + 90 * y - 42) + y * (36 - 30 * y) - 9;

      *(iter++) = {{x * (-45 * x + 60 * y - 36), diag},
                   {diag, x * (30 * x + 180 * y - 108) + y * (45 * y - 102) + 36}};

      // 3th basis function
      diag = x * ((45.0 / 2.0) * x + 30 * y - 15);

      *(iter++) = {{-15 * x * x, diag}, {diag, x * (45 * y - 18) + y * (15 * y - 27) + 15.0 / 2.0}};

      // 4th basis function
      diag = x * (3 - 15.0 / 2.0 * x);

      *(iter++) = {{0, diag}, {diag, x * (6 - 15 * y) + 3 * y - 3.0 / 2.0}};

      // 5th basis function
      diag = x * (-135.0 / 2.0 * x - 90 * y + 45);

      *(iter++) = {{45 * x * x, diag},
                   {diag, x * (-30 * x - 135 * y + 78) + y * (81 - 45 * y) - 51.0 / 2.0}};

      // 6th basis function
      diag = 0;

      *(iter++) = {{0, diag}, {diag, 0}};

      // 7th basis function
      diag = 0;

      *(iter++) = {{y * (30 * y - 24) + 3, diag}, {diag, 0}};

      // 8th basis function
      diag = x * (9 - 45.0 / 2.0 * x) + y * (30 * y - 24) + 3;

      *(iter++) = {{x * (24 - 60 * y), diag}, {diag, x * (18 - 45 * y) + 9 * y - 9.0 / 2.0}};

      // 9th basis function
      diag = x * (-45.0 / 2.0 * x - 30 * y + 15);

      *(iter++) = {{15 * x * x, diag}, {diag, x * (18 - 45 * y) + y * (27 - 15 * y) - 15.0 / 2.0}};

      // 10th basis function
      diag = x * (-135 * x - 180 * y + 90);

      *(iter++) = {{90 * x * x, diag}, {diag, x * (132 - 270 * y) + y * (186 - 90 * y) - 63}};

      // 11th basis function
      diag = x * (495 * x + 720 * y - 318);

      *(iter++) = {{x * (-360 * x - 24), diag},
                   {diag, x * (180 * x + 990 * y - 588) + y * (360 * y - 654) + 207}};

      // 12th basis function
      diag = x * (270 * x + 360 * y - 180);

      *(iter++) = {{-180 * x * x, diag}, {diag, x * (540 * y - 264) + y * (180 * y - 348) + 114}};

      // 13th basis function
      diag = x * (-990 * x - 1440 * y + 660);

      *(iter++) = {{x * (720 * x + 24), diag},
                   {diag, x * (-360 * x - 1980 * y + 1152) + y * (1284 - 720 * y) - 402}};

      // 14th basis function
      diag = x * (18 - 45 * x);

      *(iter++) = {{24 * x + y * (180 * y - 192) + 36, diag},
                   {diag, x * (36 - 90 * y) + 18 * y - 9}};

      // 15th basis function
      diag = x * (-45 * x + 180 * y - 42) + y * (180 * y - 192) + 36;

      *(iter++) = {{x * (-90 * x - 360 * y + 192), diag},
                   {diag, x * (60 - 90 * y) + y * (90 * y - 42) - 9}};

      // 16th basis function
      diag = x * (90 * x - 36);

      *(iter++) = {{-48 * x + y * (360 - 360 * y) - 60, diag},
                   {diag, x * (180 * y - 72) - 36 * y + 18}};

      // 17th basis function
      diag = x * (90 * x - 360 * y + 60) + y * (360 - 360 * y) - 60;

      *(iter++) = {{x * (180 * x + 720 * y - 360), diag},
                   {diag, x * (180 * y - 96) + y * (108 - 180 * y) + 6}};

      // 18th basis function
      diag = x * (-45 * M_SQRT2 * x - 90 * M_SQRT2 * y + 36 * M_SQRT2);

      *(iter++) = {{45 * M_SQRT2 * x * x, diag},
                   {diag, x * (-90 * M_SQRT2 * y + 48 * M_SQRT2) +
                              y * (-45 * M_SQRT2 * y + 72 * M_SQRT2) - 21 * M_SQRT2}};

      // 19th basis function
      diag = x * (-90 * M_SQRT2 * x - 90 * M_SQRT2 * y + 54 * M_SQRT2);

      *(iter++) = {{45 * M_SQRT2 * x * x, diag},
                   {diag, x * (-180 * M_SQRT2 * y + 84 * M_SQRT2) +
                              y * (-45 * M_SQRT2 * y + 90 * M_SQRT2) - 30 * M_SQRT2}};

      // 20th basis function
      diag = x * (90 * M_SQRT2 * x + 180 * M_SQRT2 * y - 72 * M_SQRT2);

      *(iter++) = {{x * (-90 * M_SQRT2 * x + 12 * M_SQRT2), diag},
                   {diag, x * (180 * M_SQRT2 * y - 84 * M_SQRT2) +
                              y * (90 * M_SQRT2 * y - 132 * M_SQRT2) + 36 * M_SQRT2}};

      // 21th basis function
      diag = x * (180 * M_SQRT2 * x + 180 * M_SQRT2 * y - 108 * M_SQRT2);

      *(iter++) = {{x * (-90 * M_SQRT2 * x - 12 * M_SQRT2), diag},
                   {diag, x * (360 * M_SQRT2 * y - 156 * M_SQRT2) +
                              y * (90 * M_SQRT2 * y - 168 * M_SQRT2) + 54 * M_SQRT2}};

      // 22th basis function
      diag = 0;

      *(iter++) = {{-24 * x, diag}, {diag, 0}};

      // 23th basis function
      diag = 24 * x;

      *(iter++) = {{-48 * x, diag}, {diag, -48 * x - 48 * y + 24}};

      // 24th basis function
      diag = 0;

      *(iter++) = {{0, diag}, {diag, -24 * x - 48 * y + 24}};
    } else
      DUNE_THROW(NotImplemented, "Invalid partial derivative");
  } else
    DUNE_THROW(NotImplemented, "Higher order derivatives are not implemented");
}

} // namespace Impl
} // namespace Dune
#endif
