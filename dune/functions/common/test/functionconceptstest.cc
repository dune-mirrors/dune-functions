// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <dune/functions/common/functionconcepts.hh>


struct TestLocalFunction
{
  double operator()(double x) const
  {
    return x;
  }

  void bind(const int& context)
  {
    context_ = context;
    bound_ = true;
  }

  void unbind()
  {
    bound_ = false;
  }

  bool bound() const
  {
    return bound_;
  }

  const int& localContext() const
  {
    return context_;
  }

private:
  bool bound_ = false;
  int context_ = 0;
};

TestLocalFunction derivative(const TestLocalFunction&)
{
  return {};
}


int main()
{
  using Dune::Functions::Concept::isDifferentiableLocalFunction;
  using Dune::Functions::Concept::isLocalFunction;

  using Signature = double(double);
  using LocalContext = int;

  static_assert(isLocalFunction<TestLocalFunction, Signature, LocalContext>());
  static_assert(isLocalFunction<TestLocalFunction&, Signature, LocalContext>());
  static_assert(isLocalFunction<const TestLocalFunction, Signature, LocalContext>());
  static_assert(isLocalFunction<const TestLocalFunction&, Signature, LocalContext>());

  static_assert(isDifferentiableLocalFunction<TestLocalFunction, Signature, LocalContext>());
  static_assert(isDifferentiableLocalFunction<TestLocalFunction&, Signature, LocalContext>());
  static_assert(isDifferentiableLocalFunction<const TestLocalFunction, Signature, LocalContext>());
  static_assert(isDifferentiableLocalFunction<const TestLocalFunction&, Signature, LocalContext>());

  return 0;
}
