// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <complex>
#include <type_traits>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/functions/common/defaultderivativetraits.hh>


int main ( int argc, char **argv )
{
  Dune::MPIHelper::instance(argc, argv);

  // check Number(Number) signature
  static_assert(std::is_same_v<typename Dune::Functions::DefaultDerivativeTraits<double(double)>::Range, double>);
  static_assert(std::is_same_v<typename Dune::Functions::DefaultDerivativeTraits<float(double)>::Range, double>);
  static_assert(std::is_same_v<typename Dune::Functions::DefaultDerivativeTraits<double(float)>::Range, double>);
  static_assert(std::is_same_v<typename Dune::Functions::DefaultDerivativeTraits<std::complex<double>(double)>::Range, std::complex<double>>);

  // check Number(Vector) signature
  static_assert(std::is_same_v<typename Dune::Functions::DefaultDerivativeTraits<double(Dune::FieldVector<double,2>)>::Range, Dune::FieldVector<double,2>>);
  static_assert(std::is_same_v<typename Dune::Functions::DefaultDerivativeTraits<float(Dune::FieldVector<double,2>)>::Range, Dune::FieldVector<double,2>>);
  static_assert(std::is_same_v<typename Dune::Functions::DefaultDerivativeTraits<double(Dune::FieldVector<float,2>)>::Range, Dune::FieldVector<double,2>>);
  static_assert(std::is_same_v<typename Dune::Functions::DefaultDerivativeTraits<std::complex<double>(Dune::FieldVector<double,2>)>::Range, Dune::FieldVector<std::complex<double>,2>>);

  // check Vector(Vector) signature
  static_assert(std::is_same_v<typename Dune::Functions::DefaultDerivativeTraits<Dune::FieldVector<double,3>(Dune::FieldVector<double,2>)>::Range, Dune::FieldMatrix<double,3,2>>);
  static_assert(std::is_same_v<typename Dune::Functions::DefaultDerivativeTraits<Dune::FieldVector<float,3>(Dune::FieldVector<double,2>)>::Range, Dune::FieldMatrix<double,3,2>>);
  static_assert(std::is_same_v<typename Dune::Functions::DefaultDerivativeTraits<Dune::FieldVector<double,3>(Dune::FieldVector<float,2>)>::Range, Dune::FieldMatrix<double,3,2>>);

  // check that the wrong dimensions are correctly detected
  static_assert(not std::is_same_v<typename Dune::Functions::DefaultDerivativeTraits<Dune::FieldVector<double,3>(Dune::FieldVector<double,2>)>::Range, Dune::FieldMatrix<double,2,3>>);
}