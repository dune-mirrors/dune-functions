// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASIXBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASIXBASIS_HH

#include <utility>

#include <basix/cell.h>

#include <basix/e-brezzi-douglas-marini.h>
#include <basix/e-crouzeix-raviart.h>
#include <basix/e-lagrange.h>
#include <basix/e-nedelec.h>
#include <basix/e-raviart-thomas.h>

#if 0
#include <basix/e-hhj.h>
#include <basix/e-regge.h>
#endif

#include <basix/element-families.h>

#include <dune/functions/functionspacebases/basix/prebasis.hh>

namespace Dune::Functions::BasisFactory {

/**
  * \brief A factory that can create a Basix pre-basis
  *
  * \ingroup FunctionSpaceBasesImplementations
  *
  * \tparam rangeClass Class of the basis range, e.g. one of {RangeClass::scalar,
  *         RangeClass::vector, RangeClass::matrix}
  * \tparam BasixType  The type of the basix finite-element
  */
template<RangeClass rangeClass, class BasixType>
auto basix (BasixType basix)
{
  return [basix=std::move(basix)]<class GridView>(const GridView& gridView) {
    auto factory = [basix=std::move(basix)](::basix::cell::type cell_type) {
      assert(basix.cell_type() == cell_type);
      return basix;
    };
    return BasixPreBasis<GridView, rangeClass, decltype(factory)>(gridView, std::move(factory));
  };
}

/** \brief A factory to create a Lagrange Basix pre-basis. */
template<class F = double>
auto basix_lagrange (int degree,
  ::basix::element::lagrange_variant variant = ::basix::element::lagrange_variant::equispaced)
{
  return [degree,variant]<class GridView>(const GridView& gridView) {
    auto factory = [degree,variant](::basix::cell::type cell_type) {
      return ::basix::element::create_lagrange<F>(cell_type, degree, variant, false);
    };
    return BasixPreBasis<GridView, RangeClass::scalar, decltype(factory)>(gridView, std::move(factory));
  };
}

/** \brief A factory to create a discontinuous Lagrange Basix pre-basis. */
template<class F = double>
auto basix_lagrangedg (int degree,
  ::basix::element::lagrange_variant variant = ::basix::element::lagrange_variant::equispaced)
{
  return [degree,variant]<class GridView>(const GridView& gridView) {
    auto factory = [degree,variant](::basix::cell::type cell_type) {
      return ::basix::element::create_lagrange<F>(cell_type, degree, variant, true);
    };
    return BasixPreBasis<GridView, RangeClass::scalar, decltype(factory)>(gridView, std::move(factory));
  };
}

/** \brief A factory to create a Nedelect Basix pre-basis of first kind. */
template<class F = double>
auto basix_nedelec (int degree,
  ::basix::element::lagrange_variant variant = ::basix::element::lagrange_variant::equispaced)
{
  return [degree,variant]<class GridView>(const GridView& gridView) {
    auto factory = [degree,variant](::basix::cell::type cell_type) {
      return ::basix::element::create_nedelec<F>(cell_type, degree, variant, false);
    };
    return BasixPreBasis<GridView, RangeClass::vector, decltype(factory)>(gridView, std::move(factory));
  };
}

/** \brief A factory to create a Raviart-Thomas Basix pre-basis */
template<class F = double>
auto basix_rt (int degree,
  ::basix::element::lagrange_variant variant = ::basix::element::lagrange_variant::equispaced)
{
  return [degree,variant]<class GridView>(const GridView& gridView) {
    auto factory = [degree,variant](::basix::cell::type cell_type) {
      return ::basix::element::create_rt<F>(cell_type, degree, variant, false);
    };
    return BasixPreBasis<GridView, RangeClass::vector, decltype(factory)>(gridView, std::move(factory));
  };
}

/** \brief A factory to create a Crouzeix-Raviart Basix pre-basis */
template<class F = double>
auto basix_cr ()
{
  return []<class GridView>(const GridView& gridView) {
    auto factory = [](::basix::cell::type cell_type) {
      return ::basix::element::create_cr<F>(cell_type, 1, false);
    };
    return BasixPreBasis<GridView, RangeClass::vector, decltype(factory)>(gridView, std::move(factory));
  };
}

/** \brief A factory to create a Brezzi-Douglas-Marini Basix pre-basis */
template<class F = double>
auto basix_bdm (int degree,
  ::basix::element::lagrange_variant variant = ::basix::element::lagrange_variant::equispaced)
{
  return [degree,variant]<class GridView>(const GridView& gridView) {
    auto factory = [degree,variant](::basix::cell::type cell_type) {
      return ::basix::element::create_bdm<F>(cell_type, degree, variant, false);
    };
    return BasixPreBasis<GridView, RangeClass::vector, decltype(factory)>(gridView, std::move(factory));
  };
}


#if 0 // not yet supported since double-piola transforms are missing

/**
 * \brief A factory to create a Regge Basix pre-basis.
 * \note This element is restricted to simplex elements.
 */
template<class F = double>
auto basix_regge (int degree)
{
  return [degree]<class GridView>(const GridView& gridView) {
    auto factory = [degree](::basix::cell::type cell_type) {
      return ::basix::element::create_regge<F>(cell_type, degree, false);
    };
    return BasixPreBasis<GridView, RangeClass::matrix, decltype(factory)>(gridView, std::move(factory));
  };
}

/**
 * \brief A factory to create a Hellan-Herrmann-Johnson Basix pre-basis.
 * \note This element is restricted to triangles.
 */
template<class F = double>
auto basix_hhj (int degree)
{
  return [degree]<class GridView>(const GridView& gridView) {
    auto factory = [degree](::basix::cell::type cell_type) {
      return ::basix::element::create_hhj<F>(cell_type, degree, false);
    };
    return BasixPreBasis<GridView, RangeClass::matrix, decltype(factory)>(gridView, std::move(factory));
  };
}
#endif

} // end namespace Dune::Functions::BasisFactory

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASIXBASIS_HH
