// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_PYTHON_FUNCTIONS_INTERPOLATE_HH
#define DUNE_PYTHON_FUNCTIONS_INTERPOLATE_HH

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/backends/istlvectorbackend.hh>

#include <dune/python/pybind11/numpy.h>
#include <dune/python/pybind11/pybind11.h>

namespace Dune {
namespace Python {
namespace Functions {

template<class Basis, typename T, typename R=T>
void interpolate(const Basis& basis, pybind11::array_t<T> x, const std::function<R(typename Basis::GridView::template Codim<0>::Geometry::GlobalCoordinate)>& f)
{
  if (x.shape(0) != basis.size())
    throw std::runtime_error("Coefficient vector has wrong size");

  auto x1 = x.template mutable_unchecked<1>();

  auto x2 = Dune::Functions::istlVectorBackend(x1);

  // TODO: The present code does not work if interpolation requires evaluating derivatives of f.
  // As an intermediate work-around, simply skip the interpolation in that case,
  // and report an error.
  constexpr bool hasHermiteDofs = std::is_same_v<Basis,Dune::Functions::DefaultGlobalBasis<Dune::Functions::ArgyrisPreBasis<typename Basis::GridView,double> > >;

  if constexpr (!hasHermiteDofs)
    interpolate(basis, x2, f);
  else
    DUNE_THROW(NotImplemented, "Method 'interpolate' not implemented for bases with Hermite degrees of freedom");
}

} /* namespace Functions */
} /* namespace Python */
} /* namespace Dune */

#endif
