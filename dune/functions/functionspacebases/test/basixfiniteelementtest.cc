// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basix/finiteelement.hh>
#include <dune/functions/functionspacebases/basix/localfiniteelement.hh>
#include <dune/functions/functionspacebases/basix/utility.hh>
#include <dune/localfunctions/test/test-localfe.hh>

#include <basix/finite-element.h>
#include <basix/e-lagrange.h>

int testUtility ()
{
  using namespace Dune;
  using namespace Dune::Functions;
  Dune::TestSuite test;

  test.check(Basix::geometryType(::basix::cell::type::point) == GeometryTypes::vertex);
  test.check(Basix::geometryType(::basix::cell::type::interval) == GeometryTypes::line);
  test.check(Basix::geometryType(::basix::cell::type::triangle) == GeometryTypes::triangle);
  test.check(Basix::geometryType(::basix::cell::type::quadrilateral) == GeometryTypes::quadrilateral);
  test.check(Basix::geometryType(::basix::cell::type::tetrahedron) == GeometryTypes::tetrahedron);
  test.check(Basix::geometryType(::basix::cell::type::hexahedron) == GeometryTypes::hexahedron);
  test.check(Basix::geometryType(::basix::cell::type::prism) == GeometryTypes::prism);
  test.check(Basix::geometryType(::basix::cell::type::pyramid) == GeometryTypes::pyramid);


  test.check(::basix::cell::type::point == Basix::cellType(GeometryTypes::vertex));
  test.check(::basix::cell::type::interval == Basix::cellType(GeometryTypes::line));
  test.check(::basix::cell::type::interval == Basix::cellType(GeometryTypes::simplex(1)));
  test.check(::basix::cell::type::interval == Basix::cellType(GeometryTypes::cube(1)));
  test.check(::basix::cell::type::triangle == Basix::cellType(GeometryTypes::triangle));
  test.check(::basix::cell::type::quadrilateral ==Basix::cellType(GeometryTypes::quadrilateral));
  test.check(::basix::cell::type::tetrahedron == Basix::cellType(GeometryTypes::tetrahedron));
  test.check(::basix::cell::type::hexahedron == Basix::cellType(GeometryTypes::hexahedron));
  test.check(::basix::cell::type::prism == Basix::cellType(GeometryTypes::prism));
  test.check(::basix::cell::type::pyramid == Basix::cellType(GeometryTypes::pyramid));

  return test.exit();
}

int main(int argc, char** argv)
{
  using namespace Dune;
  using namespace Dune::Functions;

  bool success = true;

  // BasixLocalFiniteElement<0, RangeClass::scalar> basixlfevertex{
  //   ::basix::element::create_lagrange<double>(::basix::cell::type::point, 0,
  //   ::basix::element::lagrange_variant::equispaced, false)};
  // TEST_FE(basixlfevertex);

  for (int degree = 1; degree < 5; ++degree)
  {
    std::cout << "degree " << degree << std::endl;
    BasixLocalFiniteElement<1, RangeClass::scalar> basixlfeline{
      ::basix::element::create_lagrange<double>(::basix::cell::type::interval, degree,
      ::basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfeline);

    BasixLocalFiniteElement<2, RangeClass::scalar> basixlfetri{
      ::basix::element::create_lagrange<double>(::basix::cell::type::triangle, degree,
      ::basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfetri);

    BasixLocalFiniteElement<3, RangeClass::scalar> basixlfetet{
      ::basix::element::create_lagrange<double>(::basix::cell::type::tetrahedron, degree,
      ::basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfetet);

    BasixLocalFiniteElement<2, RangeClass::scalar> basixlfequad{
      ::basix::element::create_lagrange<double>(::basix::cell::type::quadrilateral, degree,
      ::basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfequad);

    BasixLocalFiniteElement<3, RangeClass::scalar> basixlfehex{
      ::basix::element::create_lagrange<double>(::basix::cell::type::hexahedron, degree,
      ::basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfehex);

    BasixLocalFiniteElement<3, RangeClass::scalar> basixlfeprism{
      ::basix::element::create_lagrange<double>(::basix::cell::type::prism, degree,
      ::basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfeprism);

    BasixLocalFiniteElement<3, RangeClass::scalar> basixlfepyramid{
      ::basix::element::create_lagrange<double>(::basix::cell::type::pyramid, degree,
      ::basix::element::lagrange_variant::equispaced, false)};
    TEST_FE(basixlfepyramid);
  }


  auto refElem = referenceElement<double,2>(GeometryTypes::triangle);
  auto geometry = refElem.geometry<0>(0);
  using Geometry = decltype(geometry);

  BasixFiniteElement<Geometry, RangeClass::scalar> basixfetri{
    ::basix::element::create_lagrange<double>(::basix::cell::type::triangle, 3,
    ::basix::element::lagrange_variant::equispaced, false)};
  basixfetri.bind(geometry);
  TEST_FE(basixfetri);


  std::vector<std::pair<std::string, ::basix::cell::type>> cells{
    {"point", ::basix::cell::type::point},
    {"interval", ::basix::cell::type::interval},
    {"triangle", ::basix::cell::type::triangle},
    {"tetrahedron", ::basix::cell::type::tetrahedron},
    {"quadrilateral", ::basix::cell::type::quadrilateral},
    {"hexahedron", ::basix::cell::type::hexahedron},
    {"prism", ::basix::cell::type::prism},
    {"pyramid", ::basix::cell::type::pyramid}
  };

  std::cout << "Geometry:" << std::endl;
  for (auto [str,type] : cells)
  {
    std::cout << str << ":" << std::endl;

    auto [geo_data,geo_shape] = ::basix::cell::geometry<double>(type);
    auto geometry = ::basix::element::mdspan_t<double,2>{geo_data.data(),geo_shape};
    for (std::size_t i = 0; i < geo_shape[0]; ++i) {
      std::cout << "[ ";
      for (std::size_t j = 0; j < geo_shape[1]; ++j)
        std::cout << geometry(i,j) << " ";
      std::cout << "]" << std::endl;
    }
    std::cout << std::endl;
  }


  std::cout << "Topology:" << std::endl;
  for (auto [str,type] : cells)
  {
    std::cout << str << ":" << std::endl;

    auto topology = ::basix::cell::topology(type);
    for (std::size_t d = 0; d < topology.size(); ++d) {
      std::cout << "  num entities(dim=" << d << ") = " << topology[d].size() << std::endl;
      for (std::size_t s = 0; s < topology[d].size(); ++s) {
        std::cout << "    entity(" << s << ") = [ ";
        for (std::size_t v = 0; v < topology[d][s].size(); ++v) {
          std::cout << topology[d][s][v] << " ";
        }
        std::cout << "]" << std::endl;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  return (success ? 0 : 1) + testUtility();
}
