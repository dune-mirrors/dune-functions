// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/common/geometryinancestor.hh>

template<class Element, class GeometryInAncestor>
void checkGeometryInAncestor(Dune::TestSuite& suite, const GeometryInAncestor& geometryInAncestor, const Element& fineElement, const Element& coarseElement)
{
  suite.check(geometryInAncestor.fineElement() == fineElement)
    << "fineElement() does not return the expected element";
  suite.check(geometryInAncestor.coarseElement() == coarseElement)
    << "coarseElement() does not return the expected element";

  auto x = referenceElement(fineElement.geometry()).position(0, 0);
  auto xGlobal = fineElement.geometry().global(x);
  auto xInCoarse = geometryInAncestor.global(x);
  suite.check((coarseElement.geometry().global(xInCoarse) - xGlobal).two_norm() < 1e-12)
    << "global() does not map into coarse element";
  suite.check((geometryInAncestor.local(xInCoarse) - x).two_norm() < 1e-12)
    << "local() is not the inverse of global()";
}

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite suite;

  using Grid = Dune::YaspGrid<2>;
  Grid grid({1.0, 1.0}, {2, 2});
  grid.globalRefine(2);

  using Element = Grid::Codim<0>::Entity;
  auto isNotLevelZero = [](const Element& e) { return e.level() > 0; };

  for (const auto& fineElement : elements(grid.leafGridView()))
  {
    auto coarseElement = fineElement.father().father();

    // Bind to an l-value element
    Dune::Functions::GeometryInAncestor<Element> geometryInAncestor;
    suite.check(geometryInAncestor.bind(fineElement, isNotLevelZero) == coarseElement)
      << "bind() does not return the coarse element";
    checkGeometryInAncestor(suite, geometryInAncestor, fineElement, coarseElement);

    // Copies must refer to their own stored coarse element
    auto copy = geometryInAncestor;
    checkGeometryInAncestor(suite, copy, fineElement, coarseElement);
    suite.check(&copy.coarseElement() != &geometryInAncestor.coarseElement())
      << "copy refers to coarse element stored in original";

    Element otherFineElement = fineElement;
    auto copyWithElement = Dune::Functions::GeometryInAncestor<Element>(geometryInAncestor, otherFineElement);
    checkGeometryInAncestor(suite, copyWithElement, fineElement, coarseElement);
    suite.check(&copyWithElement.fineElement() == &otherFineElement)
      << "copy does not refer to passed fine element";

    // Bind to an r-value element
    Dune::Functions::GeometryInAncestor<Element> geometryInAncestor2;
    suite.check(geometryInAncestor2.bind(Element(fineElement), isNotLevelZero) == coarseElement)
      << "bind() does not return the coarse element";
    checkGeometryInAncestor(suite, geometryInAncestor2, fineElement, coarseElement);
    auto copy2 = geometryInAncestor2;
    checkGeometryInAncestor(suite, copy2, fineElement, coarseElement);
  }

  return suite.exit();
}
