// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_SIMPLEXCONTEXT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_SIMPLEXCONTEXT_HH

#include <array>
#include <cassert>
#include <bitset>
#include <cstddef>
#include <optional>
#include <vector>

#include <dune/common/rangeutilities.hh>
#include <dune/geometry/type.hh>

namespace Dune::Functions {

/**
 * \brief Bind context for simplex finite elements with per-vertex mesh-size data.
 *
 * This context extends the minimal element/geometry context by storing the
 * average mesh size associated to each local vertex.  It is intended for
 * non-affine-equivalent simplex elements such as Hermite, Morley, and Argyris
 * whose basis transformation or interpolation needs local element state beyond
 * the geometry.
 */
template<class Element>
class SimplexElementContext
{
  public:
    using Geometry = typename Element::Geometry;
    using ctype = typename Geometry::ctype;
    static constexpr int dimension = Geometry::mydimension;

    SimplexElementContext() = default;

    template<class VertexMapper>
    SimplexElementContext(Element const& element,
                          VertexMapper const& vertexMapper,
                          std::vector<ctype> const& globalAverageVertexMeshSize)
    {
      bind(element, vertexMapper, globalAverageVertexMeshSize);
    }

    template<class VertexMapper>
    void bind(Element const& element,
              VertexMapper const& vertexMapper,
              std::vector<ctype> const& globalAverageVertexMeshSize)
    {
      bind(element);
      for (auto i : Dune::range(dimension+1))
        averageVertexMeshSize_[i] = globalAverageVertexMeshSize[vertexMapper.subIndex(element, i, dimension)];
    }

    void bind(Element const& element)
    {
      element_ = &element;
      geometry_.emplace(element.geometry());
    }

    template<class ElementMapper>
    void bind(Element const& element,
              ElementMapper const& elementMapper,
              std::vector<std::bitset<dimension+1>> const& edgeOrientations)
    {
      bind(element);
      edgeOrientations_ = edgeOrientations[elementMapper.index(element)];
    }

    Element const& element() const
    {
      assert(!!element_);
      return *element_;
    }

    Geometry const& geometry() const
    {
      assert(!!geometry_);
      return *geometry_;
    }

    Dune::GeometryType type() const
    {
      return element().type();
    }

    ctype averageVertexMeshSize(std::size_t i) const
    {
      return averageVertexMeshSize_[i];
    }

    std::array<ctype,dimension+1> const& averageVertexMeshSizes() const
    {
      return averageVertexMeshSize_;
    }

    bool edgeOrientation(std::size_t i) const
    {
      return edgeOrientations_[i];
    }

    std::bitset<dimension+1> const& edgeOrientations() const
    {
      return edgeOrientations_;
    }

  private:
    Element const* element_ = nullptr;
    std::optional<Geometry> geometry_;
    std::array<ctype,dimension+1> averageVertexMeshSize_{};
    std::bitset<dimension+1> edgeOrientations_;
};

} // end namespace Dune::Functions

#endif
