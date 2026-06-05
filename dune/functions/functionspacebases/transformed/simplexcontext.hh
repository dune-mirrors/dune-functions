// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_SIMPLEXCONTEXT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_SIMPLEXCONTEXT_HH

#include <array>
#include <bitset>
#include <cstddef>
#include <vector>

#include <dune/common/rangeutilities.hh>

#include <dune/functions/functionspacebases/transformed/bindcontext.hh>

namespace Dune::Functions {

/**
 * \brief Bind context for simplex elements requiring per-vertex mesh sizes.
 *
 * Every bind operation initializes the element, geometry, and mesh-size data.
 */
template<class Element>
class SimplexVertexMeshSizeContext
{
    using ElementContext = ElementBindContext<Element>;

  public:
    using Geometry = typename ElementContext::Geometry;
    using ctype = typename Geometry::ctype;
    static constexpr int dimension = Geometry::mydimension;

    SimplexVertexMeshSizeContext() = default;

    template<class VertexMapper>
    SimplexVertexMeshSizeContext(
        Element const& element,
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
      elementContext_.bind(element);
      for (auto i : Dune::range(dimension+1))
        averageVertexMeshSize_[i] =
          globalAverageVertexMeshSize[vertexMapper.subIndex(element, i, dimension)];
    }

    Element const& element() const
    {
      return elementContext_.element();
    }

    Geometry const& geometry() const
    {
      return elementContext_.geometry();
    }

    Dune::GeometryType type() const
    {
      return elementContext_.type();
    }

    ctype averageVertexMeshSize(std::size_t i) const
    {
      return averageVertexMeshSize_[i];
    }

    std::array<ctype,dimension+1> const& averageVertexMeshSizes() const
    {
      return averageVertexMeshSize_;
    }

  private:
    ElementContext elementContext_;
    std::array<ctype,dimension+1> averageVertexMeshSize_{};
};

/**
 * \brief Bind context for simplex elements requiring edge orientations.
 *
 * Every bind operation initializes the element, geometry, and orientation data.
 */
template<class Element>
class SimplexEdgeOrientationContext
{
    using ElementContext = ElementBindContext<Element>;

  public:
    using Geometry = typename ElementContext::Geometry;
    static constexpr int dimension = Geometry::mydimension;
    static constexpr int edgeCount = dimension*(dimension+1)/2;

    SimplexEdgeOrientationContext() = default;

    template<class ElementMapper>
    SimplexEdgeOrientationContext(
        Element const& element,
        ElementMapper const& elementMapper,
        std::vector<std::bitset<edgeCount>> const& edgeOrientations)
    {
      bind(element, elementMapper, edgeOrientations);
    }

    template<class ElementMapper>
    void bind(Element const& element,
              ElementMapper const& elementMapper,
              std::vector<std::bitset<edgeCount>> const& edgeOrientations)
    {
      elementContext_.bind(element);
      edgeOrientations_ = edgeOrientations[elementMapper.index(element)];
    }

    Element const& element() const
    {
      return elementContext_.element();
    }

    Geometry const& geometry() const
    {
      return elementContext_.geometry();
    }

    Dune::GeometryType type() const
    {
      return elementContext_.type();
    }

    bool edgeOrientation(std::size_t i) const
    {
      return edgeOrientations_[i];
    }

    std::bitset<edgeCount> const& edgeOrientations() const
    {
      return edgeOrientations_;
    }

  private:
    ElementContext elementContext_;
    std::bitset<edgeCount> edgeOrientations_;
};

/**
 * \brief Bind context for simplex elements requiring vertex sizes and orientations.
 *
 * Every bind operation initializes all context data required by elements such
 * as Argyris.
 */
template<class Element>
class SimplexVertexMeshSizeAndEdgeOrientationContext
{
    using ElementContext = ElementBindContext<Element>;

  public:
    using Geometry = typename ElementContext::Geometry;
    using ctype = typename Geometry::ctype;
    static constexpr int dimension = Geometry::mydimension;
    static constexpr int edgeCount = dimension*(dimension+1)/2;

    SimplexVertexMeshSizeAndEdgeOrientationContext() = default;

    template<class VertexMapper, class ElementMapper>
    void bind(Element const& element,
              VertexMapper const& vertexMapper,
              std::vector<ctype> const& globalAverageVertexMeshSize,
              ElementMapper const& elementMapper,
              std::vector<std::bitset<edgeCount>> const& edgeOrientations)
    {
      elementContext_.bind(element);
      for (auto i : Dune::range(dimension+1))
        averageVertexMeshSize_[i] =
          globalAverageVertexMeshSize[vertexMapper.subIndex(element, i, dimension)];
      edgeOrientations_ = edgeOrientations[elementMapper.index(element)];
    }

    Element const& element() const
    {
      return elementContext_.element();
    }

    Geometry const& geometry() const
    {
      return elementContext_.geometry();
    }

    Dune::GeometryType type() const
    {
      return elementContext_.type();
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

    std::bitset<edgeCount> const& edgeOrientations() const
    {
      return edgeOrientations_;
    }

  private:
    ElementContext elementContext_;
    std::array<ctype,dimension+1> averageVertexMeshSize_{};
    std::bitset<edgeCount> edgeOrientations_;
};

} // end namespace Dune::Functions

#endif
