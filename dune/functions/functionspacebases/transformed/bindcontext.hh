// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_BINDCONTEXT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_BINDCONTEXT_HH

#include <cassert>
#include <optional>

#include <dune/geometry/type.hh>

namespace Dune::Functions {

namespace Impl {

template<class Geometry>
class GeometryContextStorage
{
  public:
    void bindGeometry(Geometry const& geometry)
    {
      geometry_.emplace(geometry);
    }

    Geometry const& geometry() const
    {
      assert(!!geometry_);
      return geometry_.value();
    }

    GeometryType type() const
    {
      return geometry().type();
    }

  private:
    std::optional<Geometry> geometry_;
};

} // namespace Impl

/**
 * \brief Basic bind context backed by a grid element and its geometry.
 *
 * This is the minimal context needed by purely geometry-dependent
 * transformations.  Transformations that need additional data, such as face
 * orientations or global vertex information, can use richer context objects
 * that refine the same concept.
 *
 * The element is borrowed and must outlive the context. Its geometry is copied
 * on bind so transformations use a stable geometry snapshot.
 */
template<class Element>
class ElementBindContext
{
  public:
    using ElementType = Element;
    using Geometry = typename Element::Geometry;

    ElementBindContext() = default;

    explicit ElementBindContext(Element const& element)
    {
      bind(element);
    }

    void bind(Element const& element)
    {
      element_ = &element;
      geometry_.bindGeometry(element.geometry());
    }

    Element const& element() const
    {
      assert(!!element_);
      return *element_;
    }

    Geometry const& geometry() const
    {
      return geometry_.geometry();
    }

    GeometryType type() const
    {
      return geometry().type();
    }

  private:
    Element const* element_ = nullptr;
    Impl::GeometryContextStorage<Geometry> geometry_;
};

/**
 * \brief Basic bind context backed directly by a geometry object.
 *
 * This context is useful for adapters that historically bind a local finite
 * element to a Geometry instead of a full grid element. The geometry is copied
 * on bind and the source object does not need to outlive the context.
 */
template<class Geometry>
class GeometryBindContext
{
  public:
    using GeometryType = Geometry;

    GeometryBindContext() = default;

    explicit GeometryBindContext(Geometry const& geometry)
    {
      bind(geometry);
    }

    void bind(Geometry const& geometry)
    {
      geometry_.bindGeometry(geometry);
    }

    Geometry const& geometry() const
    {
      return geometry_.geometry();
    }

    Dune::GeometryType type() const
    {
      return geometry_.type();
    }

  private:
    Impl::GeometryContextStorage<Geometry> geometry_;
};

} // end namespace Dune::Functions

#endif
