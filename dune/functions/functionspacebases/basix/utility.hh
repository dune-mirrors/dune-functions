// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_BASIX_UTILITY_HH
#define DUNE_FUNCTIONS_BASIX_UTILITY_HH

#if HAVE_BASIX

#include <dune/geometry/type.hh>

#include <basix/cell.h>
#include <basix/indexing.h>

namespace Dune::Functions::Basix {

  static inline constexpr int number_of_cell_types = 8;

  // Map the cell::type from basix to the Dune GeometryType
  inline GeometryType geometryType (::basix::cell::type cell_type)
  {
    switch (cell_type) {
      case ::basix::cell::type::point:          return GeometryTypes::vertex;
      case ::basix::cell::type::interval:       return GeometryTypes::line;
      case ::basix::cell::type::triangle:       return GeometryTypes::triangle;
      case ::basix::cell::type::tetrahedron:    return GeometryTypes::tetrahedron;
      case ::basix::cell::type::quadrilateral:  return GeometryTypes::quadrilateral;
      case ::basix::cell::type::hexahedron:     return GeometryTypes::hexahedron;
      case ::basix::cell::type::prism:          return GeometryTypes::prism;
      case ::basix::cell::type::pyramid:        return GeometryTypes::pyramid;
      default: return GeometryTypes::none(0);
    }
  }

  // Map the Dune::GeometryType to cell::type from basix
  inline ::basix::cell::type cellType (GeometryType type)
  {
    if (type.isVertex())
      return ::basix::cell::type::point;
    else if (type.isLine())
      return ::basix::cell::type::interval;
    else if (type.isTriangle())
      return ::basix::cell::type::triangle;
    else if (type.isTetrahedron())
      return ::basix::cell::type::tetrahedron;
    else if (type.isQuadrilateral())
      return ::basix::cell::type::quadrilateral;
    else if (type.isHexahedron())
      return ::basix::cell::type::hexahedron;
    else if (type.isPrism())
      return ::basix::cell::type::prism;
    else if (type.isPyramid())
      return ::basix::cell::type::pyramid;
    else {
      DUNE_THROW(Dune::NotImplemented, "GeometryType not representable as ::basix::cell::type");
      return ::basix::cell::type{};
    }
  }

  // Map the entity index of basix into the Dune numbering
  inline int entityIndex (::basix::cell::type cell_type, int dim, int s)
  {
    // {cell_type, dimension, entity}
    static constexpr int perm[number_of_cell_types][4][12]{
      { {0} }, // point
      { {0,1}, {0} }, // interval
      { {0,1,2}, {2,1,0}, {0} }, // triangle
      { {0,1,2,3}, {5,4,2,3,1,0}, {3,2,1,0}, {0} }, // tetrahedron
      { {0,1,2,3}, {2,0,1,3}, {0} }, // quadrilateral
      { {0,1,2,3,4,5,6,7}, {6,4,0,5,1,7,2,3,10,8,9,11}, {4,2,0,1,3,5}, {0} }, // hexahedron
      { {0,1,2,3,4,5}, {3,4,0,5,1,2,6,7,8}, {3,0,1,2,4}, {0} }, // prism
      { {0,1,2,3,4}, {2,0,4,1,5,3,6,7}, {0,3,1,2,4}, {0} } // pyramid
    };

    return perm[(int)(cell_type)][dim][s];
  }

  inline int entityIndex (GeometryType type, int dim, int s)
  {
    // {cell_type, dimension, entity}
    static constexpr int perm[number_of_cell_types][4][12]{
      { {0} }, // point
      { {0,1}, {0} }, // interval
      { {0,1,2}, {2,1,0}, {0} }, // triangle
      { {0,1,2,3}, {5,4,2,3,1,0}, {3,2,1,0}, {0} }, // tetrahedron
      { {0,1,2,3}, {1,2,0,3}, {0} }, // quadrilateral
      { {0,1,2,3,4,5,6,7}, {2,4,6,7,1,3,0,5,9,10,8,11}, {2,3,1,4,0,5}, {0} }, // hexahedron
      { {0,1,2,3,4,5}, {2,4,5,0,1,3,6,7,8}, {1,2,3,0,4}, {0} }, // prism
      { {0,1,2,3,4}, {1,3,0,5,2,4,6,7}, {0,2,3,1,4}, {0} } // pyramid
    };

    ::basix::cell::type cell_type = cellType(type);
    return perm[(int)(cell_type)][dim][s];
  }


  // Map the derivative-order tuple from the partial() method to the
  // derivative index used in basix
  template <std::size_t dim>
  int indexing (const std::array<unsigned int,dim>& orders)
  {
    if constexpr (dim == 1)
      return ::basix::indexing::idx(orders[0]);
    else if constexpr (dim == 2)
      return ::basix::indexing::idx(orders[0],orders[1]);
    else if constexpr (dim == 3)
      return ::basix::indexing::idx(orders[0],orders[1],orders[2]);
    else
      return 0;
  }

  // Map the first order derivative in direction d to the derivative index
  // used in basix
  template <std::size_t dim>
  int indexing (unsigned int d)
  {
    std::array<unsigned int,dim> orders{};
    orders[d] = 1;
    return indexing(orders);
  }

  template <class T>
  std::array<std::size_t,1> extents (const std::vector<T>& range)
  {
    return {range.size()};
  }

  template <class T, int n>
  std::array<std::size_t,1> extents (const FieldVector<T,n>& range)
  {
    return {std::size_t(n)};
  }

  template <class T, int n, int m>
  std::array<std::size_t,2> extents (const FieldMatrix<T,n,m>& range)
  {
    return {std::size_t(n), std::size_t(m)};
  }

  template <class T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
  std::array<std::size_t,1> extents (const T& range)
  {
    return {1};
  }

  template <std::size_t n>
  inline std::size_t prod (const std::array<std::size_t,n>& shape)
  {
    std::size_t result = 1;
    for (std::size_t s : shape)
      result *= s;
    return result;
  }

  inline std::size_t prod (const std::vector<std::size_t>& shape)
  {
    std::size_t result = 1;
    for (std::size_t s : shape)
      result *= s;
    return result;
  }

} // end namespace Dune::Functions::Basix

#endif // HAVE_BASIX
#endif // DUNE_FUNCTIONS_BASIX_UTILITY_HH