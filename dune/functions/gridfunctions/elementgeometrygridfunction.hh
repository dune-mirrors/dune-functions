// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_ELEMENTGEOMETRYGRIDFUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_ELEMENTGEOMETRYGRIDFUNCTION_HH

#include <tuple>
#include <optional>

#include <dune/functions/gridfunctions/facenormalgridfunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/common/defaultderivativetraits.hh>



namespace Dune {
namespace Functions {



/**
 * \brief A grid function mapping local coordinate to element and geometry
 *
 * \ingroup FunctionImplementations
 *
 * For given grid view, this grid functions maps local coordinates
 * of points to a triple of element, element-geometry, and local-coordinate.
 * While being rather strange on its own, it allows to define other
 * grid functions like, e.g., the local element size, by composition
 * with very little effort.
 *
 * \tparam GV The GridView this function is defined on
 */
template<class GV>
class ElementGeometryGridFunction
{
public:
  using GridView = GV;
  using EntitySet = Dune::Functions::GridViewEntitySet<GridView, 0>;
  using Element = typename EntitySet::Element;

  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Domain = typename EntitySet::GlobalCoordinate;
  using Range = std::tuple<const Element&, const typename Element::Geometry&, LocalDomain>;

private:

  using Traits = Dune::Functions::Imp::GridFunctionTraits<Range(Domain), EntitySet, Dune::Functions::DefaultDerivativeTraits, 16>;

  class LocalFunction
  {
    using Geometry = typename Element::Geometry;
    static const int dimension = GV::dimension;
  public:

    //! Bind the local function to a grid element
    void bind(const Element& element)
    {
      element_ = element;
      geometry_.emplace(element_.geometry());
    }

    //! Unbind the local function
    void unbind()
    {
      geometry_.reset();
    }

    //! Return if the local function is bound to a grid element
    bool bound() const
    {
      return static_cast<bool>(geometry_);
    }

    /**
     * \brief Evaluate local function.
     *
     * This return a triple of element, element-geometry, and local coordinate.
     * While element and element-geometry are returned by reference, the
     * local coordinate is returned by value.
     *
     * \returns std::tuple of element, element-geometry, and local coordinate
     */
    Range operator()(const LocalDomain& x) const
    {
      return {element_, *geometry_, x};
    }

    //! Return the bound element stored as copy in the \ref bind function.
    const Element& localContext() const
    {
      return element_;
    }

    //! There's no derivate of ElementGeometryGridFunction::LocalFunctionTraits
    friend typename Traits::LocalFunctionTraits::DerivativeInterface derivative(const LocalFunction& t)
    {
      DUNE_THROW(NotImplemented, "Not implemented");
    }

  private:
    std::optional<Geometry> geometry_;
    Element element_;
  };

public:

  //! Create an ElementGeometryGridFunction
  ElementGeometryGridFunction(const GridView& gridView) :
    entitySet_(gridView)
  {}

  //! Evaluation of ElementGeometryGridFunction in global coordinates is not implemented
  Range operator()(const Domain& x) const
  {
    DUNE_THROW(NotImplemented, "Not implemented");
  }

  //! There's no derivate of ElementGeometryGridFunction
  friend typename Traits::DerivativeInterface derivative(const ElementGeometryGridFunction& t)
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  //! Create an ElementGeometryGridFunction::LocalFunction
  friend LocalFunction localFunction(const ElementGeometryGridFunction& t)
  {
    return LocalFunction{};
  }

  //! Return the entity set this function is defined on
  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

private:
  EntitySet entitySet_;
};



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_ELEMENTGEOMETRYGRIDFUNCTION_HH
