// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_FINEGRIDVIEWFUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_FINEGRIDVIEWFUNCTION_HH

#include <optional>
#include <type_traits>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/referencehelper.hh>

#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>

namespace Dune::Functions {



/**
 * \brief A wrapper representing a fine grid function on a gridview
 *
 * \ingroup FunctionImplementations
 *
 * \tparam GridFunction Type of the wrapped grid function
 * \tparam GV Type of the target grid view this function should act on
 *
 * This wraps a grid function such that it can be used as a `GridViewFunction`
 * on a user-provided `GridView` under the following assumptions:
 * 1. The grid function's entity set and the `GridView` belong to the same grid.
 * 2. The entity set is finer than the `GridView` in the sense that any
 *    element from the entity has an ancestor in the `GridView`.
 */
template<class GridFunction, class GV, template<class> class DerivativeTraits=Dune::Functions::DefaultDerivativeTraits>
class FineGridViewFunction
{
  using RawGridFunction = Dune::ResolveRef_t<GridFunction>;

  auto&& rawFunction() const
  {
    return Dune::resolveRef(function_);
  }

public:

  using GridView = GV;
  using EntitySet = Dune::Functions::GridViewEntitySet<GridView, 0>;
  using Element = typename EntitySet::Element;
  using Domain = typename EntitySet::GlobalCoordinate;
  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Range = std::decay_t<decltype(std::declval<RawGridFunction>()(std::declval<Domain>()))>;

private:

  using FineEntitySet = std::decay_t<decltype(std::declval<RawGridFunction>().entitySet())>;
  using Traits = Dune::Functions::Imp::GridFunctionTraits<Range(Domain), EntitySet, DerivativeTraits, 56>;

  class FineGridViewLocalFunction
  {
    using Traits = typename FineGridViewFunction::Traits::LocalFunctionTraits;

  public:

    using Derivative = decltype(localFunction(derivative(std::declval<FineGridViewFunction>())));

    /**
     * \brief Construct the LocalFunction
     *
     * The LocalFunction is created from the global FineGridViewFunction.
     **/
    FineGridViewLocalFunction(typename RawGridFunction::LocalFunction&& localFunction, const FineEntitySet& fineEntitySet)
      : localFunction_(localFunction)
      , fineEntitySet_(fineEntitySet)
      , element_()
    {}

    //! Bind to an element from the GridView
    void bind(const Element& element)
    {
      element_ = element;
      forwardToFineFunction_ = fineEntitySet_.contains(*element_);
      if (forwardToFineFunction_)
        localFunction_.bind(element);
    }

    //! \brief Unbind the inner local-functions.
    void unbind()
    {
      element_.reset();
    }

    //! Return if the local function is bound to an element of the GridView
    bool bound() const
    {
      return static_cast<bool>(element_);
    }

    //! Obtain the grid element this function is bound to
    const Element& localContext() const
    {
      return *element_;
    }

    //! Obtain local derivative of this function
    friend auto derivative(const FineGridViewLocalFunction& f)
    {
      if constexpr(requires{ derivative(f.localFunction_); })
        return Derivative(derivative(f.localFunction_), f.fineEntitySet_);
      else
        return typename Traits::DerivativeInterface{};
    }

    //! Evaluate function in local coordinates
    Range operator()(LocalDomain x) const
    {
      if (forwardToFineFunction_)
        return localFunction_(x);
      return evaluateInDescendent(*element_, x);
    }

  private:

    // Find a child containing the point and evaluate there recursively
    Range evaluateInDescendent(const Element& element, LocalDomain x) const
    {
      for(const auto& childElement : descendantElements(element, element.level()+1))
      {
        auto&& geometry = childElement.geometryInFather();
        auto x_child = geometry.local(x);
        if (referenceElement(geometry).checkInside(x_child))
        {
          if (fineEntitySet_.contains(childElement))
          {
            localFunction_.bind(childElement);
            return localFunction_(x_child);
          }
          else
            return evaluateInDescendent(childElement, x_child);
        }
      }
      DUNE_THROW(Dune::Exception, "Did not find matching child for point " << x);
    }

    mutable typename RawGridFunction::LocalFunction localFunction_;
    const FineEntitySet& fineEntitySet_;
    bool forwardToFineFunction_ = false;
    std::optional<Element> element_;
  };

public:

  using LocalFunction = FineGridViewLocalFunction;

  /**
   * \brief Create FineGridViewFunction from GridFunction and GridView
   *
   * \param gridFunction The GridFunction that should be represented on gridView
   * \param gridView The GridFunction should be represented on this gridView
   */
  FineGridViewFunction(const GridFunction& function, const GridView& gridView)
    : function_(function)
    , entitySet_(gridView)
  {}

  /**
   * \brief Create FineGridViewFunction from GridFunction and GridView
   *
   * \param gridFunction The GridFunction that should be represented on gridView
   * \param gridView The GridFunction should be represented on this gridView
   */
  FineGridViewFunction(GridFunction&& function, const GridView& gridView)
    : function_(std::move(function))
    , entitySet_(gridView)
  {}

  //! Evaluate function in global coordinates
  Range operator()(const Domain& x) const
  {
    return function_(x);
  }

  //! Obtain global derivative of this function
  friend auto derivative(const FineGridViewFunction& f)
  {
    if constexpr(requires{ derivative(f.rawFunction()); })
    {
      using RawDerivative = std::decay_t<decltype(derivative(f.rawFunction()))>;
      return FineGridViewFunction<RawDerivative, GridView, DerivativeTraits>(derivative(f.rawFunction()), f.entitySet_.gridView());
    }
    else
      return typename Traits::DerivativeInterface{};
  }

  //! Create a LocalFunction for evaluation in local coordinates
  friend LocalFunction localFunction(const FineGridViewFunction& f)
  {
    return LocalFunction(localFunction(f.rawFunction()), f.rawFunction().entitySet());
  }

  //! Return the EntitySet associated to this GridViewFunction
  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

protected:

  GridFunction function_;
  EntitySet entitySet_;
};



} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_FINEGRIDVIEWFUNCTION_HH
