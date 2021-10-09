// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_HIERARCHICGRIDFUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_HIERARCHICGRIDFUNCTION_HH

#include <functional>
#include <type_traits>

#include <dune/functions/common/referencehelper.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>

namespace Dune {
namespace Functions {


/**
 * \brief A grid function defining data on a different entity level than bound to
 *
 * \ingroup FunctionImplementations
 *
 * \tparam GV  The GridView the data is bound to
 * \tparam GF  A GridFunction that is evaluated on the finer GridView
 */
template <class GV, class GF>
class HierarchicGridFunction
{
public:
  /** \brief The set of entities this function can be evaluated on */
  using EntitySet = GridViewEntitySet<GV,0>;

  /** \brief The global coordinates */
  using Domain = typename EntitySet::GlobalCoordinate;

  /** \brief The result type of the function */
  using Range = std::invoke_result_t<GF, Domain>;

  /** \brief The local coordinates */
  using LocalDomain = typename EntitySet::LocalCoordinate;

  /** \brief Type of the grid element the LocalFunction can be bound to */
  using Element = typename EntitySet::Element;


  template <class ES, class LF>
  class HierarchicLocalFunction
  {
  public:
    using Domain = LocalDomain;
    using Range = HierarchicGridFunction::Range;
    using Element = HierarchicGridFunction::Element;

  public:
    /** \brief Constructor. Stores the `localFct` by value. */
    HierarchicLocalFunction (const ES& entitySet, const LF& localFct)
      : entitySet_(entitySet)
      , localFct_(localFct)
    {}

    HierarchicLocalFunction (const ES& entitySet, LF&& localFct)
      : entitySet_(entitySet)
      , localFct_(std::move(localFct))
    {}

    /**
     * \brief Bind the wrapped local-function to grid element.
     *
     * Note: the given `element` might not be the element the wrapped local-function can be
     * bound to.
     */
    void bind (const Element& element)
    {
      element_ = &element;
      father_.emplace(element);
      local_ = [](const Domain& x) { return x; };
      while (!entitySet_.contains(*father_) && father_->hasFather()) {
        // store the coordinate transform
        local_ = [local=local_, geo=father_->geometryInFather()](const Domain& x) {
          return geo.local(local(x));
        };

        father_.emplace(father_->father());
      }

      assert(entitySet_.contains(*father_));
      localFct_.bind(*father_);
    }

    /** \brief Unbind the wrapped local-function */
    void unbind ()
    {
      localFct_.unbind();
    }

    /** \brief Evaluate LocalFunction at bound element. */
    Range operator() (const Domain& x) const
    {
      return localFct_(local_(x));
    }

    /** \brief Return the element this LocalFunction is bound to */
    const Element& localContext () const
    {
      assert(!!element_);
      return *element_;
    }

#if 0
    /** \brief Construct a derivative by wrapping the derivative of the wrapped local-function */
    auto makeDerivative () const
      -> HierarchicLocalFunction<ES,
          std::decay_t<decltype(derivative(std::declval<const LF&>()))>>
    {
      return {entitySet_, derivative(localFct_)};
    }
#endif

  private:
    const ES& entitySet_;
    LF localFct_;
    const Element* element_ = nullptr;
    std::optional<Element> father_ = std::nullopt;
    std::function<Domain(const Domain&)> local_ = nullptr;
  };


  /**
   * \brief Constructor.
   *
   * Stores the GridFunction by value. Use a reference wrapper if copy the GridFunction
   * is an expensive operation!
   **/
  HierarchicGridFunction (const GV& gridView, const GF& gridFct)
    : gridView_{gridView}
    , entitySet_{gridView}
    , gridFct_{gridFct}
  {}

  HierarchicGridFunction (const GV& gridView, GF&& gridFct)
    : gridView_{gridView}
    , entitySet_{gridView}
    , gridFct_{std::move(gridFct)}
  {}

  /** \brief Evaluate the wrapped grid-function. */
  Range operator() (const Domain& x) const
  {
    return gridFct_(x);
  }

#if 0
  /** \brief Construct a derivative by wrapping the derivative of the wrapped grid-function */
  auto makeDerivative () const
    -> HierarchicGridFunction<GV,
        std::decay_t<decltype(derivative(resolveRef(std::declval<const GF&>())))>>
  {
    return {gridView_, derivative(resolveRef(gridFct_))};
  }
#endif

  /**
   * \brief Construct local function from a DiscreteGlobalBasisFunction
   * \ingroup FunctionImplementations
   */
  friend auto localFunction (const HierarchicGridFunction& t)
  {
    using ES = std::decay_t<decltype(t.gridFct_.entitySet())>;
    using LF = std::decay_t<decltype(localFunction(resolveRef(t.gridFct_)))>;
    using LocalFunction = HierarchicLocalFunction<ES,LF>;
    return LocalFunction{t.gridFct_.entitySet(), localFunction(resolveRef(t.gridFct_))};
  }

  /** \brief Get associated EntitySet */
  const EntitySet& entitySet () const
  {
    return entitySet_;
  }

private:
  GV gridView_;
  EntitySet entitySet_;
  GF gridFct_;
};


#if 0
template <class F>
auto derivative (F const& f)
  -> decltype(f.makeDerivative())
{
  return f.makeDerivative();
}
#endif

} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_HIERARCHICGRIDFUNCTION_HH
