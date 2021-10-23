// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_COARSENEDGRIDFUNCTION_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_COARSENEDGRIDFUNCTION_HH

#include <functional>
#include <type_traits>

#include <dune/functions/common/referencehelper.hh>
#include <dune/functions/gridfunctions/gridviewentityset.hh>

namespace Dune {
namespace Functions {


/**
 * \brief A grid function defining data on a coarser entity level than bound to
 *
 * \ingroup FunctionImplementations
 *
 * This wrapper function allows to bind the wrapped grid function to elements on a finer
 * grid than this grid function was defined on. It works by traversing from the refined element
 * through the hierarchy until a coarser element is found that is in the entity-set of the
 * bound grid function.
 *
 * \tparam GV  The GridView the wrapper is bound to
 * \tparam GF  A GridFunction that is bound to a coarser GridView
 */
template <class GV, class GF>
class CoarsenedGridFunction
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
  class CoarsenedLocalFunction
  {
  public:
    using Domain = LocalDomain;
    using Range = CoarsenedGridFunction::Range;
    using Element = CoarsenedGridFunction::Element;

    using Mapping = std::function<Domain(const Domain&)>;
    struct ChildElement {
      typename Element::HierarchicIterator it;
      Mapping local;
    };

  public:
    /** \brief Constructor. Stores the `localFct` by value. */
    CoarsenedLocalFunction (const ES& entitySet, const LF& localFct, int maxLevel)
      : entitySet_(entitySet)
      , localFct_(localFct)
      , maxLevel_(maxLevel)
    {}

    CoarsenedLocalFunction (const ES& entitySet, LF&& localFct, int maxLevel)
      : entitySet_(entitySet)
      , localFct_(std::move(localFct))
      , maxLevel_(maxLevel)
    {}

    /**
     * \brief Bind the wrapped local-function to grid element.
     *
     * The given `element` might not be the element the wrapped local-function can be
     * bound to. It has to be iterated in the hierarchy down (to finer elements) until the
     * entity-set of the local function contains the found father element.
     */
    void bind (const Element& element)
    {
      element_ = &element;

      for (auto it = element.hbegin(maxLevel_); it != element.hend(maxLevel_); ++it) {
        if (entitySet_.contains(*it)) {
          childs_.emplace_back(ChildElement{.it=it, .local=makeMapping(element,*it)});
        }
      }
    }

    /** \brief Unbind the wrapped local-function */
    void unbind ()
    {
      childs_.clear();
    }

    /** \brief Evaluate LocalFunction at bound element. */
    Range operator() (const Domain& x) const
    {
      for (auto const& child : childs_) {
        auto refElem = referenceElement(*child.it);
        auto local = child.local(x);
        if (refElem.checkInside(local)) {
          localFct_.bind(*child.it);
          return localFct_(local);
        }
      }

      DUNE_THROW(Dune::Exception, "No child element found!");
      return Range(0);
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
      -> CoarsenedLocalFunction<ES,
          std::decay_t<decltype(derivative(std::declval<const LF&>()))>>
    {
      return {entitySet_, derivative(localFct_)};
    }
#endif

  protected:
    // Construct the coordinate mapping
    Mapping makeMapping (const Element& coarse, Element fine) const
    {
      Mapping map = [](const Domain& x) { return x; };
      while (coarse != fine && fine.hasFather()) {
        map = [map, geo=fine.geometryInFather()](const Domain& x) {
          return map(geo.local(x));
        };
        fine = fine.father();
      }
      return map;
    }

  private:
    const ES& entitySet_;
    mutable LF localFct_;
    int maxLevel_;
    const Element* element_ = nullptr;

    std::vector<ChildElement> childs_;
  };


  /**
   * \brief Constructor.
   *
   * Stores the GridFunction by value. Use a reference wrapper if copying the GridFunction
   * is an expensive operation!
   **/
  CoarsenedGridFunction (const GV& gridView, const GF& gridFct)
    : gridView_{gridView}
    , entitySet_{gridView}
    , gridFct_{gridFct}
  {}

  CoarsenedGridFunction (const GV& gridView, GF&& gridFct)
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
    -> CoarsenedGridFunction<GV,
        std::decay_t<decltype(derivative(resolveRef(std::declval<const GF&>())))>>
  {
    return {gridView_, derivative(resolveRef(gridFct_))};
  }
#endif

  /**
   * \brief Construct local function from a DiscreteGlobalBasisFunction
   * \ingroup FunctionImplementations
   */
  friend auto localFunction (const CoarsenedGridFunction& t)
  {
    using ES = std::decay_t<decltype(t.gridFct_.entitySet())>;
    using LF = std::decay_t<decltype(localFunction(resolveRef(t.gridFct_)))>;
    using LocalFunction = CoarsenedLocalFunction<ES,LF>;
    return LocalFunction{t.gridFct_.entitySet(), localFunction(resolveRef(t.gridFct_)),
                         t.gridView_.grid().maxLevel()};
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

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_COARSENEDGRIDFUNCTION_HH
