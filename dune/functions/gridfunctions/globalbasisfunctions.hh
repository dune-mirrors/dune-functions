// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_GRIDFUNCTIONS_GLOBALBASISFUNCTIONS_HH
#define DUNE_FUNCTIONS_GRIDFUNCTIONS_GLOBALBASISFUNCTIONS_HH

#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/functions/gridfunctions/gridviewentityset.hh>

namespace Dune {
namespace Functions {

/**
 * \brief A grid function representing all basis functions
 *
 * \ingroup FunctionImplementations
 *
 * \tparam B Type of global basis
 */
template <typename B>
class GlobalBasisFunctions
{
public:
  using Basis = B;
  using LocalView = typename Basis::LocalView;

  // Only leaf basis supported
  static_assert(LocalView::Tree::isLeaf());

  using GridView = typename Basis::GridView;
  using EntitySet = GridViewEntitySet<GridView, 0>;
  using LocalBasis = typename LocalView::Tree::FiniteElement::Traits::LocalBasisType;

  using Domain = typename EntitySet::GlobalCoordinate;
  using Range = DynamicVector<typename LocalBasis::Traits::RangeFieldType>;

  using LocalDomain = typename EntitySet::LocalCoordinate;
  using Element = typename EntitySet::Element;

  class LocalFunction
  {
  public:
    using Domain = LocalDomain;
    using Range = GlobalBasisFunctions::Range;

    LocalFunction(const std::shared_ptr<const Basis>& basis)
      : basis_(basis)
      , localView_(basis_->localView())
    {}

    /**
     * \brief Bind LocalFunction to grid element.
     */
    void bind(const Element& element)
    {
      localView_.bind(element);
    }

    void unbind()
    {
      localView_.unbind();
    }

    /**
     * \brief Construct a vector of local basis-function evaluations
     */
    Range operator()(const Domain& x) const
    {
      auto y = Range(basis_->dimension(),0);

      auto const& node = localView_.tree();
      auto const& fe = node.finiteElement();

      thread_local std::vector<typename LocalBasis::Traits::RangeType> shapeFunctionValues;
      fe.localBasis().evaluateFunction(x, shapeFunctionValues);
      for (std::size_t i = 0; i < fe.size(); ++i)
        y[localView_.index(node.localIndex(i))] += shapeFunctionValues[i];

      return y;
    }

    const Element& localContext() const
    {
      return localView_.element();
    }

  private:
    std::shared_ptr<const Basis> basis_;
    LocalView localView_;
  };

  GlobalBasisFunctions(Basis const& basis) :
    entitySet_(basis.gridView()),
    basis_(Dune::stackobject_to_shared_ptr(basis))
  {}

  GlobalBasisFunctions(std::shared_ptr<const Basis> basis) :
    entitySet_(basis->gridView()),
    basis_(basis)
  {}

  const Basis& basis() const
  {
    return *basis_;
  }

  Range operator()(const Domain& x) const
  {
    DUNE_THROW(NotImplemented,"not implemented");
  }

  /**
   * \brief Construct local function from a GlobalBasisFunctions
   */
  friend LocalFunction localFunction(const GlobalBasisFunctions& t)
  {
    return LocalFunction(t.basis_);
  }

  /**
   * \brief Get associated EntitySet
   */
  const EntitySet& entitySet() const
  {
    return entitySet_;
  }

private:
  EntitySet entitySet_;
  std::shared_ptr<const Basis> basis_;
};

} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_GRIDFUNCTIONS_GLOBALBASISFUNCTIONS_HH
