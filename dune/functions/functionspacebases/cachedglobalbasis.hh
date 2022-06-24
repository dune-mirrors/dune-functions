// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CACHEDGLOBALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CACHEDGLOBALBASIS_HH

#include <cstddef>
#include <type_traits>
#include <utility>

#include <dune/common/concept.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/concepts.hh>



namespace Dune {
namespace Functions {



/**
 * \brief Global basis with global index caching
 *
 */
template<class B>
class CachedGlobalBasis
{
  using RawBasis = B;

public:

  using PreBasis = typename RawBasis::PreBasis;
  using PrefixPath = typename RawBasis::PrefixPath;
  using GridView = typename RawBasis::GridView;
  using size_type = typename RawBasis::size_type;
  using SizePrefix = typename RawBasis::SizePrefix;
  using MultiIndex = typename RawBasis::MultiIndex;

  class LocalView
  {
  public:

    using GlobalBasis = CachedGlobalBasis<RawBasis>;
    using GridView = typename GlobalBasis::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using size_type = std::size_t;
    using Tree = typename GlobalBasis::PreBasis::Node;
    using MultiIndex = typename GlobalBasis::MultiIndex;

    /** \brief Construct local view for a given global finite element basis */
    LocalView(const GlobalBasis& globalBasis) :
      globalBasis_(&globalBasis),
      tree_(globalBasis_->preBasis().makeNode())
    {
      static_assert(models<Concept::BasisTree<GridView>, Tree>(), "Tree type passed to DefaultLocalView does not model the BasisNode concept.");
      initializeTree(tree_);
    }

    /** \brief Bind the view to a grid element
     *
     * Having to bind the view to an element before being able to actually access any of its data members
     * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
     */
    void bind(const Element& e)
    {
      element_ = e;
      bindTree(tree_, *element_);
      indexCacheOffset_ = globalBasis_->indexCacheOffset(*element_);
    }

    /** \brief Return if the view is bound to a grid element
     */
    bool bound() const
    {
      return static_cast<bool>(element_);
    }

    /** \brief Return the grid element that the view is bound to
     *
     * \throws Dune::Exception if the view is not bound to anything
     */
    const Element& element() const
    {
      return *element_;
    }

    /** \brief Unbind from the current element
     *
     * Calling this method should only be a hint that the view can be unbound.
     */
    void unbind()
    {
      element_.reset();
    }

    /** \brief Return the local ansatz tree associated to the bound entity
     *
     * \returns Tree // This is tree
     */
    const Tree& tree() const
    {
      return tree_;
    }

    /** \brief Total number of degrees of freedom on this element
     */
    size_type size() const
    {
      return tree_.size();
    }

    /**
     * \brief Maximum local size for any element on the GridView
     *
     * This is the maximal size needed for local matrices
     * and local vectors, i.e., the result is
     */
    size_type maxSize() const
    {
      return globalBasis_->preBasis().maxNodeSize();
    }

    //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
    const MultiIndex& index(size_type i) const
    {
      return globalBasis_->indexCache_[indexCacheOffset_ + i];
    }

    /** \brief Return the global basis that we are a view on
     */
    const GlobalBasis& globalBasis() const
    {
      return *globalBasis_;
    }

    const LocalView& rootLocalView() const
    {
      return *this;
    }

  protected:
    friend class DefaultGlobalBasis<PreBasis>;
    const GlobalBasis* globalBasis_;
    std::optional<Element> element_;
    Tree tree_;
    size_type indexCacheOffset_;
  };

  /**
   * \brief Constructor
   */
  CachedGlobalBasis(RawBasis& rawBasis) :
    rawBasis_(rawBasis),
    elementMapper_(rawBasis_.preBasis().gridView(), mcmgElementLayout())
  {
    static_assert(models<Concept::GlobalBasis<GridView>, RawBasis>(), "Type passed to CachedGlobalBasis does not model the GlobalBasis concept.");
    updateIndexCache();
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    return rawBasis_.gridView();
  }

  //! Obtain the pre-basis providing the implementation details
  const PreBasis& preBasis() const
  {
    return rawBasis_.preBasis();
  }

  //! Obtain the pre-basis providing the implementation details
//  PreBasis& preBasis()
//  {
//    return rawBasis_.preBasis();
//  }

  /**
   * \brief Update the stored grid view
   *
   * This will update the indexing information of the global basis.
   * It must be called if the grid has changed.
   */
  void update(const GridView & gv)
  {
    updateIndexCache();
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return rawBasis_.dimension();
  }

  //! Return number of possible values for next position in empty multi index
  size_type size() const
  {
    return rawBasis_.size();
  }

  //! Return number of possible values for next position in multi index
  size_type size(const SizePrefix& prefix) const
  {
    return rawBasis_.size(prefix);
  }

  //! Return local view for basis
  LocalView localView() const
  {
    return LocalView(*this);
  }

  //! Return *this because we are not embedded in a larger basis
  const CachedGlobalBasis& rootBasis() const
  {
    return *this;
  }

  //! Return empty path, because this is the root in the local ansatz tree
  const PrefixPath& prefixPath() const
  {
    return rawBasis_.prefixPath();
  }

protected:

  RawBasis& rawBasis_;


  friend class LocalView;

  using Element = typename GridView::template Codim<0>::Entity;

  size_type indexCacheOffset(const Element& e) const
  {
    auto elementIndex = elementMapper_.index(e);
    return elementIndexOffset_[elementIndex];
  }

  void updateIndexCache()
  {
    elementMapper_.update(preBasis().gridView());
    elementIndexOffset_.resize(elementMapper_.size());
    indexCache_.resize(0);

    auto localView = rawBasis_.localView();
    for (const auto& e : Dune::elements(preBasis().gridView()))
    {
      localView.bind(e);
      auto elementIndex = elementMapper_.index(e);
      elementIndexOffset_[elementIndex] = indexCache_.size();
      for(size_type i=0; i<localView.size(); ++i)
        indexCache_.push_back(localView.index(i));
    }
  }

  Dune::MultipleCodimMultipleGeomTypeMapper<GridView> elementMapper_;
  std::vector<size_type> elementIndexOffset_;
  std::vector<MultiIndex> indexCache_;
};



} // end namespace Functions
} // end namespace Dune



#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_CACHEDGLOBALBASIS_HH
