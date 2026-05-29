// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_QUADRATUREBASISFUNCTIONCACHE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_QUADRATUREBASISFUNCTIONCACHE_HH

#include <map>
#include <tuple>
#include <vector>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/common/typetree/traversal.hh>
#include <dune/common/typetree/treecontainer.hh>
#include <dune/common/typetree/treepath.hh>

namespace Dune::Functions
{
  /**
   * \brief Cache the evaluations of basis functions and their Jacobians in all quadrature points.
   *
   * The cache is constructed over a basis-tree and allows to store vectors of evaluations for all
   * basis functions on the leaf nodes of that tree. The whole tree cache is initialized at
   * once.
   *
   * \b Example
   * \code{.cpp}
    auto localView = basis.localView();
    auto cache = QuadratureBasisFunctionCache<LocalView::Tree, Derivatives::Value>{};
    for (auto& e : elements(basis.gridView()))
    {
      localView.bind(e);
      auto& quad = QuadratureRules<double,dim>::rule(e.type(), 4);

      cache.bind(localView, quad);
      for (std::size_t iq = 0; iq < quad.size(); ++iq)
        auto const& values = cache.evaluated(TypeTree::makeTreePath(0))[iq].get(Derivatives::Value{});
    }
   * \endcode
   */
  template <class Tree, class... CachedDerivatives>
  class QuadratureBasisFunctionCache
  {
    template <class LeafNode>
    struct LeafNodePointPrecomputeCache
    {
      using Basis = typename LeafNode::GlobalizedFiniteElement::Basis;

      //! Buffer type stored in the cache for the quantity selected by D.
      template <class D>
      using PrecomputeBuffer = typename Basis::template PrecomputeBuffer<D>;

      //! Public output range type for the quantity selected by D.
      template <class D>
      using DerivativeRange = typename Basis::template DerivativeRange<D>;

      template <class Derivative, class Domain>
      void evaluate (Derivative const& d, Domain const& x,
                     std::vector<DerivativeRange<Derivative>>& out) const
      {
        finalize(d,x,get<Derivative>(),out);
      }

      template <class Derivative>
      PrecomputeBuffer<Derivative>& get (Derivative const& d = {})
      {
        // NOTE: This does not respect the value of d, just its type.
        return std::get<PrecomputeBuffer<Derivative>>(cache_);
      }

      template <class Derivative>
      PrecomputeBuffer<Derivative> const& get (Derivative const& d = {}) const
      {
        return std::get<PrecomputeBuffer<Derivative>>(cache_);
      }

      std::tuple<PrecomputeBuffer<CachedDerivatives>...> cache_;
    };


    template <class LeafNode>
    struct LeafNodePointEvaluationCache
    {
      using Basis = typename LeafNode::GlobalizedFiniteElement::Basis;

      //! Public output range type for the quantity selected by D.
      template <class D>
      using DerivativeRange = typename Basis::template DerivativeRange<D>;

      //! Container type stored in the cache for the quantity selected by D.
      template <class D>
      using EvaluationContainer = std::vector<DerivativeRange<D>>;

      template <class Derivative>
      EvaluationContainer<Derivative>& get (Derivative const& d = {})
      {
        return std::get<EvaluationContainer<Derivative>>(cache_);
      }

      template <class Derivative>
      EvaluationContainer<Derivative> const& get (Derivative const& d = {}) const
      {
        return std::get<EvaluationContainer<Derivative>>(cache_);
      }

      std::tuple<EvaluationContainer<CachedDerivatives>...> cache_;
    };


    template <class LeafNode>
    using LeafNodeQuadraturePrecomputeCache = std::vector<LeafNodePointPrecomputeCache<LeafNode>>;

    template <class LeafNode>
    using LeafNodeQuadratureEvaluationCache = std::vector<LeafNodePointEvaluationCache<LeafNode>>;

    struct LeafNodeQuadraturePrecomputeCacheFactory
    {
      template <class LeafNode>
      auto operator() (LeafNode const& /*leafNode*/) const
      {
        return LeafNodeQuadraturePrecomputeCache<LeafNode>{};
      }
    };

    struct LeafNodeQuadratureEvaluationCacheFactory
    {
      template <class LeafNode>
      auto operator() (LeafNode const& /*leafNode*/) const
      {
        return LeafNodeQuadratureEvaluationCache<LeafNode>{};
      }
    };

    using Key = std::tuple<GeometryType, int, int>;

    template <class Element, class ct, int dim>
    Key getKey (Element const& e, QuadratureRule<ct, dim> const& quadRule)
    {
      GeometryType gt = e.type();
      int order = quadRule.order();
      return Key{gt,order,dim};
    }


  public:
    QuadratureBasisFunctionCache () = default;

    QuadratureBasisFunctionCache (CachedDerivatives const&... d)
      : cachedDerivatives_(d...)
    {}

    QuadratureBasisFunctionCache (QuadratureBasisFunctionCache const& other)
      : cachedDerivatives_(other.cachedDerivatives_)
      , precomputeCache_(other.precomputeCache_)
      , evaluationCache_(other.evaluationCache_)
    {}

    QuadratureBasisFunctionCache& operator= (QuadratureBasisFunctionCache const& other)
    {
      cachedDerivatives_ = other.cachedDerivatives_;
      precomputeCache_ = other.precomputeCache_;
      evaluationCache_ = other.evaluationCache_;
      elementPrecomputeCacheIterator_ = precomputeCache_.end();
      elementEvaluationCacheIterator_ = evaluationCache_.end();
    }

    QuadratureBasisFunctionCache (QuadratureBasisFunctionCache&& other)
      : cachedDerivatives_(std::move(other.cachedDerivatives_))
      , precomputeCache_(std::move(other.precomputeCache_))
      , evaluationCache_(std::move(other.evaluationCache_))
    {}

    QuadratureBasisFunctionCache& operator= (QuadratureBasisFunctionCache&& other)
    {
      cachedDerivatives_ = std::move(other.cachedDerivatives_);
      precomputeCache_ = std::move(other.precomputeCache_);
      evaluationCache_ = std::move(other.evaluationCache_);
      elementPrecomputeCacheIterator_ = precomputeCache_.end();
      elementEvaluationCacheIterator_ = evaluationCache_.end();
    }

    template <class TreePath>
    auto const& precomputed (TreePath const& tp) const
    {
      assert(elementPrecomputeCacheIterator_ != precomputeCache_.end());
      return (elementPrecomputeCacheIterator_->second)[tp];
    }

    template <class TreePath>
    auto const& evaluated (TreePath const& tp) const
    {
      assert(elementEvaluationCacheIterator_ != evaluationCache_.end());
      return (elementEvaluationCacheIterator_->second)[tp];
    }

    template <class LocalView, class ct, int dim>
    void bind (LocalView const& localView, QuadratureRule<ct, dim> const& quadRule)
    {
      auto key = getKey(localView.element(), quadRule);
      auto pit = precomputeCache_.find(key);
      if (pit == precomputeCache_.end())
      {
        bool inserted = false;
        std::tie(pit, inserted) = precomputeCache_.emplace(key,
          TypeTree::makeTreeContainer(localView.tree(), LeafNodeQuadraturePrecomputeCacheFactory{}));
        assert(inserted);

        // precompute the values
        TypeTree::forEachLeafNode(localView.tree(), [&](auto&& node, auto&& tp)
        {
          auto& data = (pit->second)[tp];
          data.resize(quadRule.size());
          for (std::size_t iq = 0; iq < quadRule.size(); ++iq) {
            std::apply([&](auto... d) {
              (basis(node).precompute(d, quadRule[iq].position(), data[iq].get(d)),...);
            }, cachedDerivatives_);
          }
        });
      }
      elementPrecomputeCacheIterator_ = pit;

      // prepare the evaluation container
      auto eit = evaluationCache_.find(key);
      if (eit == evaluationCache_.end())
      {
        bool inserted = false;
        std::tie(eit, inserted) = evaluationCache_.emplace(key,
          TypeTree::makeTreeContainer(localView.tree(), LeafNodeQuadratureEvaluationCacheFactory{}));
        assert(inserted);
      }
      elementEvaluationCacheIterator_ = eit;

      // Finalize the evaluation on the current element
      TypeTree::forEachLeafNode(localView.tree(), [&](auto&& node, auto&& tp)
      {
        auto& cache = (pit->second)[tp];
        auto& data = (eit->second)[tp];
        data.resize(quadRule.size());
        for (std::size_t iq = 0; iq < quadRule.size(); ++iq) {
          std::apply([&](auto... d) {
            (basis(node).finalize(d, quadRule[iq].position(), cache[iq].get(d), data[iq].get(d)),...);
          }, cachedDerivatives_);
        }
      });
    }

  private:
    template <class Node>
    auto const& basis (const Node& node) const
    {
      return node.globalizedFiniteElement().basis();
    }

  private:
    std::tuple<CachedDerivatives...> cachedDerivatives_;

    using PrecomputeCache = TypeTree::TreeContainer<LeafNodeQuadraturePrecomputeCache, Tree>;
    std::map<Key,PrecomputeCache> precomputeCache_ = {};

    using EvaluationCache = TypeTree::TreeContainer<LeafNodeQuadratureEvaluationCache, Tree>;
    std::map<Key,EvaluationCache> evaluationCache_ = {};

    typename std::map<Key,PrecomputeCache>::const_iterator elementPrecomputeCacheIterator_ = {};
    typename std::map<Key,EvaluationCache>::const_iterator elementEvaluationCacheIterator_ = {};
  };

} // end namespace Dune::Functions

#endif