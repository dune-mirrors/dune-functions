// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_QUADRATUREBASISFUNCTIONCACHE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_QUADRATUREBASISFUNCTIONCACHE_HH

#include <concepts>
#include <map>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>
#include <dune/common/overloadset.hh>
#include <dune/common/typetree/traversal.hh>
#include <dune/common/typetree/treecontainer.hh>
#include <dune/common/typetree/treepath.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

namespace Dune::Functions
{
  /**
   * \brief Cache selected physical-basis evaluations at quadrature points.
   *
   * Each leaf cache exposes the bound physical basis and caches its reference
   * precomputations independently for every derivative tag. Quadrature rules are
   * identified approximately by geometry type, delivered order, and point count.
   * Finalization is repeated after every bind and therefore always uses the
   * currently bound physical element.
   *
   * \b Example
   * \code{.cpp}
    auto localView = basis.localView();
    auto cache = QuadratureBasisFunctionCache<LocalView::Tree, Derivatives::Value>{};
    for (auto& e : elements(basis.gridView()))
    {
      localView.bind(e);
      cache.bind(localView);

      auto& quad = QuadratureRules<double,dim>::rule(e.type(), 4);
      auto& values = cache.evaluate(Derivatives::Value{}, quad);
      for (std::size_t iq = 0; iq < quad.size(); ++iq)
        auto const& point_values = values[iq];
    }
   * \endcode
   */
  template <class Tree, class... Derivatives>
  class QuadratureBasisFunctionCache
  {
    template <class Node>
    class QuadratureNodeCache
    {
      using Basis = typename Node::FiniteElement::PhysicalBasis;

      using QuadratureKey = std::tuple<
        GeometryType,
        int,           // delivered order
        std::size_t>;  // number of points

      //! Buffer type stored in the cache for the quantity selected by D.
      template <class D>
      using PrecomputeBuffer = std::vector<typename Basis::template PrecomputeBuffer<D>>;

      template<class D>
      struct PrecomputeCache
      {
        std::map<QuadratureKey,PrecomputeBuffer<D>> buffers;
      };

      using PrecomputeCaches = std::tuple<PrecomputeCache<Derivatives>...>;

      //! Public output range type for the quantity selected by D.
      template <class D>
      using DerivativeRange = typename Basis::template DerivativeRange<D>;

      template <class D>
      using EvaluationBuffer = std::vector<std::vector<DerivativeRange<D>>>;

      //! Container type stored in the cache for the quantity selected by D.
      template <class D>
      struct EvaluationCache
      {
        EvaluationBuffer<D> buffers;
      };

      using EvaluationBuffers = std::tuple<EvaluationCache<Derivatives>...>;

    public:

      QuadratureNodeCache () = default;

      explicit QuadratureNodeCache (Node const& node, Derivatives const&... d)
        : node_(&node)
        , derivatives_(d...)
      {}

      Basis const& basis () const
      {
        assert(!!node_);
        return node_->finiteElement().physicalBasis();
      }

      //! Return the number of shape functions in the bound physical basis.
      std::size_t size() const
      {
        return basis().size();
      }

      //! Return the polynomial order of the bound physical basis.
      int order() const
      {
        return basis().order();
      }

      void bind(Node const& node)
      {
        node_ = &node;
      }

      template <class Derivative, class ct, int dim>
      auto const& evaluate (Derivative const& d, QuadratureRule<ct, dim> const& quad)
      {
        auto [precomputeBuffer, inserted] = getPrecomputeBuffer(d,quadratureKey(quad));
        if (inserted) {
          precomputeBuffer.resize(quad.size());
          for (std::size_t iq = 0; iq < quad.size(); ++iq)
            basis().precompute(d, quad[iq].position(), precomputeBuffer[iq]);
        }

        auto& evaluationBuffer = getEvaluationBuffer(d);
        evaluationBuffer.resize(quad.size());
        for (std::size_t iq = 0; iq < quad.size(); ++iq)
          basis().finalize(d, quad[iq].position(), precomputeBuffer[iq], evaluationBuffer[iq]);

        return evaluationBuffer;
      }

    private:

      template <class ct, int dim>
      QuadratureKey quadratureKey(QuadratureRule<ct,dim> const& quad) const
      {
        // QuadratureRule does not expose its QuadratureType. Use stable metadata
        // as an approximate identity. Weights are irrelevant because precompute()
        // only evaluates the basis.
        return {
          quad.type(),
          quad.order(),
          quad.size()
        };
      }

      template <class Derivative>
      std::pair<PrecomputeBuffer<Derivative>&,bool>
      getPrecomputeBuffer (Derivative const& d, QuadratureKey const& key)
      {
        auto& buffers = std::get<PrecomputeCache<Derivative>>(precomputeCaches_).buffers;
        auto [it,inserted] = buffers.try_emplace(key);
        return {it->second,inserted};
      }

      template <class Derivative>
      EvaluationBuffer<Derivative>& getEvaluationBuffer (Derivative const& d)
      {
        return std::get<EvaluationCache<Derivative>>(evaluationBuffers_).buffers;
      }

    private:

      Node const* node_ = nullptr;
      std::tuple<Derivatives...> derivatives_ = {};

      PrecomputeCaches precomputeCaches_;
      EvaluationBuffers evaluationBuffers_;
    };

    struct QuadratureNodeCacheFactory
    {
      explicit QuadratureNodeCacheFactory (std::tuple<Derivatives...> const& cachedDerivatives)
        : cachedDerivatives_(cachedDerivatives)
      {}

      template <class LeafNode>
      auto operator() (LeafNode const& leafNode) const
      {
        return std::apply([&](auto const&... d) {
          return QuadratureNodeCache<LeafNode>{leafNode,d...};
        }, cachedDerivatives_);
      }

      std::tuple<Derivatives...> cachedDerivatives_;
    };

  public:
    QuadratureBasisFunctionCache ()
      : cachedDerivatives_(Derivatives{}...)
    {}

    explicit QuadratureBasisFunctionCache (Derivatives const&... d)
      : cachedDerivatives_(d...)
    {}

    template <class... II>
    auto& get (TypeTree::TreePath<II...> const& tp)
    {
      assert(cache_);
      return (*cache_)[tp];
    }

    template <std::convertible_to<std::size_t>... II>
    auto& get (II const&... ii)
    {
      return get(TypeTree::treePath(ii...));
    }

    template <class LocalView>
    void bind (LocalView const& localView)
    {
      if (!cache_)
        cache_.emplace(TypeTree::makeTreeContainer(localView.tree(),
          QuadratureNodeCacheFactory{cachedDerivatives_}));

      TypeTree::forEachLeafNode(localView.tree(), [&](auto const& node, auto const& treePath) {
        (*cache_)[treePath].bind(node);
      });
    }

  private:
    std::tuple<Derivatives...> cachedDerivatives_;

    using Cache = TypeTree::TreeContainer<QuadratureNodeCache, Tree>;
    std::optional<Cache> cache_;
  };

} // end namespace Dune::Functions

#endif
