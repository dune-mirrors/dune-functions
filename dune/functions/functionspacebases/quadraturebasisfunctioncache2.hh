// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_QUADRATUREBASISFUNCTIONCACHE2_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_QUADRATUREBASISFUNCTIONCACHE2_HH

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <memory>
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

#include <dune/functions/functionspacebases/transformed/precomputeidentity.hh>

namespace Dune::Functions
{
  /**
   * \brief Cache selected physical-basis evaluations at quadrature points.
   *
   * Each leaf cache exposes the bound physical basis and caches its reference
   * precomputations independently for every derivative tag, reference-basis
   * identity, and quadrature-rule object. The basis identity distinguishes, for
   * example, Raviart-Thomas orientation variants and hp finite elements. Rules
   * obtained from QuadratureRules::rule() have process-lifetime stable addresses
   * and are therefore ideal cache keys. Other rules must outlive the cache.
   *
   * Finalization is repeated after every bind and therefore always uses the
   * currently bound physical element. Call clear() after updating a basis unless
   * its reference finite elements provide generation-aware precompute identities.
   *
   * \b Example
   * \code{.cpp}
    auto localView = basis.localView();
    auto cache = QuadratureBasisFunctionCache2<LocalView::Tree, Derivatives::Value>{};
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
  class QuadratureBasisFunctionCache2
  {
    template <class Node>
    class QuadratureNodeCache
    {
      using Basis = typename Node::FiniteElement::PhysicalBasis;

      struct PrecomputeKey
      {
        LocalBasisPrecomputeIdentity basis;
        void const* quadratureRule = nullptr;

        bool operator==(PrecomputeKey const&) const = default;
      };

      //! Buffer type stored in the cache for the quantity selected by D.
      template <class D>
      using PrecomputeBuffer = std::vector<typename Basis::template PrecomputeBuffer<D>>;

      template<class D>
      struct PrecomputeCache
      {
        std::vector<std::pair<PrecomputeKey, PrecomputeBuffer<D>>> buffers;
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
        auto [precomputeBuffer, inserted] =
          getPrecomputeBuffer(d,precomputeKey(quad));
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
      PrecomputeKey precomputeKey(QuadratureRule<ct,dim> const& quad) const
      {
        return {
          node_->finiteElement().precomputeIdentity(),
          std::addressof(quad)
        };
      }

      template <class Derivative>
      std::pair<PrecomputeBuffer<Derivative>&,bool>
      getPrecomputeBuffer (Derivative const& d, PrecomputeKey const& key)
      {
        using Cache = PrecomputeCache<Derivative>;
        using Buffer = PrecomputeBuffer<Derivative>;
        auto& buffers = std::get<Cache>(precomputeCaches_).buffers;
        if (auto it = std::find_if(buffers.begin(), buffers.end(),
              [&](auto& v) { return v.first == key; }); it != buffers.end())
          return {it->second, false};
        else
          return {buffers.emplace_back(key,Buffer{}).second, true};
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
    QuadratureBasisFunctionCache2 ()
      : cachedDerivatives_(Derivatives{}...)
    {}

    explicit QuadratureBasisFunctionCache2 (Derivatives const&... d)
      : cachedDerivatives_(d...)
    {}

    /**
     * \brief Remove all cached precomputations and bound-node references.
     *
     * Call this after updating the global basis if its reference finite-element
     * identities are based on object addresses without generation tracking.
     */
    void clear()
    {
      cache_.reset();
    }

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
