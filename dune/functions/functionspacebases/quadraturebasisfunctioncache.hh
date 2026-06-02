// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_QUADRATUREBASISFUNCTIONCACHE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_QUADRATUREBASISFUNCTIONCACHE_HH

#include <concepts>
#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include <dune/common/overloadset.hh>
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
   * basis functions on the leaf nodes of that tree. The tree cache is initialized at
   * evaluation at quadrature points.
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
  template <class Tree, class... CachedDerivatives>
  class QuadratureBasisFunctionCache
  {
    template <class Node>
    class QuadratureNodeCache
    {
      using Basis = typename Node::FiniteElement::PhysicalBasis;

      //! Buffer type stored in the cache for the quantity selected by D.
      template <class D>
      using PrecomputeBuffer = std::vector<typename Basis::template PrecomputeBuffer<D>>;

      using PrecomputeBuffers = std::tuple<PrecomputeBuffer<CachedDerivatives>...>;

      //! Public output range type for the quantity selected by D.
      template <class D>
      using DerivativeRange = typename Basis::template DerivativeRange<D>;

      //! Container type stored in the cache for the quantity selected by D.
      template <class D>
      using EvaluationBuffer = std::vector<std::vector<DerivativeRange<D>>>;

      using EvaluationBuffers = std::tuple<EvaluationBuffer<CachedDerivatives>...>;

    public:

      QuadratureNodeCache () = default;

      explicit QuadratureNodeCache (Node const& node, CachedDerivatives const&... d)
        : node_(&node)
        , derivatives_(d...)
      {}

      Basis const& basis () const
      {
        assert(!!node_);
        return node_->finiteElement().physicalBasis();
      }

      template <class Derivative, class ct, int dim>
      auto const& evaluate (Derivative const& d, QuadratureRule<ct, dim> const& quad)
      {
        auto& precomputeBuffer = getPrecomputeBuffer(d);
        if (precomputeBuffer.size() != quad.size())
        {
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

      template <class Derivative>
      PrecomputeBuffer<Derivative>& getPrecomputeBuffer (Derivative const& d)
      {
        PrecomputeBuffer<Derivative>* result = nullptr;
        auto setResult = overload(
          [&](Derivative const& d_ii, PrecomputeBuffer<Derivative>& buffer_ii) { d_ii == d && (result = &buffer_ii); },
          [&](auto const& d_ii, auto& buffer_ii) {}
        );
        unpackIntegerSequence([&](auto... ii) {
          (setResult(std::get<ii>(derivatives_), std::get<ii>(precomputeBuffers_)), ...);
        }, std::index_sequence_for<CachedDerivatives...>{});

        assert(!!result);
        return *result;
      }

      template <class Derivative>
      EvaluationBuffer<Derivative>& getEvaluationBuffer (Derivative const& d)
      {
        EvaluationBuffer<Derivative>* result = nullptr;
        auto setResult = overload(
          [&](Derivative const& d_ii, EvaluationBuffer<Derivative>& buffer_ii) { d_ii == d && (result = &buffer_ii); },
          [&](auto const& d_ii, auto& buffer_ii) {}
        );
        unpackIntegerSequence([&](auto... ii) {
          (setResult(std::get<ii>(derivatives_), std::get<ii>(evaluationBuffers_)), ...);
        }, std::index_sequence_for<CachedDerivatives...>{});

        assert(!!result);
        return *result;
      }

    private:

      Node const* node_ = nullptr;
      std::tuple<CachedDerivatives...> derivatives_ = {};

      PrecomputeBuffers precomputeBuffers_;
      EvaluationBuffers evaluationBuffers_;
    };

    struct QuadratureNodeCacheFactory
    {
      explicit QuadratureNodeCacheFactory (std::tuple<CachedDerivatives...> const& cachedDerivatives)
        : cachedDerivatives_(cachedDerivatives)
      {}

      template <class LeafNode>
      auto operator() (LeafNode const& leafNode) const
      {
        return std::apply([&](auto const&... d) {
          return QuadratureNodeCache<LeafNode>{leafNode,d...};
        }, cachedDerivatives_);
      }

      std::tuple<CachedDerivatives...> cachedDerivatives_;
    };

  public:
    QuadratureBasisFunctionCache ()
      : cachedDerivatives_(CachedDerivatives{}...)
    {}

    explicit QuadratureBasisFunctionCache (CachedDerivatives const&... d)
      : cachedDerivatives_(d...)
    {}

    template <class... II>
    auto& get (TypeTree::TreePath<II...> const& tp)
    {
      assert(cacheIterator_ != cache_.end());
      return (cacheIterator_->second)[tp];
    }

    template <std::convertible_to<std::size_t>... II>
    auto& get (II const&... ii)
    {
      return get(TypeTree::treePath(ii...));
    }

    template <class LocalView>
    void bind (LocalView const& localView)
    {
      auto key = localView.element().type();
      auto it = cache_.find(key);
      if (it == cache_.end())
      {
        bool inserted = false;
        std::tie(it, inserted) = cache_.emplace(key,
          TypeTree::makeTreeContainer(localView.tree(),
            QuadratureNodeCacheFactory{cachedDerivatives_}));
        assert(inserted);
      }
      cacheIterator_ = it;
    }

  private:
    std::tuple<CachedDerivatives...> cachedDerivatives_;

    using Key = GeometryType;
    using Cache = TypeTree::TreeContainer<QuadratureNodeCache, Tree>;
    std::map<Key,Cache> cache_ = {};

    typename std::map<Key,Cache>::iterator cacheIterator_ = {};
  };

} // end namespace Dune::Functions

#endif
