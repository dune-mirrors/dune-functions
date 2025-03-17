// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NODES_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NODES_HH

#include <cassert>
#include <memory>
#include <vector>
#include <array>
#include <optional>

#include <dune/common/indices.hh>
#include <dune/common/tuplevector.hh>

#include <dune/typetree/nodetags.hh>
#include <dune/typetree/traversal.hh>
#include <dune/typetree/visitor.hh>

namespace Dune {
  namespace Functions {


    namespace Impl {


      struct ClearSizeVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename Node, typename TreePath>
        void pre(Node& node, TreePath treePath)
        {
          leaf(node,treePath);
          node.setSize(0);
        }

        template<typename Node, typename TreePath>
        void leaf(Node& node, TreePath treePath)
        {
          node.setOffset(offset_);
        }

        ClearSizeVisitor(std::size_t offset)
          : offset_(offset)
        {}

        const std::size_t offset_;

      };


      template<typename Entity>
      struct BindVisitor
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename Node, typename TreePath>
        void pre(Node& node, TreePath)
        {
          node.setOffset(offset_);
        }

        template<typename Node, typename TreePath>
        void post(Node& node, TreePath)
        {
          node.setSize(offset_ - node.offset());
        }

        template<typename Node, typename TreePath>
        void leaf(Node& node, TreePath)
        {
          node.setOffset(offset_);
          node.bind(entity_);
          offset_ += node.size();
        }

        BindVisitor(const Entity& entity, std::size_t offset = 0)
          : entity_(entity)
          , offset_(offset)
        {}

        const Entity& entity_;
        std::size_t offset_;

      };


      struct InitializeTreeVisitor :
        public TypeTree::TreeVisitor,
        public TypeTree::DynamicTraversal
      {
        template<typename Node, typename TreePath>
        void pre(Node& node, TreePath)
        {
          node.setTreeIndex(treeIndex_);
          ++treeIndex_;
        }

        template<typename Node, typename TreePath>
        void leaf(Node& node, TreePath)
        {
          node.setTreeIndex(treeIndex_);
          ++treeIndex_;
        }

        InitializeTreeVisitor(std::size_t treeIndexOffset = 0) :
          treeIndex_(treeIndexOffset)
        {}

        std::size_t treeIndex_;
      };

    } // end namespace Impl


    class BasisNodeMixin
    {

      friend struct Impl::ClearSizeVisitor;

      template<typename>
      friend struct Impl::BindVisitor;

      friend struct Impl::InitializeTreeVisitor;

    public:

      using size_type = std::size_t;

      BasisNodeMixin() :
        offset_(0),
        size_(0),
        treeIndex_(0)
      {}

      size_type localIndex(size_type i) const
      {
        assert(i < size_);
        return offset_ + i;
      }

      size_type size() const
      {
        return size_;
      }

      size_type treeIndex() const
      {
        return treeIndex_;
      }

    protected:

      size_type offset() const
      {
        return offset_;
      }

      void setOffset(const size_type offset)
      {
        offset_ = offset;
      }

      void setSize(const size_type size)
      {
        size_ = size;
      }

      void setTreeIndex(size_type treeIndex)
      {
        treeIndex_ = treeIndex;
      }

    private:

      size_type offset_;
      size_type size_;
      size_type treeIndex_;

    };


    class LeafBasisNode :
        public BasisNodeMixin
    {
    public:

      // Begin of node interface

      static const bool isLeaf = true;
      static const bool isPower = false;
      static const bool isComposite = false;
      using NodeTag = Dune::TypeTree::LeafNodeTag;

      static constexpr auto degree()
      {
        return Dune::index_constant<0>{};
      }

      // End of node interface

    };



    // A mixin class for generalized child access from
    // multiple indices or a tree path. The derived class
    // only has to provide the child(i) method with
    // a single index for accessing direct children.
    template<class Impl>
    class ChildAccessMixIn
    {

      Impl& asImpl()
      {
        return static_cast<Impl&>(*this);
      }

      const Impl& asImpl() const
      {
        return static_cast<const Impl&>(*this);
      }

    public:

      // Access with tree path
      template<class I>
      const auto& child(I i) const
      requires (not std::is_convertible_v<I, std::size_t>)
      {
        if constexpr (I::size()==0)
          return asImpl();
        else if constexpr (I::size()==1)
          return asImpl().child(i.front());
        else
          return asImpl().child(i.front()).child(pop_front(i));
      }

      template<class I>
      auto& child(I i)
      requires (not std::is_convertible_v<I, std::size_t>)
      {
        if constexpr (I::size()==0)
          return asImpl();
        else if constexpr (I::size()==1)
          return asImpl().child(i.front());
        else
          return asImpl().child(i.front()).child(pop_front(i));
      }

      // Access without path
      const auto& child() const
      {
        return asImpl();
      }

      auto& child()
      {
        return asImpl();
      }

      // Access without path longer than 1
      template<class I0, class... II>
      const auto& child(I0 i0, II... ii) const
      requires (sizeof...(II)>0)
      {
        return asImpl().child(i0).child(ii...);
      }

      // Access without path longer than 1
      template<class I0, class... II>
      auto& child(I0 i0, II... ii)
      requires (sizeof...(II)>0)
      {
        return asImpl().child(i0).child(ii...);
      }

    };



    template<class Node, class I>
    requires (not std::is_convertible_v<I, std::size_t>)
    auto& child(Node&& node, I i)
    {
      if constexpr(I::size()==0)
        return node;
      else
        return node.child(i);
    }

    template<class Node, class... II>
    requires ((std::is_convertible_v<II, std::size_t> && ...))
    auto& child(Node&& node, II... ii)
    {
      if constexpr(sizeof...(II)==0)
        return node;
      else
        return node.child(ii...);
    }


    template<typename T, std::size_t n>
    class PowerBasisNode
      : public BasisNodeMixin
      , public ChildAccessMixIn<PowerBasisNode<T, n>>
    {
    public:

      // Begin of node interface

      static const bool isLeaf = false;
      static const bool isPower = true;
      static const bool isComposite = false;
      using NodeTag = Dune::TypeTree::PowerNodeTag;

      using ChildType = T;

      static constexpr auto degree()
      {
        return Dune::index_constant<n>{};
      }

      template<class Index>
      requires (std::is_convertible_v<Index, std::size_t>)
      const auto& child(Index i) const
      {
        return children_[i].value();
      }

      template<class Index>
      requires (std::is_convertible_v<Index, std::size_t>)
      auto& child(Index i)
      {
        return children_[i].value();
      }

      using ChildAccessMixIn<PowerBasisNode<T, n>>::child;

      // End of node interface

      using Element = typename T::Element;

      PowerBasisNode() = default;

      const Element& element() const
      {
        return child(Dune::Indices::_0).element();
      }

      template<class Index, class TT>
      void setChild(Index i, TT&& t)
      {
        children_[i].emplace(std::forward<TT>(t));
      }

    private:
      std::array<std::optional<T>, n> children_;
    };


    template<typename T>
    class DynamicPowerBasisNode
      : public BasisNodeMixin
      , public ChildAccessMixIn<DynamicPowerBasisNode<T>>
    {
    public:

      // Begin of node interface

      static const bool isLeaf = false;
      static const bool isPower = true;
      static const bool isComposite = false;
      using NodeTag = Dune::TypeTree::DynamicPowerNodeTag;

      using ChildType = T;

      std::size_t degree() const
      {
        return children_.size();
      }

      template<class Index>
      requires (std::is_convertible_v<Index, std::size_t>)
      const auto& child(Index i) const
      {
        return children_[i].value();
      }

      template<class Index>
      requires (std::is_convertible_v<Index, std::size_t>)
      auto& child(Index i)
      {
        return children_[i].value();
      }

      using ChildAccessMixIn<DynamicPowerBasisNode<T>>::child;

      // End of node interface

      using Element = typename T::Element;

      DynamicPowerBasisNode (std::size_t children)
        : children_(children)
      {}

      const Element& element() const
      {
        return child(Dune::Indices::_0).element();
      }

      template<class Index, class TT>
      void setChild(Index i, TT&& t)
      {
        children_[i].emplace(std::forward<TT>(t));
      }

    private:
      std::vector<std::optional<T>> children_;
    };


    template<typename... T>
    class CompositeBasisNode
      : public BasisNodeMixin
      , public ChildAccessMixIn<CompositeBasisNode<T...>>
    {
    public:

      // Begin of node interface

      static const bool isLeaf = false;
      static const bool isPower = false;
      static const bool isComposite = true;
      using NodeTag = Dune::TypeTree::CompositeNodeTag;

      static constexpr auto degree()
      {
        return Dune::index_constant<sizeof...(T)>{};
      }

      using ChildTypes = std::tuple<T...>;

      template<std::size_t k>
      struct Child {
        static_assert((k < degree()), "child index out of range");

        //! The type of the child.
        using Type = typename std::tuple_element_t<k,ChildTypes>;

        using type = Type;
      };

      template<std::size_t i>
      const auto& child(Dune::index_constant<i> ii) const
      {
        return children_[ii].value();
      }

      template<std::size_t i>
      auto& child(Dune::index_constant<i> ii)
      {
        return children_[ii].value();
      }

      using ChildAccessMixIn<CompositeBasisNode<T...>>::child;

      // End of node interface

      using Element = typename Child<0>::Type::Element;

      CompositeBasisNode() = default;

      const Element& element() const
      {
        return child(Dune::Indices::_0).element();
      }

      template<std::size_t i, class TT>
      void setChild (TT&& t, Dune::index_constant<i> ii = {})
      {
        children_[ii].emplace(std::forward<TT>(t));
      }

    private:
      Dune::TupleVector<std::optional<T>...> children_;
    };


    template<typename Tree>
    void clearSize(Tree& tree, std::size_t offset)
    {
      TypeTree::applyToTree(tree,Impl::ClearSizeVisitor(offset));
    }

    template<typename Tree, typename Entity>
    void bindTree(Tree& tree, const Entity& entity, std::size_t offset = 0)
    {
      Impl::BindVisitor<Entity> visitor(entity,offset);
      TypeTree::applyToTree(tree,visitor);
    }

    template<typename Tree>
    void initializeTree(Tree& tree, std::size_t treeIndexOffset = 0)
    {
      Impl::InitializeTreeVisitor visitor(treeIndexOffset);
      TypeTree::applyToTree(tree,visitor);
    }


  } // namespace Functions

} // namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NODES_HH
