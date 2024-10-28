// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NODES_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_NODES_HH

#include <cassert>
#include <memory>

#include <dune/common/indices.hh>

#include <dune/typetree/leafnode.hh>
#include <dune/typetree/powernode.hh>
#include <dune/typetree/dynamicpowernode.hh>
#include <dune/typetree/compositenode.hh>
#include <dune/typetree/traversal.hh>
#include <dune/typetree/visitor.hh>

namespace Dune {
  namespace Functions {

    namespace Concept {

      struct HasElementBind
      {
        template<class N>
        auto require(N&& node) -> decltype(
          node.bind(std::declval<typename N::Element>())
          );
      };

      struct HasElementBind2
      {
        template<class N>
        auto require(N&& node) -> decltype(
          node.bind(std::declval<typename N::Element>(), std::declval<std::size_t&>())
          );
      };

    } // end namespace Concept

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
          node.setSize(0);
        }

        ClearSizeVisitor(std::size_t offset)
          : offset_(offset)
        {}

        const std::size_t offset_;

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

      template<typename Node, typename Element>
      void callNodeBind(Node& node, const Element& entity, std::size_t& offset)
      {
        using namespace Dune::Functions::Concept;
        static_assert(
          models<HasElementBind2,Node>()
          || (Node::isLeaf && models<HasElementBind,Node>()),
          "");

        // if the node does not implement the new 2-arg bind
        // *and* the node is a leaf node, we can do the leaf-bind directly here.
        if constexpr (Node::isLeaf && !models<HasElementBind2,Node>())
        {
          // we directly implement the leaf bind here, until the leaf
          // interface is updated to also handle the offset
          node.setOffset(offset);
          node.bind(entity);
          offset += node.size();
        }
        else {
          node.bind(entity, offset);
        }
      }

    } // end namespace Impl


    class BasisNodeMixin
    {

      friend struct Impl::ClearSizeVisitor;

      friend struct Impl::InitializeTreeVisitor;

      template<typename Node, typename Element>
      friend
      void Impl::callNodeBind(Node& node, const Element& entity, std::size_t& offset);

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
        public BasisNodeMixin,
        public TypeTree::LeafNode
    {
    public:
      constexpr bool active () const { return true; }
    };


    template<typename Node, typename Element>
    class BranchNodeMixin :
      public BasisNodeMixin
    {
    public:

      void bind(const Element& entity, std::size_t& offset)
      {
        // cast to actual implementation
        Node& node = *static_cast<Node*>(this);

        node.setOffset(offset);

        // iterate over child-nodes
        Dune::Hybrid::forEach(Dune::range(node.degree()), [&](auto i) {
          Impl::callNodeBind(node.child(i), entity, offset);
        });

        node.setSize(offset - node.offset());
      }

    };

    template<typename T, std::size_t n>
    class PowerBasisNode :
      public BranchNodeMixin<PowerBasisNode<T,n>,
                             typename T::Element>,
      public TypeTree::PowerNode<T,n>
    {

      using Node = TypeTree::PowerNode<T,n>;

    public:

      using Element = typename T::Element;

      PowerBasisNode() = default;

      PowerBasisNode(const typename Node::NodeStorage& children) :
        Node(children)
      {}

      const Element& element() const
      {
        return this->child(Dune::Indices::_0).element();
      }

    };


    template<typename T>
    class DynamicPowerBasisNode :
      public BranchNodeMixin<DynamicPowerBasisNode<T>,
                             typename T::Element>,
      public TypeTree::DynamicPowerNode<T>
    {

      using Node = TypeTree::DynamicPowerNode<T>;

    public:

      using Element = typename T::Element;

      DynamicPowerBasisNode (std::size_t children)
        : Node(children)
      {}

      DynamicPowerBasisNode (typename Node::NodeStorage children)
        : Node(std::move(children))
      {}

      const Element& element() const
      {
        return this->child(0).element();
      }

    };


    template<typename... T>
    class CompositeBasisNode :
      public BranchNodeMixin<CompositeBasisNode<T...>,
                             typename TypeTree::CompositeNode<T...>::template Child<0>::Type::Element>,
      public TypeTree::CompositeNode<T...>
    {

      using Node = TypeTree::CompositeNode<T...>;

    public:

      using Element = typename Node::template Child<0>::Type::Element;

      CompositeBasisNode() = default;

      CompositeBasisNode(const typename Node::NodeStorage& children) :
        Node(children)
      {}

      template<typename... Children>
      CompositeBasisNode(const std::shared_ptr<Children>&... children) :
        Node(children...)
      {}

      const Element& element() const
      {
        return this->child(Dune::Indices::_0).element();
      }

    };


    template<typename Tree>
    void clearSize(Tree& tree, std::size_t offset)
    {
      TypeTree::applyToTree(tree,Impl::ClearSizeVisitor(offset));
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
