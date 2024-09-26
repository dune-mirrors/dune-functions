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

      // forward declaration
      template<typename Node, typename Entity>
      void bindTreeForward(Node& node, const Entity& entity, std::size_t& offset);

      // generic bind, i.e. propagte size information and collect offsets
      // should match any non-specialized tag,
      // i.e. PowerNodeTag, DynamicPowerNodeTag and CompositeNodeTag
      template<typename Node, typename Entity, typename Tag>
      void bindTree(Tag, Node& node, const Entity& entity, std::size_t& offset)
      {
          node.setOffset(offset);

          // iterate over child-nodes
          Dune::Hybrid::forEach(Dune::range(node.degree()), [&](auto i) {
            bindTreeForward(node.child(i), entity, offset);
          });

          node.setSize(offset - node.offset());
      }

      // specialization for leaf node
      template<typename Node, typename Entity>
      void bindTree(TypeTree::LeafNodeTag, Node& node, const Entity& entity, std::size_t& offset)
      {
        node.setOffset(offset);
        node.bind(entity);
        offset += node.size();
      }

      // special handling for restricted nodes
      // template<typename Tree, typename Entity>
      // void bindTree(RestrictedNodeTag, Tree& tree, const Entity& entity, std::size_t offset = 0)
      // {
      // }

      // taking a node we extract the tag and forward the call using a tag dispatch
      template<typename Node, typename Entity>
      void bindTreeForward(Node& node, const Entity& entity, std::size_t& offset)
      {
        // static_assert(hasNodeTag_v<Tree>, "tree is not an instance of a typetree node");
        bindTree(typename Node::NodeTag{}, node, entity, offset);
      }

    } // end namespace Impl


    class BasisNodeMixin
    {

      friend struct Impl::ClearSizeVisitor;

      // template<typename>
      // friend struct Impl::BindVisitor;

      friend struct Impl::InitializeTreeVisitor;

      template<typename Node, typename Entity, typename Tag>
      friend void Impl::bindTree(Tag, Node&, const Entity&, std::size_t&);

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

      // friend declaration does not work ... temporarily make protected methods public
      // protected:

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
    {};


    template<typename T, std::size_t n>
    class PowerBasisNode :
      public BasisNodeMixin,
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
      public BasisNodeMixin,
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
      public BasisNodeMixin,
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

    template<typename Tree, typename Entity>
    void bindTree(Tree& tree, const Entity& entity, std::size_t offset = 0)
    {
      Impl::bindTreeForward(tree, entity, offset);
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
