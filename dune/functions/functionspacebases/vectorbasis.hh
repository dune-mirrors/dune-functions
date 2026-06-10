// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_BASIS_VECTORBASIS_HH
#define DUNE_FUNCTIONS_BASIS_VECTORBASIS_HH

#include <dune/functions/common/utility.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/containerdescriptors.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/concepts.hh>

#include <dune/common/typetree/nodeconcepts.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/typeutilities.hh>
#include <dune/common/indices.hh>

namespace Dune {
namespace Functions {


// *****************************************************************************
// This is the reusable part of the dynamic vector bases. It contains
//
//   VectorPreBasis
//
// The pre-basis allows to create the others and is the owner of possible shared
// state. These components do _not_ depend on the global basis and local view
// and can be used without a global basis.
// *****************************************************************************

/**
 * \brief A pre-basis for dynamic vector bases
 *
 * This pre-basis represents a vector of pre-bases.
 * Its node type is a DynamicPowerBasisNodes for the given subnode.
 *
 * \tparam IMS An IndexMergingStrategy used to merge the global indices of the child factories
 * \tparam SPB  The child pre-basis
 */
template<class IMS, class SPB>
class VectorPreBasis
{
  static const bool isBlocked = std::is_same_v<IMS,Functions::BasisFactory::BlockedLexicographic>;

public:

  //! The child pre-basis
  using SubPreBasis = SPB;

  //! The grid view that the FE basis is defined on
  using GridView = typename SPB::GridView;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Strategy used to merge the global indices of the child factories
  using IndexMergingStrategy = IMS;

  //! Template mapping root tree path to type of created tree node
  using Node = Functions::DynamicPowerBasisNode<typename SubPreBasis::Node>;

  static constexpr size_type maxMultiIndexSize = SubPreBasis::maxMultiIndexSize + isBlocked;
  static constexpr size_type minMultiIndexSize = SubPreBasis::minMultiIndexSize + isBlocked;
  static constexpr size_type multiIndexBufferSize = SubPreBasis::multiIndexBufferSize + isBlocked;

  //! Constructor for a collection of child pre-basis objects
  VectorPreBasis(const std::vector<SubPreBasis>& subPreBases)
  {
    for (const auto& spb : subPreBases)
      subPreBases_.emplace_back(spb);
  }

  //! Initialize the global indices
  void initializeIndices()
  {
    for (auto& spb : subPreBases_)
      spb->initializeIndices();
  }

  //! Obtain the grid view that the basis is defined on
  const GridView& gridView() const
  {
    for (auto& spb : subPreBases_)
      return spb->gridView();
    DUNE_THROW(InvalidStateException, "No subPreBasis available to get grid view from.");
  }

  //! Update the stored grid view, to be called if the grid has changed
  void update(const GridView& gv)
  {
    for (auto& spb : subPreBases_)
      spb->update(gv);
  }

  /**
   * \brief Create tree node
   */
  Node makeNode() const
  {
    auto node = Node{subPreBases_.size()};
    for (std::size_t i=0; i<subPreBases_.size(); ++i)
      node.setChild(i, subPreBases_[i]->makeNode());
    return node;
  }

  //! Get the number of child pre-bases
  std::size_t children() const
  {
    return subPreBases_.size();
  }

  //! Same as size(prefix) with empty prefix
  size_type size() const
  {
    return size(Dune::ReservedVector<size_type, multiIndexBufferSize>{});
  }

  //! Return number of possible values for next position in multi index
  template<class SizePrefix>
  size_type size(const SizePrefix& prefix) const
  {
    return sizeImpl(prefix, IndexMergingStrategy{});
  }

protected:

  template<class SizePrefix>
  size_type sizeImpl(SizePrefix prefix, Functions::BasisFactory::FlatLexicographic) const
  {
    if (prefix.size() == 0) {
      std::size_t s = 0;
      for (const auto& spb : subPreBases_)
        s += spb->size();
      return s;
    }

    for (std::size_t i = 0; i < subPreBases_.size(); ++i) {
      if (prefix[0] < subPreBases_[i]->size())
        return subPreBases_[i]->size(prefix);
      prefix[0] -= subPreBases_[i]->size();
    }
    DUNE_THROW(RangeError, "Prefix index out of range in VectorPreBasis::size().");
  }

  template<class MultiIndex>
  static void multiIndexPopFront(MultiIndex& M)
  {
    for(std::size_t i=0; i<M.size()-1; ++i)
      M[i] = M[i+1];
    M.resize(M.size()-1);
  }

  template<class SizePrefix>
  size_type sizeImpl(SizePrefix prefix, Functions::BasisFactory::BlockedLexicographic) const
  {
    if (prefix.size() == 0)
      return subPreBases_.size();
    auto child = prefix[0];
    multiIndexPopFront(prefix);
    return subPreBases_[child]->size(prefix);
  }

public:

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    std::size_t dim = 0;
    for (const auto& spb : subPreBases_)
      dim += spb->dimension();
    return dim;
  }

  //! Get the maximal number of DOFs associated to node for any element
  size_type maxNodeSize() const
  {
    std::size_t mnz = 0;
    for (const auto& spb : subPreBases_)
      mnz += spb->maxNodeSize();
    return mnz;
  }

  //! Const access to the stored prebasis of the factor in the power space
  const SubPreBasis& subPreBasis(std::size_t i) const
  {
    return *subPreBases_[i];
  }

  //! Mutable access to the stored prebasis of the factor in the power space
  SubPreBasis& subPreBasis(std::size_t i)
  {
    return *subPreBases_[i];
  }

  //! Maps from subtree index set [0..size-1] to a globally unique multi index in global basis
  template<class NodeType, typename It>
  requires Dune::TypeTree::Concept::UniformInnerTreeNode<NodeType>
  It indices(const NodeType& node, It it) const
  {
    return indicesImpl(node, it, IndexMergingStrategy{});
  }

  //! Return the associated container descriptor
  auto containerDescriptor() const
  {
    return containerDescriptorImpl(subPreBases_.size());
  }

protected:

  template<class NodeType, typename It>
  It indicesImpl(const NodeType& node, It multiIndices, Functions::BasisFactory::FlatLexicographic) const
  {
    std::size_t sb_offset = 0;
    for (std::size_t child = 0; child<subPreBases_.size(); ++child)
    {
      auto next = subPreBasis(child).indices(node.child(child), multiIndices);
      for (std::size_t i = 0; i<node.child(child).size(); ++i)
        multiIndices[i][0] += sb_offset;
      sb_offset += subPreBasis(child).size();
      multiIndices = next;
    }
    return multiIndices;
  }

  template<class MultiIndex>
  static void multiIndexPushFront(MultiIndex& M, size_type M0)
  {
    M.resize(M.size()+1);
    for(std::size_t i=M.size()-1; i>0; --i)
      M[i] = M[i-1];
    M[0] = M0;
  }

  template<class NodeType, typename It>
  It indicesImpl(const NodeType& node, It multiIndices, Functions::BasisFactory::BlockedLexicographic) const
  {
    for (std::size_t child = 0; child<subPreBases_.size(); ++child)
    {
      auto next = subPreBasis(child).indices(node.child(child), multiIndices);
      for (std::size_t i = 0; i<node.child(child).size(); ++i)
        multiIndexPushFront(multiIndices[i], child);
      multiIndices = next;
    }
    return multiIndices;
  }

  template<class Children>
  auto containerDescriptorImpl(Children children) const
  {
    auto subTree0 = Dune::Functions::containerDescriptor(*subPreBases_[0]);
    using SubTree = decltype(subTree0);
    if constexpr (std::same_as<SubTree, ContainerDescriptors::Unknown>)
      return ContainerDescriptors::Unknown{};
    else if constexpr(std::is_same_v<IMS, Functions::BasisFactory::FlatLexicographic>) {
      using SubSubTree = typename ContainerDescriptors::Impl::BlockType<SubTree>::type;
      if constexpr (std::same_as<SubTree, ContainerDescriptors::Unknown>)
        return ContainerDescriptors::Unknown{};
      else {
        ContainerDescriptors::Vector<SubSubTree> vector;
        for (std::size_t i = 0; i != subPreBases_.size(); ++i) {
          auto desc = Dune::Functions::containerDescriptor(*subPreBases_[i]);
          Hybrid::forEach(range(Hybrid::size(desc)), [&](auto j){
            vector.emplace_back(std::move(desc[j]));
          });
        }
        return vector;
      }
    } else if constexpr(std::is_same_v<IMS, Functions::BasisFactory::BlockedLexicographic>) {
      ContainerDescriptors::Vector<SubTree> vector;
      for (std::size_t i = 0; i != subPreBases_.size(); ++i)
        vector.emplace_back(Dune::Functions::containerDescriptor(*subPreBases_[i]));
      return vector;
    }
    else
      return ContainerDescriptors::Unknown{};
  }

protected:
  std::vector<std::optional<SubPreBasis>> subPreBases_;
};

} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_BASIS_VECTORBASIS_HH
