// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PRISMATICPRODUCTBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PRISMATICPRODUCTBASIS_HH

#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dune/localfunctions/meta/product.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/leafprebasismixin.hh>

namespace Dune::Functions {

// forward declaration
template <class N1, class N2>
class PrismaticProductNode;

template<class GB>
class PrismaticProductLocalView;

template<class PB1, class PB2>
class PrismaticProductGlobalBasis;

/**
 * \brief A pre-basis for a prismatic product of two prebases
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam PB1  Type of the first prebasis
 * \tparam PB2  Type of the second prebasis
 */
template<class PB1, class PB2>
class PrismaticProductPreBasis
    : public LeafPreBasisMixin<PrismaticProductPreBasis<PB1,PB2>>
{
public:
  using PreBasis1 = PB1;
  using PreBasis2 = PB2;

  // using GridView = PrismaticProductGridView<typename PreBasis1::GridView1, typename PreBasis2::GridView2>;
  using GridView1 = typename PreBasis1::GridView;
  using GridView2 = typename PreBasis2::GridView;

  using Node = PrismaticProductNode<typename PreBasis1::Node, typename PreBasis2::Node>;
  using size_type = std::common_type_t<typename PreBasis1::size_type, typename PreBasis2::size_type>;

  // using MI1 = typename DefaultGlobalBasis<PreBasis1>::MultiIndex;
  // using MI2 = typename DefaultGlobalBasis<PreBasis2>::MultiIndex;

  template<class>
  friend class PrismaticProductLocalView;

  template<class, class>
  friend class PrismaticProductGlobalBasis;

public:
  PrismaticProductPreBasis (const PreBasis1& preBasis1, const PreBasis2& preBasis2)
    : preBasis1_(preBasis1)
    , preBasis2_(preBasis2)
    // , indices1_(preBasis1_.maxNodeSize())
    // , indices2_(preBasis2_.maxNodeSize())
  {}

  //! Export the stored GridView
  const GridView1& gridView1() const
  {
    return preBasis1_.gridView();
  }

  const GridView2& gridView2() const
  {
    return preBasis2_.gridView();
  }

  //! Update the stored GridView
  void update(const GridView1& gv1, const GridView2& gv2)
  {
    preBasis1_.update(gv1);
    preBasis2_.update(gv2);
    // indices1_.resize(preBasis1_.maxNodeSize());
    // indices2_.resize(preBasis2_.maxNodeSize());
  }

  void initializeIndices()
  {
    preBasis1_.initializeIndices();
    preBasis2_.initializeIndices();
  }

  Node makeNode() const
  {
    return Node{preBasis1_.makeNode(), preBasis2_.makeNode()};
  }

  //! Return total number of basis functions
  size_type dimension() const
  {
    return preBasis1_.dimension() * preBasis2_.dimension();
  }

  //! Return maximal number of basis functions per element
  size_type maxNodeSize() const
  {
    return preBasis1_.maxNodeSize() * preBasis2_.maxNodeSize();
  }

  //! Fill cache with global indices of DOFs associated to the given bound node
  // template<class Node, class It>
  // It indices(const Node& node, It it) const
  // {
  //   const auto& node1 = node.subNode1();
  //   const auto& node2 = node.subNode2();

  //   indices1_.resize(node1.size());
  //   indices2_.resize(node2.size());
  //   preBasis1_.indices(node1, indices1_.begin());
  //   preBasis2_.indices(node2, indices2_.begin());

  //   for (const auto& idx1 : indices1_) {
  //     for (const auto& idx2 : indices2_) {
  //       *it++ = {{ idx1 * preBasis2_.dimension() + idx2 }};
  //     }
  //   }
  //   return it;
  // }

protected:
  PreBasis1 preBasis1_;
  PreBasis2 preBasis2_;
  // mutable std::vector<MI1> indices1_;
  // mutable std::vector<MI2> indices2_;
};


template <class N1, class N2>
class PrismaticProductNode
    : public LeafBasisNode
{
  static_assert(N1::isLeaf && N2::isLeaf);

public:
  using Node1 = N1;
  using Node2 = N2;

  using size_type = std::common_type_t<typename Node1::size_type, typename Node2::size_type>;
  using Element1 = typename Node1::Element;
  using Element2 = typename Node2::Element;

  using LFE1 = typename Node1::FiniteElement;
  using LFE2 = typename Node2::FiniteElement;
  using FiniteElement = Dune::PrismaticProduct<LFE1,LFE2>;

  template<class>
  friend class PrismaticProductLocalView;

public:
  //! Constructor from two individual prebasis nodes
  PrismaticProductNode (const Node1& node1, const Node2& node2)
    : node1_(node1)
    , node2_(node2)
  {}

  //! Return current element, throw if unbound
  auto element() const
  {
    return std::tie(node1_.element(), node2_.element());
  }

  //! Return the LocalFiniteElement for the element we are bound to
  const FiniteElement& finiteElement() const
  {
    assert(!!finiteElement_);
    return *finiteElement_;
  }

  void bind1(const Element1& e1)
  {
    node1_.bind(e1);
  }

  void bind2(const Element2& e2)
  {
    node2_.bind(e2);
  }

  void setup()
  {
    finiteElement_.emplace(node1_.finiteElement(), node2_.finiteElement());
    this->setSize(finiteElement_->size());
  }

  //! Bind to element.
  void bind(const Element1& e1, const Element2& e2)
  {
    bind1(e1);
    bind2(e2);
    setup();
  }

  template <class Element>
  void bind(const Element& e)
  {
    bind(std::get<0>(e), std::get<1>(e));
  }

  // const Node1& subNode1() const { return node1_; }
  // const Node2& subNode2() const { return node2_; }

private:
  Node1 node1_;
  Node2 node2_;
  std::optional<FiniteElement> finiteElement_ = std::nullopt;
};


template<class PB1, class PB2>
class PrismaticProductGlobalBasis
{
public:
  using PreBasis1 = PB1;
  using PreBasis2 = PB2;

  using PreBasis = PrismaticProductPreBasis<PreBasis1, PreBasis2>;

  //! The empty prefix path that identifies the root in the local ansatz tree
  using PrefixPath = TypeTree::HybridTreePath<>;

  //! The grid view that the FE space is defined on
  using GridView1 = typename PreBasis1::GridView;
  using GridView2 = typename PreBasis2::GridView;

  //! Type used for indices and size information
  using size_type = std::size_t;

  //! Type of the local view on the restriction of the basis to a single element
  using LocalView = PrismaticProductLocalView<PrismaticProductGlobalBasis>;

  //! Type used for global numbering of the basis vectors
  using MultiIndex = typename LocalView::MultiIndex;

  //! Type used for prefixes handed to the size() method
  using SizePrefix = Dune::ReservedVector<std::size_t, PreBasis::multiIndexBufferSize>;

public:
  PrismaticProductGlobalBasis(const GridView1& gridView1, const GridView2& gridView2)
    : PrismaticProductGlobalBasis(PreBasis1{gridView1}, PreBasis2{gridView2})
  {}

  PrismaticProductGlobalBasis(const PreBasis1& preBasis1, const PreBasis2& preBasis2)
    : preBasis_(preBasis1, preBasis2)
    , prefixPath_()
  {
    // static_assert(models<Concept::PreBasis<GridView>, PreBasis>(), "Type passed to DefaultGlobalBasis does not model the PreBasis concept.");
    preBasis_.initializeIndices();
  }

  //! Obtain the grid view that the basis is defined on
  const GridView1& gridView1() const
  {
    return preBasis_.gridView1();
  }

  //! Obtain the grid view that the basis is defined on
  const GridView2& gridView2() const
  {
    return preBasis_.gridView2();
  }

  //! Obtain the pre-basis providing the implementation details
  const PreBasis& preBasis() const
  {
    return preBasis_;
  }

  //! Obtain the pre-basis providing the implementation details
  PreBasis& preBasis()
  {
    return preBasis_;
  }

  const PreBasis1& preBasis1() const { return preBasis_.preBasis1_; }
  PreBasis1& preBasis1() { return preBasis_.preBasis1_; }

  const PreBasis2& preBasis2() const { return preBasis_.preBasis2_; }
  PreBasis2& preBasis2() { return preBasis_.preBasis2_; }

  /**
   * \brief Update the stored grid view
   *
   * This will update the indexing information of the global basis.
   * It must be called if the grid has changed.
   */
  void update(const GridView1& gv1, const GridView2& gv2)
  {
    preBasis_.update(gv1, gv2);
    preBasis_.initializeIndices();
  }

  //! Get the total dimension of the space spanned by this basis
  size_type dimension() const
  {
    return preBasis_.dimension();
  }

  //! Return number of possible values for next position in empty multi index
  size_type size() const
  {
    return preBasis_.size();
  }

  //! Return number of possible values for next position in multi index
  size_type size(const SizePrefix& prefix) const
  {
    return preBasis_.size(prefix);
  }

  //! Return local view for basis
  LocalView localView() const
  {
    return LocalView(*this);
  }

  //! Return *this because we are not embedded in a larger basis
  const PrismaticProductGlobalBasis& rootBasis() const
  {
    return *this;
  }

  //! Return empty path, because this is the root in the local ansatz tree
  const PrefixPath& prefixPath() const
  {
    return prefixPath_;
  }

protected:
  PreBasis preBasis_;
  PrefixPath prefixPath_;
};


/** \brief The restriction of a finite element basis to a single element */
template<class GB>
class PrismaticProductLocalView
{
public:
  //! The global FE basis that this is a view on
  using GlobalBasis = GB;

  //! The grid view the global FE basis lives on
  using GridView1 = typename GlobalBasis::GridView1;
  using GridView2 = typename GlobalBasis::GridView2;

  //! Type of the grid element we are bound to
  using Element1 = typename GridView1::template Codim<0>::Entity;
  using Element2 = typename GridView2::template Codim<0>::Entity;

  //! The type used for sizes
  using size_type = std::size_t;

  //! Tree of local finite elements / local shape function sets
  using Tree = typename GlobalBasis::PreBasis::Node;
  using Tree1 = typename GlobalBasis::PreBasis1::Node;
  using Tree2 = typename GlobalBasis::PreBasis2::Node;

protected:

  using PreBasis = typename GlobalBasis::PreBasis;

  template <class PB>
  using MultiIndexStorage_t =
      std::conditional_t<(PB::minMultiIndexSize == PB::maxMultiIndexSize),
        OverflowArray<StaticMultiIndex<size_type, PB::maxMultiIndexSize>, PB::multiIndexBufferSize>,
        Dune::ReservedVector<size_type, PB::multiIndexBufferSize>>;

  using MultiIndexStorage = MultiIndexStorage_t<typename GlobalBasis::PreBasis>;
  using MultiIndexStorage1 = MultiIndexStorage_t<typename GlobalBasis::PreBasis1>;
  using MultiIndexStorage2 = MultiIndexStorage_t<typename GlobalBasis::PreBasis2>;

public:

  /** \brief Type used for global numbering of the basis vectors */
  using MultiIndex =
      std::conditional_t<(PreBasis::minMultiIndexSize == PreBasis::maxMultiIndexSize),
        StaticMultiIndex<size_type, PreBasis::maxMultiIndexSize>,
        Dune::ReservedVector<size_type, PreBasis::multiIndexBufferSize>>;


  /** \brief Construct local view for a given global finite element basis */
  PrismaticProductLocalView(const GlobalBasis& globalBasis) :
    globalBasis_(&globalBasis),
    tree_(globalBasis_->preBasis().makeNode())
  {
    // static_assert(models<Concept::BasisTree<GridView>, Tree>(), "Tree type passed to DefaultLocalView does not model the BasisNode concept.");
    initializeTree(tree_);
  }

  void bind1(const Element1& e1)
  {
    element1_ = e1;
    bindTree(tree_.node1_, *element1_);
    indices1_.resize(tree1().size());
    globalBasis_->preBasis1().indices(tree1(), indices1_.begin());
  }

  void bind2(const Element2& e2)
  {
    element2_ = e2;
    bindTree(tree_.node2_, *element2_);
    indices2_.resize(tree2().size());
    globalBasis_->preBasis2().indices(tree2(), indices2_.begin());
  }

  void setup()
  {
    assert(!!element1_);
    assert(!!element2_);
    tree_.setup();
    indices_.resize(size());
    assert(indices_.size() == indices1_.size() * indices2_.size());
    std::size_t stride = indices2_.size();
    std::size_t i = 0;
    for (const auto& idx1 : indices1_) {
      for (const auto& idx2 : indices2_) {
        indices_[i++] = MultiIndexStorage{ idx1[0] * stride + idx2[0] };
      }
    }
  }

  /** \brief Bind the view to a grid element
   *
   * Having to bind the view to an element before being able to actually access any of its data members
   * offers to centralize some expensive setup code in the 'bind' method, which can save a lot of run-time.
   */
  void bind(const Element1& e1, const Element2& e2)
  {
    bind1(e1);
    bind2(e1);
    setup();
  }

  template <class Element>
  void bind(const Element& e)
  {
    bind(std::get<0>(e), std::get<1>(e));
  }

  /** \brief Return if the view is bound to a grid element
   */
  bool bound() const
  {
    return !!element1_ && !!element2_;
  }

  /** \brief Return the grid element that the view is bound to
   *
   * \throws Dune::Exception if the view is not bound to anything
   */
  auto element() const
  {
    return std::tie(*element1_, *element2_);
  }

  /** \brief Unbind from the current element
   *
   * Calling this method should only be a hint that the view can be unbound.
   */
  void unbind()
  {
    element1_.reset();
    element2_.reset();
  }

  /** \brief Return the local ansatz tree associated to the bound entity
   *
   * \returns Tree // This is tree
   */
  const Tree& tree() const
  {
    return tree_;
  }

  const Tree1& tree1() const { return tree_.node1_; }
  const Tree2& tree2() const { return tree_.node2_; }

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
    return indices_[i];
  }

  /** \brief Return the global basis that we are a view on
   */
  const GlobalBasis& globalBasis() const
  {
    return *globalBasis_;
  }

  const PrismaticProductLocalView& rootLocalView() const
  {
    return *this;
  }

protected:
  const GlobalBasis* globalBasis_;
  std::optional<Element1> element1_;
  std::optional<Element2> element2_;
  Tree tree_;

  std::vector<MultiIndexStorage1> indices1_;
  std::vector<MultiIndexStorage2> indices2_;
  std::vector<MultiIndexStorage> indices_;
};


namespace BasisFactory {

  /**
   * \brief A factory that can create a prismatic product of pre-bases
   *
   * \ingroup FunctionSpaceBasesImplementations
   *
   * \tparam PBF1  Type of the first prebasis factory
   * \tparam PBF2  Type of the second prebasis factory
   */
  template<class PBF1, class PBF22>
  auto prismaticProduct(const PBF1& pbf1, const PBF22& pbf2)
  {
    return [&](const auto& gridView1, const auto& gridView2) {
      return PrismaticProductPreBasis{pbf1(gridView1), pbf2(gridView2)};
    };
  }

} // end namespace BasisFactory


/** \brief Basis of a prismatic product of two bases
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV The GridView that the space is defined on
 * \tparam B1  Type of the first finite-element basis
 * \tparam B2  Type of the second finite-element basis
 */
template<class B1, class B2>
using PrismaticProductBasis
  = DefaultGlobalBasis<PrismaticProductPreBasis<typename B1::PreBasis, typename B1::PreBasis> >;

} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PRISMATICPRODUCTBASIS_HH
