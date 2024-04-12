// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HERMITEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HERMITEBASIS_HH

#include <algorithm>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <type_traits>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/genericlocalfiniteelements/hermite.hh>
#include <dune/functions/functionspacebases/globalvaluedlocalfiniteelement.hh>
#include <dune/functions/functionspacebases/leafprebasismappermixin.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/lineartransformedlocalfiniteelement.hh>
#include <dune/istl/bcrsmatrix.hh>


namespace Dune
{
namespace Functions
{
namespace Impl
{


// Helper function returning an unordered range
// of global indices associated to the element.
// This could be implemented cheaper internally in
// the MCMGMapper by storing a precomputed
// container of all subsentities addressed by the layout.
template<class GridView>
auto subIndexSet(Dune::MultipleCodimMultipleGeomTypeMapper<GridView> const &mapper,
                 const typename GridView::template Codim<0>::Entity &element)
{
  using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
  using Index = typename Mapper::Index;
  constexpr auto dimension = GridView::dimension;
  auto subIndices = std::vector<Index>();
  auto referenceElement = Dune::referenceElement<double, dimension>(element.type());
  for (auto codim : Dune::range(dimension + 1)) {
    for (auto subEntity : Dune::range(referenceElement.size(codim))) {
      std::size_t c = mapper.layout()(referenceElement.type(subEntity, codim), dimension);
      if (c > 0) {
        std::size_t firstIndex = mapper.subIndex(element, subEntity, codim);
        for (auto j : Dune::range(firstIndex, firstIndex + c)) {
          subIndices.push_back(j);
        }
      }
    }
  }
  return subIndices;
}

  // Helper function computing an average mesh size per subentity
  // by averaging over the adjacent elements. This only considers
  // the subentities handled by the given mapper and returns a
  // vector of mesh sizes indixed according to the mapper.
  template<class FieldType = double, class Mapper>
  auto computeAverageSubEntityMeshSize(const Mapper& mapper)
  {
    constexpr auto dimension = Mapper::GridView::dimension;

    std::vector<unsigned int> adjacentElements(mapper.size(), 0);
    std::vector<FieldType> subEntityMeshSize(mapper.size(), 0.0);
    for(const auto& element : Dune::elements(mapper.gridView()))
    {
      auto A = element.geometry().volume();
      for(auto i : Impl::subIndexSet(mapper, element))
      {
        subEntityMeshSize[i] += A;
        ++(adjacentElements[i]);
      }
    }
    for(auto i : Dune::range(mapper.size()))
      subEntityMeshSize[i] = std::pow(subEntityMeshSize[i]/adjacentElements[i], 1./dimension);
    return subEntityMeshSize;
  }

/** \brief Helper struct describing the Traits of the Global State of an Hermite Element*/
template<class M, class R>
struct HermiteGlobalStateTraits {
    using Mapper = M;
    using GlobalState = std::tuple<Mapper, std::vector<R>>;
    using LocalState = std::vector<R>;
};


/** \brief This class implements the Transformation for Hermite finite elements, that 'corrects' the non affine-equivalence.
*   It is bindable to an Element, stateful and offers a transform(...) method.
*   Its state can change upon binding and can be accessed via localState().
*   Additionally, it offers a GlobalValuedInterpolation class to be used for interpolation.
*   \tparam Element   The Grid Element
*   \tparam R         The Fieldtype of the Finite Element
*   \tparam reduced   If True, use the reduced Hermite element aka Discrete Kirchhoff Triangle
*/
template<class Element, class R, class GlobalStateTraits, bool reduced = false>
struct HermiteTransformator {

  private:
    static constexpr int dim = Element::mydimension;
    static_assert(dim > 0 && dim < 4);
    static_assert(!(reduced && (dim != 2))); // TODO is there a reduced 3d version ?
    static constexpr std::size_t nDofs = (dim == 1) ? 4 : (dim == 2) ? 10 : 20;
    static constexpr std::size_t numberOfVertices = dim + 1;
    static constexpr std::size_t numberOfInnerDofs = reduced ? 0 : (dim - 1) * (dim - 1);
    static constexpr std::size_t numberOfVertexDofs = numberOfVertices * numberOfVertices;
  public:
    using GlobalState = typename GlobalStateTraits::GlobalState;
    using LocalState = typename GlobalStateTraits::LocalState;
    using size_type = std::size_t;

    HermiteTransformator(GlobalState const &globalState) : globalState_(&globalState) {}

    /** Binds the Transformator to an element.
     */
    void bind(Element const &e)
    {
      localState_ = std::apply(
          [&e](auto &&mapper, auto &&data) {
            LocalState localState;
            for (auto const &index : range(e.subEntities(dim)))
              localState.push_back(data[mapper.subIndex(e, index, dim)]);
            return localState;
          },
          *globalState_);

      fillMatrix(e.geometry(), localState_);
    }

    LocalState const &localState() const { return localState_; }

    //! The size of the transformed finite element.
    static constexpr size_type size()
    {
      if constexpr (dim == 1)
        return 4;
      else if constexpr (dim == 2) {
        if constexpr (reduced)
          return 9;
        else
          return 10;
      } else // dim == 3
        return 20;
    }

    /** Applies the transformation. Note that we do not distinguish for
    Scalar/Vector/Matrix Type,
    * but only assume the Values to be Elements of a Vectorspace.
      We assume random access containers. */
    template<typename InputValues, typename OutputValues, class LocalCoordinate>
    void transform(InputValues const &inValues, OutputValues &outValues, LocalCoordinate const& x) const
    {
      assert(inValues.size() == numberOfVertexDofs + numberOfInnerDofs);
      assert(reduced || (outValues.size() == inValues.size()));
      auto inIt = inValues.begin();
      auto outIt = outValues.begin();

      for (auto vertex : Dune::range(dim + 1)) {
        *outIt = *inIt; // value dof is not transformed
        outIt++, inIt++;
        // transform the gradient dofs together
        for (auto &&[row_i, i] : sparseRange(subMatrices_[vertex])) {
          outIt[i] = 0.;
          for (auto &&[val_i_j, j] : sparseRange(row_i))
            outIt[i] += val_i_j * inIt[j];
        }
        // increase pointer by size of gradient = dim
        outIt += dim, inIt += dim;
      }

      if constexpr (not reduced)
        // copy all remaining inner dofs
        std::copy(inIt, inValues.end(), outIt);
    }

  private:
    /**
     * \brief Create the transformationmatrix m
     *
     * \tparam Geometry
     * \param geometry
     * \return BCRSMatrix with size nDofs x nDofs
     *
     *      |1          0|} repeat
     *      |   J        |} dim + 1 times
     * m =  |     1      |
     *      |       J ...|
     *      |0          1|  } repeat (dim-1)^2 times
     *
     */
    template<class Geometry>
    void fillMatrix(Geometry const &geometry, LocalState const &averageSubEntityMeshSize)
    {
      auto const &refElement = Dune::ReferenceElements<typename Geometry::ctype, dim>::simplex();
      for (std::size_t i = 0; i < numberOfVertices; ++i) // dim + 1 vertices
      {
        subMatrices_[i] = geometry.jacobian(refElement.position(i, dim));
        subMatrices_[i] /= averageSubEntityMeshSize[i];
      }
    }

    // one transformation per vertex
    std::array<Dune::FieldMatrix<R, dim, dim>, numberOfVertices> subMatrices_;
    GlobalState const *globalState_;
    LocalState localState_;

  public:
    /**
     * \brief Class that evaluates the push forwards of the global nodes of a
     * LocalFunction. It stretches the LocalInterpolation interface, because we
     * evaluate the derivatives of f.
     *
     */
    class GlobalValuedInterpolation
    {

        using size_type = std::size_t;
        using ctype = typename Element::Geometry::ctype;
        static constexpr size_type numberOfVertices = dim + 1;
        static constexpr size_type innerDofCodim = (dim == 2) ? 0 : 1; // probably wrong for dim > 3
        static constexpr size_type numberOfInnerDofs =
            (dim - 1) * (dim - 1); // probably wrong for dim > 3

      public:
        GlobalValuedInterpolation(HermiteTransformator const &t) : transformator_(&t) {}

        /** \brief bind the Interpolation to an element and a localInterpolation.*/
        template<class LocalValuedLocalInterpolation>
        void bind([[maybe_unused]] Element const &element,
                  [[maybe_unused]] LocalValuedLocalInterpolation const &localInterpolation)
        {
          localState_ =&( transformator_->localState());
        }
      public:
        /** \brief Evaluate a given function and its derivatives at the nodes
         *
         * \tparam F Type of function to evaluate
         * \tparam C Type used for the values of the function
         * \param[in] f Function to evaluate
         * \param[out] out Array of function values
         */
        template<typename F, typename C>
        void interpolate(const F &f, std::vector<C> &out) const
        {
          auto df = derivative(f);
          out.resize(transformator_->size());

          auto const &refElement = Dune::ReferenceElements<ctype, dim>::simplex();
          // Iterate over vertices, dim dofs per vertex
          for (int i = 0; i < dim + 1; ++i) {
            auto x = refElement.position(i, dim);

            auto derivativeValue = df(x);
            out[i * numberOfVertices] = f(x);
            for (int d = 0; d < dim; ++d)
              out[i * numberOfVertices + d + 1] = getPartialDerivative(derivativeValue,d) * (*localState_)[i];
          }

          if constexpr (not reduced)
            for (size_type i = 0; i < numberOfInnerDofs; ++i) {
              out[numberOfVertices * numberOfVertices + i] =
                  f(refElement.position(i, innerDofCodim));
            }
        }

      protected:
        template<class DerivativeType, class FieldType>
        FieldType getPartialDerivative(DerivativeType const &df, std::size_t i) const
        {
          DUNE_THROW(Dune::NotImplemented, "Derivative Type is neither FieldMatrix<double,1,d> nor "
                                          "FieldVector<double,d>");
        }

        template<class FieldType, int d>
        FieldType getPartialDerivative(Dune::FieldVector<FieldType, d> const &df, std::size_t i) const
        {
          return df[i];
        }

        template<class FieldType, int d>
        FieldType getPartialDerivative(Dune::FieldMatrix<FieldType, 1, d> const &df,
                                      std::size_t i) const
        {
          if (df.N() == 1)
            return df[0][i];
          else if (df.M() == 1)
            return df[i][0];
          else
            DUNE_THROW(Dune::NotImplemented, "Derivative of scalar function is a matrix!");
        }
        HermiteTransformator const *transformator_;
        LocalState const* localState_;
    };

};

} // namespace Impl

/**
 * \brief A pre-basis for a Hermitebasis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam R   Range type used for shape function values
 * \note This only works for simplex grids
 */
template<typename GV, typename R, bool reduced = false>
class HermitePreBasis : public LeafPreBasisMapperMixin<GV>
{
    using Base = LeafPreBasisMapperMixin<GV>;
    using SubEntityMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;

    static constexpr auto cubicHermiteMapperLayout(Dune::GeometryType type, int gridDim)
    {
      if (type.isVertex())
        return 1 + gridDim;
      if (gridDim == 1)
        return 0;
      if ((type.isTriangle()) and (not reduced))
        return 1;
      else
        return 0;
    }

  public:
    //! The grid view that the FE basis is defined on
    using GridView = GV;

    //! Type used for indices and size information
    using size_type = std::size_t;

  private:
    static const size_type dim = GV::dimension;
    static constexpr size_type innerDofCodim = (dim == 2) ? 0 : 1;
    using Element = typename GridView::template Codim<0>::Entity;

    // the following typedefs configure the whole transformation framework
    //! Type used for the generic LocalFiniteElement
    using LocalFE = HermiteLocalFiniteElement<typename GridView::ctype, R, dim, reduced>;
    //! The Traits for the global state of the Hermite elements
    using GlobalStateTraits =
        Impl::HermiteGlobalStateTraits<Dune::MultipleCodimMultipleGeomTypeMapper<GridView>, R>;
    using GlobalState = typename GlobalStateTraits::GlobalState;
    //! Type used for the transformator which turns the generic LocalFiniteElement
    //! into a LocalFiniteElement that is affine equivalent to the global
    //! FiniteElement
    using HermiteTrafo = Impl::HermiteTransformator<Element, R, GlobalStateTraits, reduced>;

  public:
    //! Template mapping root tree path to type of created tree node
    using Node = Impl::TransformedNode<GridView, HermiteTrafo, LocalFE>;

    static constexpr size_type maxMultiIndexSize = 1;
    static constexpr size_type minMultiIndexSize = 1;
    static constexpr size_type multiIndexBufferSize = 1;

  public:
    //! Constructor for a given grid view object
    HermitePreBasis(const GV &gv)
        : Base(gv, cubicHermiteMapperLayout), globalState_({gv, mcmgVertexLayout()}, {})
    {
      updateState(gv);
      if (dim > 3)
        DUNE_THROW(Dune::NotImplemented, "HermitePreBasis only implemented for dim <= 3");
    }

    //! Update the stored grid view, to be called if the grid has changed
    void update(GridView const &gv)
    {
      Base::update(gv);
      updateState(gv);
    }

    /**
     * \brief Create tree node
     */
    Node makeNode() const { return Node(globalState_); }

  protected:
    void updateState(GridView const &gridView)
    {
      auto &[mapper, data] = globalState_;

      mapper.update(gridView);
      data = Impl::computeAverageSubEntityMeshSize<R>(mapper);
    }

    GlobalState globalState_;

}; // class HermitePreBasis

namespace BasisFactory
{

/**
 * \brief construct a PreBasisFactory
 *
 * \tparam R RangeFieldType
 * \return the PreBasisFactory
 */
template<typename R = double>
auto hermite()
{
  return [=](auto const &gridView) {
    return HermitePreBasis<std::decay_t<decltype(gridView)>, R>(gridView);
  };
}

/**
 * \brief construct a PreBasisFactory
 *
 * \tparam R RangeFieldType
 * \return the PreBasisFactory
 */
template<typename R = double>
auto reducedHermite()
{
  return [=](auto const &gridView) {
    return HermitePreBasis<std::decay_t<decltype(gridView)>, R, true>(gridView);
  };
}

} // namespace BasisFactory

} // namespace Functions
} // namespace Dune

#endif
