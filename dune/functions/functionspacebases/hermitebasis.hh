// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HERMITEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HERMITEBASIS_HH

#include <dune/common/exceptions.hh>
#include <type_traits>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/localfunctions/hermite.hh>

#include <dune/functions/functionspacebases/nodes.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/globalvaluedlocalfiniteelement.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dune/localfunctions/tensormatvec.hh>
#include <dune/functions/functionspacebases/lineartransformedlocalfiniteelement.hh>

namespace Dune {
namespace Functions {
namespace Impl {

/**
 * \brief Implements a transformation from pullbacks of reference basis functions to global
 * basisfunctions for Hermite elements. For more information, see
 * dune-functions/dune/functions/functionspacebases/globalvaluedlocalfiniteelement.hh and
 * c1element/dune/functions/functionspacebases/lineartransformedlocalfinitelement.hh
 * Primary Template uses a sparse Matrix approach and the class
 * LinearTransformedLocalFiniteElement and is valid for all dimensions. Spezialized template
 * for dim == 1 and dim == 2 perform the transformation "by hand" using the class
 * GlobalValuedLocalFiniteElement. They are approx 1.4 times faster than primary template
 * (for average binding + evaluation time of shapefunction and local gradients with usual quad
 * rule).
 * \tparam dim dimension of the underlying grid
 * \tparam useSpezialization whether to
 * use the spezialized versions for dim 1 and 2. Ignored if dim == 3
 */

template <unsigned int dim, bool useSpecialization, class R = double>
struct HermiteTransformator {

    /** \brief ElementInformation object
     *  Stores information that specify the directions of derivative dofs of an element
     *  Also stores, if those directions are modified (i.e. not the coordinate axes) and if they
     * are to be used for dirichlet interpolation (i.e. tangential)
     * // TODO include a check, whether strong enforcement of Dirichlet conditions is possible
     * and throw Error on call to isDirichlet if not
     */
    template <class Element>
    class ElementInformation {
        using GlobalCoordinate = typename Element::Geometry::GlobalCoordinate;
        using ctype = typename Element::Geometry::ctype;

        static_assert(std::is_same<GlobalCoordinate, FieldVector<ctype, dim>>::
                          value); // not sure my code works for other types of global coordinates
      public:
        ElementInformation() {
          for (std::size_t i = 0; i < dim + 1; ++i) {
            // default directions are global coordinates
            derivativeDirections_[i] = 0;
            for (std::size_t j = 0; j < dim; ++j)
              derivativeDirections_[i][j][j] = 1.;
            setTangential_[i] = 0; // Nothing set, we dont known which dof is dirichlet
          }
        }
        ElementInformation(std::array<FieldMatrix<ctype, dim, dim>, dim + 1> directions,
                           std::array<std::bitset<2 * dim>, dim + 1> bitset)
            : derivativeDirections_(directions), setTangential_(bitset) {}

        /**
         * \brief return whether the some direction of vertex"th vertex was modified
         * \param vertex
         * \return bool
         * */
        bool isModified(std::size_t vertex) const {
          bool ret = false;
          for (std::size_t i = 0; i < dim; ++i)
            ret |= setTangential_[vertex][i];
          return ret;
        }
        /**
         * \brief return whether the "direction"th direction of vertex"th vertex is a Dirichlet
         * dof (i.e. is a tangential)
         * \param vertex
         * \param direction
         * \return bool
         */
        bool isDirichlet(std::size_t vertex, std::size_t direction) const {
          if (dim == 1)
            return false;
          return setTangential_[vertex][dim + direction];
        }

        /**
         * \brief return the dof corresponding to localKey is a Dirichlet
         * dof (i.e. is a tangential)
         * \param localKey
         *
         * \return bool
         */
        bool isDirichlet(LocalKey localKey) const {
          if (localKey.codim() == dim) {
            if (localKey.index() == 0) // value -> always Dirichlet
              return true;
            else if (localKey.index() == 1) // first direction
              return isDirichlet(localKey.subEntity(), 0);
            else if (localKey.index() == 2) // second direction
              return isDirichlet(localKey.subEntity(), 1);
          }

          DUNE_THROW(NotImplemented, "Invalid LocalKey");
        }
        /**
         * \brief Get the Derivative Directions object, and array storing the directions of
         * derivative dofs as rows of a fieldmatrix for each vertex
         *
         * \return std::array<FieldMatrix<ctype, dim, dim>, dim + 1> const&
         */
        std::array<FieldMatrix<ctype, dim, dim>, dim + 1> const &getDerivativeDirections() const {
          return derivativeDirections_;
        }

      private:
        std::array<FieldMatrix<ctype, dim, dim>, dim + 1> derivativeDirections_;
        std::array<std::bitset<2 * dim>, dim + 1>
            setTangential_; // first dim bits encode whether direction was modified, second dim
                            // bits encode whether to use this dof for dirichlet information
    };

  private:
    static_assert(dim > 0 && dim < 4);
    static constexpr std::size_t nDofs = (dim == 1) ? 4 : (dim == 2) ? 10 : 20;
    static constexpr std::size_t numberOfVertices = dim + 1;
    static constexpr std::size_t numberOfInnerDofs = (dim - 1) * (dim - 1);
    static constexpr std::size_t numberOfVertexDofs = numberOfVertices * numberOfVertices;

  public:
    HermiteTransformator() : mat_(nDofs, nDofs, BCRSMatrix<R>::random) { setupMatrix(); }

    template <class Element>
    void bind(Element const &e, ElementInformation<Element> const &elementInformation) {
      fillMatrix(e.geometry(), elementInformation);
    }

    template <typename Values>
    void apply(Values &values) const {
      Values tmp = values;
      mat_.mv(tmp, values);
    }

  private:
    void setupMatrix() {
      for (std::size_t i = 0; i < numberOfVertices; ++i) {
        mat_.setrowsize(numberOfVertices * i, 1);
        for (std::size_t j = 1; j <= dim; ++j)
          mat_.setrowsize(numberOfVertices * i + j, dim);
      }
      for (std::size_t i = 0; i < numberOfInnerDofs; ++i)
        mat_.setrowsize(numberOfVertexDofs + i, 1);
      mat_.endrowsizes();

      for (std::size_t i = 0; i < numberOfVertices; ++i) {
        mat_.addindex(numberOfVertices * i, numberOfVertices * i);
        for (std::size_t j = 1; j <= dim; ++j) {
          for (std::size_t k = 1; k <= dim; ++k)
            mat_.addindex(numberOfVertices * i + j, numberOfVertices * i + k);
        }
      }
      for (std::size_t i = 0; i < numberOfInnerDofs; ++i)
        mat_.addindex(numberOfVertexDofs + i, numberOfVertexDofs + i);
      mat_.endindices();
    }

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
    template <class Geometry, typename Element>
    void fillMatrix(Geometry const &geometry, ElementInformation<Element> const &elementInfo) {
      auto const &refElement = Dune::ReferenceElements<typename Geometry::ctype, dim>::simplex();
      auto const &directions = elementInfo.getDerivativeDirections();
      // auto jacobianTransposed = geometry.jacobianTransposed(xi);
      for (std::size_t i = 0; i < numberOfVertices; ++i) // dim + 1 vertices
      {
        Dune::FieldMatrix<R, dim, dim> directionalJacobianTransposed =
            geometry.jacobianTransposed(refElement.position(i, dim));
        if (elementInfo.isModified(i)) {
          FieldMatrix<R, dim, dim> dir = directions[i];
          dir.invert();
          directionalJacobianTransposed = directionalJacobianTransposed * dir;
        }
        mat_[numberOfVertices * i][numberOfVertices * i] = 1.;
        for (std::size_t j = 0; j < dim; ++j)
          for (std::size_t k = 0; k < dim; ++k) {
            mat_[numberOfVertices * i + 1 + j][numberOfVertices * i + 1 + k] =
                directionalJacobianTransposed[k][j]; // jacobianTransposed[k][j];
          }
      }
      for (std::size_t i = 0; i < (dim - 1) * (dim - 1); ++i) // inner dofs
      {
        mat_[numberOfVertexDofs + i][numberOfVertexDofs + i] = 1.;
      }
    }

    BCRSMatrix<R> mat_;

  public:
    /**
     * \brief Class that evaluates the push forwards of the global nodes of a LocalFunction.
     *        It stretches the LocalInterpolation interface, because we evaluate the
     *        derivatives of f
     * \tparam LocalBasis Basis of the reference element, to get the Domain and RangeFieldType
     * \tparam Element
     */
    template <class LocalBasis, class Element>
    class GlobalValuedInterpolation {
        using size_type = std::size_t;
        using LocalCoordinate = typename LocalBasis::Traits::DomainType;
        using RFT = typename LocalBasis::Traits::RangeFieldType;
        using ctype = typename Element::Geometry::ctype;
        static constexpr size_type numberOfVertices = dim + 1;
        static constexpr size_type innerDofCodim = (dim == 2) ? 0 : 1; // probably wrong for dim > 3
        static constexpr size_type numberOfInnerDofs =
            (dim - 1) * (dim - 1); // probably wrong for dim > 3
      public:
        GlobalValuedInterpolation() {}
        void bind(const Element &element, const ElementInformation<Element> &elementInfo) {
          elementInfo_ = elementInfo;
        }

        /** \brief Evaluate a given function at the Lagrange nodes
         *
         * \tparam F Type of function to evaluate
         * \tparam C Type used for the values of the function
         * \param[in] ff Function to evaluate
         * \param[out] out Array of function values
         */
        template <typename F, typename C>
        void interpolate(const F &f, std::vector<C> &out) const {


          auto df = derivative(f);
          out.resize(LocalBasis::size());

          auto const &refElement = Dune::ReferenceElements<ctype, dim>::simplex();
          Dune::FieldVector<RFT, dim> derivativeValue;
          // Iterate over vertices, dim dofs per vertex
          for (size_type i = 0; i < dim + 1; ++i) {
            LocalCoordinate x = refElement.position(i, dim);

            // matrix storing directions as rows (!)
            auto const &directionMatrix = elementInfo_.getDerivativeDirections()[i];
            // auto const &derivativeValue = vectorToMatrix(df(x)) * transpose(directionMatrix);
            directionMatrix.mv(matrixToVector(df(x)), derivativeValue);
            out[i * numberOfVertices] = f(x);
            for (size_type d = 0; d < dim; ++d)
              out[i * numberOfVertices + d + 1] = derivativeValue[d];
          }

          for (size_type i = 0; i < numberOfInnerDofs; ++i) {
            out[numberOfVertices * numberOfVertices + i] = f(refElement.position(i, innerDofCodim));
          }
        }

      protected:
        ElementInformation<Element> elementInfo_;
    };
};

template <class R>
struct HermiteTransformator<1, true, R> {

    template <typename Values, typename LocalCoordinate, typename Geometry>
    static auto apply(Values &values, const LocalCoordinate &xi, const Geometry &geometry) {
      auto jacobianTransposed = geometry.jacobianTransposed(xi);
      // transform the values of basisfunctions that correspond to a "derivative" node
      for (std::size_t i = 0; i < 2; ++i) {
        Dune::FieldVector<R, 1> localGrad = {values[i * 2 + 1]};
        Dune::FieldVector<R, 1> globalGrad(0.);
        jacobianTransposed.mtv(localGrad, globalGrad);
        values[i * 2 + 1] = globalGrad[0];
      }
    }

    template <typename Gradients, typename LocalCoordinate, typename Geometry>
    static auto applyJacobian(Gradients &values, const LocalCoordinate &xi,
                              const Geometry geometry) {
      auto jacobianTransposed = geometry.jacobianTransposed(xi);
      // transform the values of basisfunctions that correspond to a "derivative" node
      for (std::size_t i = 0; i < 2; ++i) {
        Dune::FieldVector<R, 1> globalDerivative(0.);
        jacobianTransposed.mtv(Dune::FieldVector<R, 1>{values[i * 2 + 1][0][0]}, globalDerivative);

        // TODO Assumes Gradients are a container of FieldMatrix<F,1,1>. Make this viable for
        // FieldVector<F,1> as well
        values[i * 2 + 1][0][0] = globalDerivative[0];
      }
    }

    template <class Function, class LocalCoordinate, class Element>
    class LocalValuedFunction {
        Function f_;
        const Element &element_;

      public:
      template<class F>
        LocalValuedFunction( F&& f, const Element &e) : f_(std::forward<F>(f)), element_(e) {}

        auto operator()(const LocalCoordinate &xi) const {
          auto &&f = Dune::Impl::makeFunctionWithCallOperator<LocalCoordinate>(f_);
          return f(xi); // no transformation for functionvalues needed
        }

        friend auto derivative(LocalValuedFunction const &t) {
          // I would like to do this
          // return LocalValuedFunction<typename std::decay<decltype(derivative(t.f_))>::type,
          // LocalCoordinate, Element>(derivative(t.f_), t.element_); But it leads to a
          // DifferentiableFunction which causes a segmentation violation when calling asInterface()
          // ....
          return LocalValuedFunction{derivative(t.f_), t.element_};
        }
    };
};

template <class R>
struct HermiteTransformator<2, true, R> {
    template <typename Values, typename LocalCoordinate, typename Geometry>
    static auto apply(Values &values, const LocalCoordinate &xi, const Geometry &geometry) {
      auto jacobianTransposed = geometry.jacobianTransposed(xi);
      // transform the values of basisfunctions that correspond to a "derivative" node
      for (std::size_t i = 0; i < 3; ++i) {
        Dune::FieldVector<R, 2> localGrad = {values[i * 3 + 1], values[i * 3 + 2]};
        Dune::FieldVector<R, 2> globalGrad(0.);
        jacobianTransposed.mtv(localGrad, globalGrad);
        // TODO make this nice by removing globalGrad and write directly into values
        values[i * 3 + 1] = globalGrad[0];
        values[i * 3 + 2] = globalGrad[1];
      }
    }

    template <typename Gradients, typename LocalCoordinate, typename Geometry>
    static auto applyJacobian(Gradients &values, const LocalCoordinate &xi,
                              const Geometry geometry) {
      auto jacobianTransposed = geometry.jacobianTransposed(xi);
      // transform the values of basisfunctions that correspond to a "derivative" node
      for (std::size_t i = 0; i < 3; ++i) {
        Dune::FieldVector<R, 2> localValued_dx = {values[i * 3 + 1][0][0], values[i * 3 + 2][0][0]},
                                localValued_dy = {values[i * 3 + 1][0][1], values[i * 3 + 2][0][1]};
        auto globalValued_dx = localValued_dx, globalValued_dy = localValued_dy;
        jacobianTransposed.mtv(localValued_dx, globalValued_dx);
        jacobianTransposed.mtv(localValued_dy, globalValued_dy);

        values[i * 3 + 1][0][0] = globalValued_dx[0];
        values[i * 3 + 1][0][1] = globalValued_dy[0];

        values[i * 3 + 2][0][0] = globalValued_dx[1];
        values[i * 3 + 2][0][1] = globalValued_dy[1];
      }
    }

    template <class Function, class LocalCoordinate, class Element>
    class LocalValuedFunction {
        const Function &f_;
        const Element &element_;

      public:
        LocalValuedFunction(const Function &f, const Element &e) : f_(f), element_(e) {}

        auto operator()(const LocalCoordinate &xi) const {
          auto &&f = Dune::Impl::makeFunctionWithCallOperator<LocalCoordinate>(f_);
          return f(xi); // no transformation for functionvalues needed
        }

        friend auto derivative(LocalValuedFunction const &t) {
          return derivative(t.f_); // TODO return a LocalValuedFunction here
        }
    };
};

// primary template, not implemented
template <typename GV, int dim, typename R>
class HermiteElementInformationMap;

// specialization for dim == 2
/**
 * \brief Class that creates a Mapping from Elements of a GridView to their ElementInformation
 *
 * \tparam GV GridView
 * \tparam R Rangetype of the finite element
 */
template <typename GV, typename R>
class HermiteElementInformationMap<GV, 2, R> {
    static constexpr int dim = 2;
    using D = typename GV::ctype;
    using Element = typename GV::template Codim<0>::Entity;
    static_assert(GV::dimension == dim);
    using ElementInformation =
        typename HermiteTransformator<dim, false, R>::template ElementInformation<Element>;
    using GlobalCoordinate = typename Element::Geometry::GlobalCoordinate;
    using IndexType = typename GV::IndexSet::IndexType;

  public:
    HermiteElementInformationMap(GV const &gv, bool useTangentials = true)
        : elementMapper_(gv, mcmgElementLayout()), elementInformation_(),
          useTangentials_(useTangentials) {
      elementInformation_.resize(gv.size(0));
      if (useTangentials_)
        fill(gv);
    }

    template <class Function>
    HermiteElementInformationMap(GV const &gv, Function const &f)
        : elementMapper_(gv, mcmgElementLayout()), elementInformation_(gv.size(0)),
          tangentialMap_(f) {
      fill(gv);
      useTangentials_ = true;
    }

    void update(GV const &gv) {
      if (!useTangentials_)
        return;
      elementInformation_.resize(gv.size(0));
      elementMapper_.update(gv);
      fill(gv);
    }

    template <class Element>
    const auto &find(const Element &element) const {
      if (useTangentials_)
        return elementInformation_[elementMapper_.index(element)];
      else
        return defaultInfo_;
    }

  private:
    void fill(const GV &gv, std::enable_if_t<dim == 2, int> = 0) {

      // compute orientation for all elements
      const auto &indexSet = gv.indexSet();
      // vector with directions and bools to indicate whether this dir was already set
      FieldMatrix<D, dim, dim> eye = 0;
      for (std::size_t i = 0; i < dim; ++i)
        eye[i][i] = 1.;
      std::tuple<FieldMatrix<D, dim, dim>, std::bitset<2 * dim>> defaultInfo;
      std::get<0>(defaultInfo) = eye;
      std::get<1>(defaultInfo) = 0;
      std::vector<std::tuple<FieldMatrix<D, dim, dim>, std::bitset<2 * dim>>> directionPerVertex(
          indexSet.size(dim), defaultInfo);
      // iterate over intersections
      for (const auto &element : elements(gv)) {
        for (const auto &intersection : intersections(gv, element)) {
          // fill vertex Data with linear independent tangentials
          if (intersection.boundary()) {
            auto facet = intersection.inside().template subEntity<1>(intersection.indexInInside());

            auto startIndex = indexSet.subIndex(facet, 0, dim);
            auto endIndex = indexSet.subIndex(facet, 1, dim);

            if (tangentialMap_) {
              setVertexData(directionPerVertex[startIndex],
                            tangentialMap_(facet.geometry().corner(0)));
              setVertexData(directionPerVertex[endIndex],
                            tangentialMap_(facet.geometry().corner(1)));
            } else {
              auto tangential =
                  intersection.geometry().corner(1) - intersection.geometry().corner(0);
              setVertexData(directionPerVertex[startIndex], tangential);
              setVertexData(directionPerVertex[endIndex], tangential);
            }
          }
          // inner dofs, keep global axes
        }
      }

      // iterate over vector and set the remaining direction to normals
      for (auto &[dir, isSet] : directionPerVertex)
        if (isSet[0]) {
          if (isSet[1])
            continue; // both are tangetials, i.e. vertex is corner
          else        //, set second direction normal to first
          {
            dir[1] = Dune::FieldVector<D, dim>{-dir[0][1], dir[0][0]};
            isSet[1] = true;
            isSet[dim + 1] = false;
          }
        }
      // do nothing if non was set

      // iterate over elements and construct elementinfomap
      for (const auto &element : elements(gv)) {
        auto elementIndex = elementMapper_.index(element);
        std::array<Dune::FieldMatrix<D, dim, dim>, dim + 1> directions;
        std::array<std::bitset<2 * dim>, 3> booleans;
        for (std::size_t i = 0; i < element.subEntities(dim); ++i) {
          auto vertexIndex = indexSet.subIndex(element, i, dim);
          auto const &[dir, b] = directionPerVertex[vertexIndex];
          directions[i] = dir;
          booleans[i] = b;
        }

        elementInformation_[elementIndex] = ElementInformation(directions, booleans);
      }
    }

    Dune::MultipleCodimMultipleGeomTypeMapper<GV> elementMapper_;
    std::vector<ElementInformation> elementInformation_;
    std::function<GlobalCoordinate(GlobalCoordinate)> tangentialMap_;
    ElementInformation defaultInfo_;
    bool useTangentials_;

    // helper functions
    bool linearIndependent(Dune::FieldVector<D, 2> a,
                           Dune::FieldVector<D, 2> b) { // write formula explicitly
      return std::abs(a.dot(Dune::FieldVector<D, dim>{-b[1], b[0]})) > 1e-14;
    }

    // add tangential to vertexData, if there is no direction already that is linear dependent
    void
    setVertexData(std::tuple<FieldMatrix<D, dim, dim>, std::bitset<2 * dim>> &directionPerVertex,
                  GlobalCoordinate const &tangential) {
      auto &[dir, isSet] = directionPerVertex;

      if (!isSet[0]) // none is set
      {              // one could check for linear dependency here, but this prevents scaling
        dir[0] = tangential / tangential.two_norm();
        isSet[0] = true;   // bit for modified
        isSet[dim] = true; // bit for Dirichlet
      } else if (linearIndependent(dir[0], tangential)) {
        if (isSet[1]) {
          if (linearIndependent(dir[1], tangential))
            DUNE_THROW(Dune::NotImplemented,
                       "More than two linear independent tangentials at vertex");
          else
            isSet[dim + 1] = true; // linear dependent to second, dont modify, but set dirichlet bit
        } else                     // second not set
        {
          dir[1] = tangential / tangential.two_norm();
          isSet[1] = true;
          isSet[dim + 1] = true;
        }
      } else
        isSet[dim] = true; // linear dependent to first, dont modify, but set dirichlet bit
    }
};

// spezialization for dim == 1, no need to do anything here
template <typename GV, typename R>
class HermiteElementInformationMap<GV, 1, R> {
    using D = typename GV::ctype;
    using Element = typename GV::template Codim<0>::Entity;
    static_assert(GV::dimension == 1);
    using ElementInformation =
        typename HermiteTransformator<1, false, R>::template ElementInformation<Element>;
    using GlobalCoordinate = typename Element::Geometry::GlobalCoordinate;
    using IndexType = typename GV::IndexSet::IndexType;

  public:
    HermiteElementInformationMap(GV const &gv, bool dummy) {
      std::array<Dune::FieldMatrix<D, 1, 1>, 2> dir;
      std::array<std::bitset<2>, 2> b;
      for (std::size_t i = 0; i < 2; ++i) {
        dir[i] = Dune::FieldMatrix<D, 1, 1>({{1.}});
        b[i] = (unsigned short)1;
      }

      defaultElementInfo_ = ElementInformation(dir, b);
    }
    void update(GV const &gv) {}

    const auto &find(const Element &element) const { return defaultElementInfo_; }

  private:
    ElementInformation defaultElementInfo_;
};

// spezialization for dim == 3, no default solution here, since one vertex on boundary has
// arbitrarily many boundary tangentials
template <typename GV, typename R>
class HermiteElementInformationMap<GV, 3, R> {
    using D = typename GV::ctype;
    using Element = typename GV::template Codim<0>::Entity;
    static_assert(GV::dimension == 3);
    using ElementInformation =
        typename HermiteTransformator<3, false, R>::template ElementInformation<Element>;
    using GlobalCoordinate = typename Element::Geometry::GlobalCoordinate;
    using IndexType = typename GV::IndexSet::IndexType;

  public:
    HermiteElementInformationMap(GV const &gv, bool dummy) {
      std::array<FieldMatrix<D, 3, 3>, 4> data;
      std::array<std::bitset<6>, 4> b;
      for (std::size_t vertex = 0; vertex < 4; ++vertex) {

        data[vertex] = 0;
        for (std::size_t i = 0; i < 3; ++i)
          data[vertex][i][i] = 1.;
        b[vertex] = 0;
      }
      defaultElementInfo_ = ElementInformation(data, b);
    }
    void update(GV const &gv) {}

    const auto &find(const Element &element) const { return defaultElementInfo_; }

  private:
    ElementInformation defaultElementInfo_;
};
} // namespace Impl

// primary template, not implemented
template <typename GV, typename R, bool useSpezialization>
class HermiteNode;

// specialization for GlobalValuedLocalFiniteElement
// Does not make use of ElementInformation structure. No additional structure for Dirichlet
// interpolation.
template <typename GV, typename R>
class HermiteNode<GV, R, true> : public LeafBasisNode {
    static constexpr unsigned int dim = GV::dimension;
    using LocalValuedFE = HermiteLocalFiniteElement<typename GV::ctype, R, dim>;

  public:
    using size_type = std::size_t;
    using Element = typename GV::template Codim<0>::Entity;

    using FiniteElement =
        Impl::GlobalValuedLocalFiniteElement<Impl::HermiteTransformator<dim, true>, LocalValuedFE,
                                             Element>;
    HermiteNode()
        : localValuedFiniteElement_(std::make_shared<LocalValuedFE>()),
          finiteElement_(std::make_shared<FiniteElement>()) {
      this->setSize(finiteElement_->size());
    }

    ~HermiteNode() {}

    //! Return current element, throw if unbound
    const Element &element() const { return element_; }

    /** \brief Return the LocalFiniteElement for the element we are bound to
     *
     * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions
     * module
     */
    const FiniteElement &finiteElement() const { return *finiteElement_; }

    //! Bind to element.
    void bind(const Element &e) {
      if (not e.type().isSimplex())
        DUNE_THROW(Dune::NotImplemented, "HermiteBasis can only be bound to simplex elements");
      element_ = e;

      finiteElement_->bind(*localValuedFiniteElement_, element_);
    }

  protected:
    unsigned int order() const { return 3; }
    std::shared_ptr<LocalValuedFE> localValuedFiniteElement_;
    std::shared_ptr<FiniteElement> finiteElement_;
    Element element_;
};

// spezialization for LinearTransformedLocalFiniteElement
template <typename GV, typename R>
class HermiteNode<GV, R, false> : public LeafBasisNode {
    static constexpr unsigned int dim = GV::dimension;
    using LocalValuedFE = HermiteLocalFiniteElement<typename GV::ctype, R, dim>;
    using ElementInformationMap = Impl::HermiteElementInformationMap<GV, dim, R>;

  public:
    using size_type = std::size_t;
    using Element = typename GV::template Codim<0>::Entity;

    using FiniteElement =
        Impl::LinearTransformedLocalFiniteElement<Impl::HermiteTransformator<dim, false, R>,
                                                  LocalValuedFE, Element>;
    HermiteNode(const ElementInformationMap *elementInformationMap)
        : localValuedFiniteElement_(std::make_shared<LocalValuedFE>()),
          finiteElement_(std::make_shared<FiniteElement>()),
          elementInformationMap_(elementInformationMap) {
      this->setSize(finiteElement_->size());
    }

    ~HermiteNode() {}

    //! Return current element, throw if unbound
    const Element &element() const { return element_; }

    /** \brief Return the LocalFiniteElement for the element we are bound to
     *
     * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions
     * module
     */
    const FiniteElement &finiteElement() const { return *finiteElement_; }

    //! Bind to element.
    void bind(const Element &e) {
      if (not e.type().isSimplex())
        DUNE_THROW(Dune::NotImplemented, "HermiteBasis can only be bound to simplex elements");
      element_ = e;

      finiteElement_->bind(*localValuedFiniteElement_, element_,
                           elementInformationMap_->find(element_));
      this->setSize(finiteElement_->size());
    }

    unsigned int order() const { return 3; }

  protected:
    std::shared_ptr<LocalValuedFE> localValuedFiniteElement_;
    std::shared_ptr<FiniteElement> finiteElement_;
    const ElementInformationMap *elementInformationMap_;
    Element element_;
};

/**
 * \brief A pre-basis for a Hermitebasis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam GV  The grid view that the FE basis is defined on
 * \tparam R   Range type used for shape function values
 * \note This only works for simplex grids
 */
template <typename GV, typename R, bool useSpecialization>
class HermitePreBasis {
  public:
    //! The grid view that the FE basis is defined on
    using GridView = GV;

    //! Type used for indices and size information
    using size_type = std::size_t;

    //! Template mapping root tree path to type of created tree node
    using Node = HermiteNode<GridView, R, useSpecialization>;

    static constexpr size_type maxMultiIndexSize = 1;
    static constexpr size_type minMultiIndexSize = 1;
    static constexpr size_type multiIndexBufferSize = 1;

  private:
    static const size_type dim = GV::dimension;
    static constexpr size_type innerDofCodim = (dim == 2) ? 0 : 1;
    using ElementInformationMap = Impl::HermiteElementInformationMap<GV, dim, R>;

  public:
    //! Constructor for a given grid view object
    HermitePreBasis(const GV &gv, bool useTangential)
        : gridView_(gv), elementInformationMap_(gv, useTangential), useDirections_(useTangential) {
      if (dim > 3)
        DUNE_THROW(Dune::NotImplemented, "HermitePreBasis only implemented for dim <= 3");
    }

    //! Constructor for a given grid view object and function specifying direction at boundary
    //! vertices
    template <class Function>
    HermitePreBasis(const GV &gv, Function const &f)
        : gridView_(gv), elementInformationMap_(gv, f), useDirections_(true) {
      if (dim > 3)
        DUNE_THROW(Dune::NotImplemented, "HermitePreBasis only implemented for dim <= 3");
    }

    //! Initialize the global indices
    void initializeIndices() {}

    //! Obtain the grid view that the basis is defined on
    const GridView &gridView() const { return gridView_; }

    //! Update the stored grid view, to be called if the grid has changed
    void update(const GridView &gv) {
      gridView_ = gv;
      elementInformationMap_.update(gv);
    }

    /**
     * \brief Create tree node
     */
    Node makeNode() const {
      if constexpr (useSpecialization)
        return Node{};
      else
        return Node(&elementInformationMap_);
    }

    //! Same as size(prefix) with empty prefix
    size_type size() const {
      if (dim != 0 && dim <= 3)
        return (dim + 1) * gridView_.size(dim) + ((dim == 1) ? 0 : gridView_.size(innerDofCodim));
      else
        DUNE_THROW(Dune::NotImplemented, "No size method for " << dim << "d grids available yet!");
    }

    //! Return number of possible values for next position in multi index
    template <class SizePrefix>
    size_type size(const SizePrefix prefix) const {
      // this basically means this is a leaf node with dimrange 1 right?
      assert(prefix.size() == 0 || prefix.size() == 1);
      return (prefix.size() == 0) ? size() : 0;
    }

    //! Get the total dimension of the space spanned by this basis
    size_type dimension() const { return size(); }

    //! Get the maximal number of DOFs associated to node for any element
    size_type maxNodeSize() const { return dim == 1 ? 4 : dim == 2 ? 10 : 20; }

    template <typename It>
    It indices(const Node &node, It it) const {
      const auto &gridIndexSet = gridView().indexSet();
      const auto &element = node.element();

      // throw if Element is not simplex
      if (not(element.type().isSimplex()))
        DUNE_THROW(Dune::NotImplemented, "Hermite Basis only implemented for simplex elements");
      for (size_type i = 0, end = node.finiteElement().size(); i < end; ++it, ++i) {
        Dune::LocalKey localKey = node.finiteElement().localCoefficients().localKey(i);
        // TODO probably unnecessary
        if (!gridIndexSet.contains(element))
          DUNE_THROW(Dune::RangeError, "Element is not in gridIndexSet!");

        if constexpr (dim <= 3) {
          *it = {{(size_type)((localKey.codim() == dim)
                                  ? ((dim + 1) *
                                         gridIndexSet.subIndex(element, localKey.subEntity(), dim) +
                                     localKey.index())
                                  : (dim + 1) * gridIndexSet.size(dim) +
                                        gridIndexSet.subIndex(element, localKey.subEntity(),
                                                              innerDofCodim) +
                                        localKey.index())}};
        } else
          DUNE_THROW(Dune::NotImplemented,
                     "indices() for Hermite Elements not implemented for dim > 3");
      }
      return it;
    }

  protected:
    GridView gridView_;
    ElementInformationMap elementInformationMap_;
    bool useDirections_ = true;

}; // class HermitePreBasis

namespace BasisFactory {

/**
 * @brief construct a PreBasisFactory
 *
 * @tparam R RangeFieldType
 * @tparam useSpecialization whether to use the GlobalValuedLocalFiniteElement implementation
 * (faster) or LinearTransformedLocalFiniteElement implementation (possilibily to strongly enforce
 * BC)
 * @param useTangentials whether to have directions of derivative DOFs on the boundary in tangential
 * and normal direction. Necessary for strong enforcement of BC. Ignored if \p useSpecialization is
 * true
 * @return the PreBasisFactory
 */
template <typename R = double, bool useSpecialization = false>
auto hermite(bool useTangentials = true) {
  return [=](const auto &gridView) {
    return HermitePreBasis<std::decay_t<decltype(gridView)>, R, useSpecialization>(gridView,
                                                                                   useTangentials);
  };
}

/**
 * @brief construct a PreBasisFactory
 *
 * @tparam R RangeFieldType
 * @tparam F Function mapping from GlobalCoordinates to GlobalCoordinates
 * @param f Indicates for each point where to point the first derivative to. Second derivative is
 * normal to first.
 * @return the PreBasisFactory
 */

template <class F, typename R = double>
auto hermite(F const &f) {
  return [=](const auto &gridView) {
    return HermitePreBasis<std::decay_t<decltype(gridView)>, R, false>(gridView, f);
  };
}

} // namespace BasisFactory

} // namespace Functions
} // namespace Dune

#endif
