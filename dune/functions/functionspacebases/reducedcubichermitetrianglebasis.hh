// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REDUCEDCUBICHERMITETRIANGLEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REDUCEDCUBICHERMITETRIANGLEBASIS_HH

/** \file
 * \brief A reduced cubic Hermite triangle global function space basis
 *
 * Based on a standard cubic Lagrange basis on each triangle
 * the basis is constructed as the dual basis to the set of
 * degrees of freedom consisting of: function evaluations in
 * each vertex and evaluations of the partial derivatives in
 * each vertex of the triangle.
 *
 * The (transposed) basis transformation matrix C_ is computed
 * in the 'bind'-method of the LocalFiniteElement class.
 *
 * This function space is involved in the definition of
 * discrete Kirchhoff Triangles (DKT - Elements) for the
 * discretization of Kirchhoff plate bending problems.
 * See "Braess - Finite Elements 5.th Edition p.335" for more details.
 */

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/istl/matrix.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

namespace Dune::Functions
{

    template <typename GV, typename R>
    class ReducedCubicHermiteTriangleLocalFiniteElement;

    template <typename GV, typename R = double>
    class ReducedCubicHermiteTrianglePreBasis;

    /** \brief LocalBasis class in the sense of dune-localfunctions, presenting the restriction
     * of a ReducedCubicHermiteTriangleBasis to an element
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     * \tparam GV Grid view that the basis is defined on
     * \tparam R Number type used for function values
     */
    template <class GV, class R>
    class ReducedCubicHermiteTriangleLocalBasis
    {
        friend class ReducedCubicHermiteTriangleLocalFiniteElement<GV, R>;

        typedef typename GV::ctype D;
        enum
        {
            dim = GV::dimension
        };

    public:
        using Element = typename GV::template Codim<0>::Entity;
        using grid_view = GV;
        using Range = R;

        //! \brief export type traits for function signature
        typedef LocalBasisTraits<D, dim, FieldVector<D, dim>, R, 1, FieldVector<R, 1>, FieldMatrix<R, 1, dim>> Traits;

        //! \brief Constructor with a given PreBasis and LocalFiniteElement
        ReducedCubicHermiteTriangleLocalBasis(const ReducedCubicHermiteTriangleLocalFiniteElement<GV, R> &lFE)
            : lFE_(&lFE)
        {}

        /** \brief  Evaluate all shape functions
         *
         * There is a total of 9 shape functions and the indices for the Output vector are as follows:
         * Indices 0,1,2 return the values of the shape function corresponding to the function evaluation at the vertices 0,1,2, respectively.
         * Indices 3,4,5 return the values of the shape function corresponding to the partial derivative w.r.t x at the vertices 0,1,2, respectively.
         * Indices 6,7,8 return the values of the shape function corresponding to the partial derivative w.r.t y at the vertices 0,1,2, respectively.
         *
         *  \param in Coordinates where to evaluate the functions, in local coordinates
         *  \param out Output vector of values of shape function values
         */
        void evaluateFunction(const typename Traits::DomainType &in,
                              std::vector<typename Traits::RangeType> &out) const
        {
            auto element = lFE_->element();
            auto geometry = element.geometry();

            std::vector<FieldVector<double, 1>> values;
            P3LocalFiniteElement_.localBasis().evaluateFunction(in, values);

            FieldVector<R, 10> tmp(0);

            for(std::size_t j=0; j<P3LocalFiniteElement_.size(); j++)
                tmp[j] = values[j][0];

            out.resize(size());

            // multiply Lagrange Basis-evaluation with basis transformation matrix.
            (*lFE_).C_.mtv(tmp, out);
        }

        /** \brief Evaluate Jacobian of all shape functions
         *
         * There is a total of 9 shape functions and the indices for the Output vector are as follows:
         * Indices 0,1,2 return the Jacobian of the shape function corresponding to the function evaluation at the vertices 0,1,2, respectively.
         * Indices 3,4,5 return the Jacobian of the shape function corresponding to the partial derivative w.r.t x at the vertices 0,1,2, respectively.
         * Indices 6,7,8 return the Jacobian of the shape function corresponding to the partial derivative w.r.t y at the vertices 0,1,2, respectively.
         *
         *  \param in Coordinates where to evaluate the Jacobian, in local coordinates
         *  \param out The Jacobians of all shape functions at the point x
         */
        void evaluateJacobian(const typename Traits::DomainType &in,
                              std::vector<typename Traits::JacobianType> &out) const
        {

            auto element = lFE_->element();
            auto geometry = element.geometry();

            std::vector<FieldMatrix<double, 1, dim>> referenceGradients;
            P3LocalFiniteElement_.localBasis().evaluateJacobian(in, referenceGradients);

            const auto jacobian = element.geometry().jacobianInverseTransposed(in);

            std::vector<FieldVector<R, dim>> gradients(referenceGradients.size());
            for (size_t i = 0; i<gradients.size(); i++)
                jacobian.mv(referenceGradients[i][0], gradients[i]);

            auto transformationMatrix = (*lFE_).C_.transposed();

            // Multiply Jacobian evaluation with basis transformation matrix
            out.resize(size());
            for(std::size_t i=0; i<transformationMatrix.N(); i++)
            for(std::size_t j=0; j<transformationMatrix.M(); j++)
            {
                out[i][0][0] += transformationMatrix[i][j]*gradients[j][0];
                out[i][0][1] += transformationMatrix[i][j]*gradients[j][1];
            }
        }

        //! \brief Polynomial order of the shape functions
        unsigned int order() const
        {
            return 3;
        }

        //! \brief Return the number of basis functions
        std::size_t size() const
        {
            return lFE_->size();
        }

        //! Return current element, throw if unbound
        const Element &element() const
        {
            return *(*lFE_).element_;
        }

        //! Update local finite element
        void updateLFE(ReducedCubicHermiteTriangleLocalFiniteElement<GV, R> &localFiniteElement)
        {
            lFE_ = &localFiniteElement;
        }

    private:
        const ReducedCubicHermiteTriangleLocalFiniteElement<GV, R> *lFE_;

        Dune::LagrangeSimplexLocalFiniteElement<R, double, dim, 3> P3LocalFiniteElement_;
    };

    //! \brief Associations of degrees of freedom to subentities
    class ReducedCubicHermiteTriangleLocalCoefficients
    {
    public:
        //! number of coefficients
        static constexpr std::size_t size()
        {
            return 9;
        }

        //! get i'th index
        const LocalKey localKey(std::size_t i) const
        {
            return LocalKey(i % 3, 2, static_cast<int>(std::floor(i / 3)));
        }
    };

    /** \brief Evaluate the degrees of freedom of a reduced cubic Hermite triangle basis
     * \tparam LocalBasis The corresponding set of shape functions
     */
    template <class LocalBasis>
    class ReducedCubicHermiteTriangleLocalInterpolation
    {
    public:
        using GV = typename LocalBasis::grid_view;
        using R = typename LocalBasis::Range;
        using Element = typename LocalBasis::Element;

        ReducedCubicHermiteTriangleLocalInterpolation(const ReducedCubicHermiteTriangleLocalFiniteElement<GV, R> &lFE)
            : lFE_(&lFE)
        {}

        /** \brief Local interpolation of a function
         *
         * The evaluation of f happens in world coordinates!
         */
        template <typename F, typename C>
        void interpolate(const F &f, std::vector<C> &out) const
        {
            auto&& geometry = lFE_->element().geometry();

            if (geometry.type() != GeometryTypes::triangle)
              DUNE_THROW(NotImplemented, "ReducedCubicHermiteTriangleLocalInterpolation is only implemented for triangle elements!");

            out.resize(9);

            auto df = derivative(f);

            // Loop over the vertices
            for (std::size_t i=0; i<3; i++)
            {
                // Value-type degrees of freedom
                out[i+0] = f(geometry.corner(i));

                // Partial-derivative-type degrees of freedom
                auto dfv = df(geometry.corner(i));
                out[i+3] = dfv[0];
                out[i+6] = dfv[1];
            }
        }

        //! Update local finite element
        void updateLFE(ReducedCubicHermiteTriangleLocalFiniteElement<GV, R> &localFiniteElement)
        {
            lFE_ = &localFiniteElement;
        }

    private:
        const ReducedCubicHermiteTriangleLocalFiniteElement<GV, R> *lFE_;
    };

    /** \brief LocalFiniteElement in the sense of dune-localfunctions, for the reduced cubic Hermite triangle basis
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     * This class ties together the implementation classes
     * ReducedCubicHermiteTriangleLocalBasis, ReducedCubicHermiteTriangleLocalCoefficients,
     * and ReducedCubicHermiteTriangleLocalInterpolation
     *
     * \tparam D Number type used for domain coordinates
     * \tparam R Number type used for spline function values
     */
    template <class GV, class R>
    class ReducedCubicHermiteTriangleLocalFiniteElement
    {
        typedef typename GV::ctype D;
        enum
        {
            dim = GV::dimension
        };

        friend class ReducedCubicHermiteTriangleLocalBasis<GV, R>;

    public:
        using Element = typename GV::template Codim<0>::Entity;

        /** \brief Export various types related to this LocalFiniteElement
         */
        typedef LocalFiniteElementTraits<ReducedCubicHermiteTriangleLocalBasis<GV, R>,
                                         ReducedCubicHermiteTriangleLocalCoefficients,
                                         ReducedCubicHermiteTriangleLocalInterpolation<ReducedCubicHermiteTriangleLocalBasis<GV, R>>>
            Traits;

        /** \brief Constructor with a given ReducedCubicHermiteTriangle basis
         */
        ReducedCubicHermiteTriangleLocalFiniteElement()
            : localBasis_(*this),
              localInterpolation_(*this),
              element_(nullptr)
              {}

        /** \brief Bind LocalFiniteElement to a specific element
         *
         *  Compute the basis transformation matrix on each element for the
         *  transformation of the standard cubic Lagrange basis to the dual basis
         *  corresponding to the degrees of freedom: function values and gradients at the nodes.
         */
        void bind(const Element &e)
        {
            element_ = &e;
            auto geometry = element_->geometry();

            /*
             * Left-hand side matrix N are the degrees of freedom
             * applied to the cubic Lagrange basis.
             */
            FieldMatrix<double, 10, 10> N(0);

            for (int i = 0; i < geometry.corners(); i++)
            {
                auto c = geometry.corner(i);
                auto center = geometry.center();
                auto diff = center - c;
                std::vector<FieldVector<double, 1>> values;
                P3LocalFiniteElement_.localBasis().evaluateFunction(geometry.local(c), values);

                std::vector<FieldVector<double, 1>> valuesCenter;
                P3LocalFiniteElement_.localBasis().evaluateFunction(geometry.local(center), valuesCenter);

                const auto jacobian = e.geometry().jacobianInverseTransposed(geometry.local(c));

                std::vector<FieldMatrix<double, 1, dim>> referenceGradients;
                P3LocalFiniteElement_.localBasis().evaluateJacobian(geometry.local(c), referenceGradients);

                std::vector<FieldVector<R, dim>> gradients(referenceGradients.size());
                for (size_t i = 0; i<gradients.size(); i++)
                   jacobian.mv(referenceGradients[i][0], gradients[i]);

                for(std::size_t j=0; j<P3LocalFiniteElement_.size(); j++)
                {
                    // first three degrees of freedom: point-evaluations at corners
                    N[i][j] = values[j][0];
                    // next three degrees of freedom: point-evaluations of partial_x derivative at corners
                    N[i+3][j] = gradients[j][0];
                    // next three degrees of freedom: point-evaluations of partial_y derivative at corners
                    N[i+6][j] = gradients[j][1];
                    //  The following "kinematic condition" removes degree of freedom in the center of the element.
                    N[9][j] += valuesCenter[j][0] - (values[j][0] + gradients[j]*diff);
                }
            }

            // Right-Hand side is identity matrix (of size 9x9) with added zero-row
            FieldMatrix<double, 10, 9> b(0);
            for (std::size_t i = 0; i < size(); i++)
                b[i][i] = 1.0;

            N.invert(); // Todo: improve

            C_ = b.leftmultiply(N);

            // Update local finite element
            localBasis_.updateLFE(*this);
            localInterpolation_.updateLFE(*this);
        }

        /** \brief Hand out a LocalBasis object */
        const ReducedCubicHermiteTriangleLocalBasis<GV, R> &localBasis() const
        {
            return localBasis_;
        }

        /** \brief Hand out a LocalCoefficients object */
        const ReducedCubicHermiteTriangleLocalCoefficients &localCoefficients() const
        {
            return localCoefficients_;
        }

        /** \brief Hand out a LocalInterpolation object */
        const ReducedCubicHermiteTriangleLocalInterpolation<ReducedCubicHermiteTriangleLocalBasis<GV, R>> &localInterpolation() const
        {
            return localInterpolation_;
        }

        // //! \brief Number of shape functions in this finite element
        std::size_t size() const
        {
            return localCoefficients_.size();
        }

        //! Return current element, throw if unbound
        const Element &element() const
        {
            return *element_;
        }

        /** \brief Return the geometry type of element the local finite element is defined on (here, a Simplex)
         */
        GeometryType type() const
        {
            return GeometryTypes::simplex(dim);
        }

    private:
        ReducedCubicHermiteTriangleLocalBasis<GV, R> localBasis_;
        ReducedCubicHermiteTriangleLocalCoefficients localCoefficients_;
        ReducedCubicHermiteTriangleLocalInterpolation<ReducedCubicHermiteTriangleLocalBasis<GV, R>> localInterpolation_;

        const Element *element_;
        FieldMatrix<double, 10, 9> C_; // (transposed) basis transformation matrix

        Dune::LagrangeSimplexLocalFiniteElement<R, double, dim, 3> P3LocalFiniteElement_;
    };

    template <typename GV>
    class ReducedCubicHermiteTriangleNode;

    /** \brief Pre-basis for the Reduced cubic Hermite triangle basis
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     * \tparam GV The GridView that the space is defined on
     * \tparam R   Range type used for shape function values
     *
     * The ReducedCubicHermiteTrianglePreBasis can be used to embed a ReducedCubicHermiteTriangleBasis
     * in a larger basis for the construction of product spaces.
     */
    template <typename GV, typename R>
    class ReducedCubicHermiteTrianglePreBasis
    {

    public:
        /** \brief The grid view that the FE space is defined on */
        using GridView = GV;
        using size_type = std::size_t;
        using Node = ReducedCubicHermiteTriangleNode<GV>;

        static constexpr size_type maxMultiIndexSize = 1;
        static constexpr size_type minMultiIndexSize = 1;
        static constexpr size_type multiIndexBufferSize = 1;

        //! Constructor for a given grid view object
        ReducedCubicHermiteTrianglePreBasis(const GridView &gv) : gridView_(gv), mcmgMap_(gv, Layout())
        {
            if (GV::dimension != 2)
                DUNE_THROW(NotImplemented, "ReducedCubicHermiteTriangleBasis only implemented for two-dimensional grids");
        }

        //! Initialize the global indices
        void initializeIndices()
        {
        }

        //! Obtain the grid view that the basis is defined on
        const GridView &gridView() const
        {
            return gridView_;
        }

        //! Update the stored grid view & `MultipleCodimMultipleGeomTypeMapper`, to be called if the grid has changed
        void update(const GridView &gv)
        {
            gridView_ = gv;
            mcmgMap_.update(gridView_);
        }

        /**
         * \brief Create tree node
         */
        Node makeNode() const
        {
            return Node{this};
        }

        //! \brief Total number of reduced cubic Hermite triangle basis functions
        size_type size() const
        {
            return mcmgMap_.size();
        }

        //! Return number of possible values for next position in multi index
        template <class SizePrefix>
        size_type size(const SizePrefix prefix) const
        {
            assert(prefix.size() == 0 || prefix.size() == 1);
            return (prefix.size() == 0) ? size() : 0;
        }

        //! Get the total dimension of the space spanned by this basis
        size_type dimension() const
        {
            return size();
        }

        //! Get the maximal number of DOFs associated to node for any element
        size_type maxNodeSize() const
        {
            return 9;
        }

        template <typename It>
        It indices(const Node &node, It it) const
        {
            const auto &element = node.element();

            // throw if Element is not of simplicial type
            if (not(element.type().isSimplex()))
                DUNE_THROW(NotImplemented, "ReducedCubicHermiteTriangleBasis only implemented for simplicial elements!");

            for (size_type i = 0, end = node.finiteElement().size(); i < end; ++it, ++i)
            {
                Dune::LocalKey localKey = node.finiteElement().localCoefficients().localKey(i);
                *it = {{(size_type)(mcmgMap_.subIndex(element, localKey.subEntity(), localKey.codim())) + (size_type)localKey.index()}};
            }
            return it;
        }

    protected:
        GridView gridView_;

        unsigned int order() const
        {
            return 3;
        }

        //! Obtain the `MultipleCodimMultipleGeomTypeMapper`
        const auto &mcmgMap() const
        {
            return mcmgMap_;
        }

        MultipleCodimMultipleGeomTypeMapper<GridView> mcmgMap_;

    private:
        //! layout function for `MultipleCodimMultipleGeomTypeMapper`
        static auto Layout()
        {
            return [](Dune::GeometryType type, int gridDim)
            {
                return (type.isVertex()) ? 3 : 0;
            };
        }
    };

    template <typename GV>
    class ReducedCubicHermiteTriangleNode : public LeafBasisNode
    {

    public:
        using size_type = std::size_t;
        using Element = typename GV::template Codim<0>::Entity;
        using FiniteElement = ReducedCubicHermiteTriangleLocalFiniteElement<GV, double>;

        ReducedCubicHermiteTriangleNode(const ReducedCubicHermiteTrianglePreBasis<GV> *preBasis)
        : element_(nullptr)
        {
        }

        //! Return current element, throw if unbound
        //     static const Element& element()
        const Element &element() const
        {
            return *element_;
        }

        //! \brief Return the LocalFiniteElement for the element we are bound to
        const FiniteElement &finiteElement() const
        {
            return finiteElement_;
        }

        //! Bind to element.
        void bind(const Element &e)
        {
            element_ = &e;

            if (e.type() != finiteElement_.type())
                DUNE_THROW(Dune::Exception,
                           "Reduced cubic Hermite-elements do not exist for elements of type " << e.type());

            // Here the basis transformation matrix is determined for the current element
            finiteElement_.bind(*element_);

            this->setSize(finiteElement_.size());
        }

    protected:
        unsigned int order() const
        {
            return 3;
        }

        FiniteElement finiteElement_;
        const Element *element_;
    };

    namespace BasisFactory
    {

        /**
         * \brief Create a pre-basis factory that can create a reduced cubic Hermite triangle pre-basis
         *
         * \ingroup FunctionSpaceBasesImplementations
         *
         * \tparam R   The range type of the local basis
         */
        template <typename R = double>
        auto reducedCubicHermiteTriangle()
        {
            return [](const auto &gridView)
            {
                return ReducedCubicHermiteTrianglePreBasis<std::decay_t<decltype(gridView)>, R>(gridView);
            };
        }

    } // end namespace BasisFactory

    /** \brief A global reduced cubic Hermite triangle basis
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     * \tparam GV The GridView that the space is defined on
     * \tparam R The range type of the local basis
     */
    template <typename GV, typename R = double>
    using ReducedCubicHermiteTriangleBasis = DefaultGlobalBasis<ReducedCubicHermiteTrianglePreBasis<GV, R>>;

} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REDUCEDCUBICHERMITETRIANGLEBASIS_HH