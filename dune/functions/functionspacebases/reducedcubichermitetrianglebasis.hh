// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REDUCEDCUBICHERMITETRIANGLEBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REDUCEDCUBICHERMITETRIANGLEBASIS_HH

/** \file
 * \brief A reduced cubic hermite triangle global function space basis
 *
 * This function space is involved in the definition of
 * discrete Kirchhoff Triangles (DKT - Elements) for the
 * discretization of Kirchhoff plate bending problems.
 * See "Braess - Finite Elements 5.th Edition p.335" for more details.
 */

#include <array>
#include <numeric>

/** \todo Don't use this matrix */
#include <dune/common/dynmatrix.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/io.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/common/diagonalmatrix.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/geometry/type.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/grid/common/mcmgmapper.hh>

namespace Dune {

    namespace Functions {

    // A maze of dependencies between the different parts of this.  We need a few forward declarations
    template<typename GV, typename R>
    class ReducedCubicHermiteTriangleLocalFiniteElement;

    template<typename GV, typename R=double>
    class ReducedCubicHermiteTrianglePreBasis;

    /** \brief LocalBasis class in the sense of dune-localfunctions, presenting the restriction
    * of a ReducedCubicHermiteTriangleBasis to an element
    *
    * \ingroup FunctionSpaceBasesImplementations
    *
    * \tparam GV Grid view that the basis is defined on
    * \tparam R Number type used for function values
    */
    template<class GV, class R>
    class ReducedCubicHermiteTriangleLocalBasis
    {
    friend class ReducedCubicHermiteTriangleLocalFiniteElement<GV,R>;

    typedef typename GV::ctype D;
    enum {dim = GV::dimension};
    public:
        using Element = typename GV::template Codim<0>::Entity;
        using grid_view = GV;
        using Range = R;

        //! \brief export type traits for function signature
        typedef LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,FieldMatrix<R,1,dim> > Traits;

        //! \brief Constructor with a given preBasis
        ReducedCubicHermiteTriangleLocalBasis(const ReducedCubicHermiteTrianglePreBasis<GV>& preBasis,
                                              const ReducedCubicHermiteTriangleLocalFiniteElement<GV,R>& lFE)
        : preBasis_(preBasis),
            lFE_(lFE)
        {}

        /** \brief  Evaluate all shape functions
        *  \param in Coordinates where to evaluate the functions, in local coordinates
        *  \param out Output vector of values of shape function values
        */
        void evaluateFunction(const typename Traits::DomainType& in,
                              std::vector<typename Traits::RangeType>& out) const
        {
            auto element = lFE_.element();
            auto geometry = element.geometry();

            // map input to global coordinates
            auto x = geometry.global(in);

            out.resize(size());

            // Evaluate the standard Monomial cubic basis {1,x,y,x^2,...,x^3} at the given input
            FieldVector<R, 10> tmp(0);

            tmp[0] = 1.0;
            tmp[1] = x[0];
            tmp[2] = x[1];
            tmp[3] = x[0]*x[0];
            tmp[4] = x[0]*x[1];
            tmp[5] = x[1]*x[1];
            tmp[6] = x[0]*x[0]*x[0];
            tmp[7] = x[0]*x[0]*x[1];
            tmp[8] = x[0]*x[1]*x[1];
            tmp[9] = x[1]*x[1]*x[1];

            // Multiply Monomial Basis-evaluation with basis transformation matrix.
            // Note that we need to multiply the transposed of C_
            lFE_.C_mtv(tmp, out);
            return;
        }

        /** \brief Evaluate Jacobian of all shape functions
        *  \param x Coordinates where to evaluate the Jacobian, in local coordinates
        *  \param out The Jacobians of all shape functions at the point x
        */
        void evaluateJacobian (const typename Traits::DomainType& in,
                               std::vector<typename Traits::JacobianType>& out) const
        {

            auto element = lFE_.element();
            auto geometry = element.geometry();

            // map input to global coordinates
            auto x = geometry.global(in);

        //       std::vector<typename Traits::JacobianType> tmp;
        //       tmp.resize(10);
            FieldVector<double,10> tmp_x(0);
            FieldVector<double,10> tmp_y(0);

            // Evaluate the partial derivatives of the standard Monomial cubic basis {1,x,y,x^2,...,x^3} at the given input
            tmp_x[0] = 0.0;                   tmp_y[0] = 0.0;
            tmp_x[1] = 1.0;                   tmp_y[1] = 0.0;
            tmp_x[2] = 0.0;                   tmp_y[2] = 1.0;
            tmp_x[3] = 2.0*x[0];              tmp_y[3] = 0.0;
            tmp_x[4] = x[1];                  tmp_y[4] = x[0];
            tmp_x[5] = 0.0;                   tmp_y[5] = 2.0*x[1];
            tmp_x[6] = 3.0*x[0]*x[0];         tmp_y[6] = 0.0;
            tmp_x[7] = 2.0*x[0]*x[1];         tmp_y[7] = x[0]*x[0];
            tmp_x[8] = x[1]*x[1];             tmp_y[8] = 2.0*x[0]*x[1];
            tmp_x[9] = 0.0 ;                  tmp_y[9] = 3.0*x[1]*x[1];

            // Multiply Monomial Basis-evaluation with C....
        //        for (size_t i=0; i<10; i++)
        //         C.mv(tmp[i][0], out[i]);

            // Multiply Jacobian evaluation with basis transformation matrix
            FieldVector<double,9> out_x;
            FieldVector<double,9> out_y;

            lFE_.C_.mtv(tmp_x, out_x);
            lFE_.C_.mtv(tmp_y, out_y);

            // fill output vector
            out.resize(size());
            out[0][0][0] = out_x[0];                   out[0][0][1] = out_y[0];
            out[1][0][0] = out_x[1];                   out[1][0][1] = out_y[1];
            out[2][0][0] = out_x[2];                   out[2][0][1] = out_y[2];
            out[3][0][0] = out_x[3];                   out[3][0][1] = out_y[3];
            out[4][0][0] = out_x[4];                   out[4][0][1] = out_y[4];
            out[5][0][0] = out_x[5];                   out[5][0][1] = out_y[5];
            out[6][0][0] = out_x[6];                   out[6][0][1] = out_y[6];
            out[7][0][0] = out_x[7];                   out[7][0][1] = out_y[7];
            out[8][0][0] = out_x[8];                   out[8][0][1] = out_y[8];

            return;
        }

        //! \brief Polynomial order of the shape functions
        unsigned int order () const
        {
            return preBasis_.order();
        }

        //! \brief Return the number of basis functions
        std::size_t size() const
        {
            return lFE_.size();
        }

        //! Return current element, throw if unbound
        const Element& element() const
        {
            return *lFE_.element_;
        }

    private:
        const ReducedCubicHermiteTrianglePreBasis<GV,R>& preBasis_;
        const ReducedCubicHermiteTriangleLocalFiniteElement<GV,R>& lFE_;
    };



    /** \brief Associations of degrees of freedom to subentities of the reference simplex
    *  \tparam dim Dimension of the reference simplex
    */
    template<unsigned int dim>
    class ReducedCubicHermiteTriangleLocalCoefficients
    {
    public:
        //! \brief Default constructor
        ReducedCubicHermiteTriangleLocalCoefficients ()
        : localKeys_(size())
        {
            if (dim != 2)
                DUNE_THROW(NotImplemented, "ReducedCubicHermiteTriangleLocalCoefficients only implemented for dim=2");

            localKeys_[0] = LocalKey(0,2,0);    // function evaluation - vertex dof
            localKeys_[1] = LocalKey(1,2,0);    // function evaluation - vertex dof
            localKeys_[2] = LocalKey(2,2,0);    // function evaluation - vertex dof
            localKeys_[3] = LocalKey(0,2,1);    // partial_x - derivative-vertex dof
            localKeys_[4] = LocalKey(1,2,1);    // partial_x - derivative-vertex dof
            localKeys_[5] = LocalKey(2,2,1);    // partial_x - derivative-vertex dof
            localKeys_[6] = LocalKey(0,2,2);    // partial_y - derivative-vertex dof
            localKeys_[7] = LocalKey(1,2,2);    // partial_y - derivative-vertex dof
            localKeys_[8] = LocalKey(2,2,2);    // partial_y - derivative-vertex dof
            return;

        }

        //! number of coefficients
        static constexpr std::size_t size ()
        {
            return 9;
        }

        //! get i'th index
        const LocalKey& localKey (std::size_t i) const
        {
        return localKeys_[i];
        }

    private:
        std::vector<LocalKey> localKeys_;
    };


    template<typename GV>
    class ReducedCubicHermiteTriangleNode;


    /** \brief Evaluate the degrees of freedom of a reduced cubic Hermite triangle basis
    * \tparam LocalBasis The corresponding set of shape functions
    */
    template<class LocalBasis>
    class ReducedCubicHermiteTriangleLocalInterpolation
    {
    public:
        using GV      = typename LocalBasis::grid_view;
        using R       = typename LocalBasis::Range;
        using Element = typename LocalBasis::Element;

        ReducedCubicHermiteTriangleLocalInterpolation(const ReducedCubicHermiteTriangleLocalFiniteElement<GV,R>& lFE)      // BSpline nutzt preBasis_.evaluateFunction
        : lFE_(lFE)
        {}

        //! \brief Local interpolation of a function
        template <typename F, typename C>
        void interpolate(const F &f, std::vector<C> &out) const
        {
                DUNE_THROW(NotImplemented, "ReducedCubicHermiteTriangleLocalInterpolation::interpolate");
        }

    private:
        const ReducedCubicHermiteTriangleLocalFiniteElement<GV,R>& lFE_;
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
    template<class GV, class R>
    class ReducedCubicHermiteTriangleLocalFiniteElement
    {
        typedef typename GV::ctype D;
        enum {dim = GV::dimension};

        friend class ReducedCubicHermiteTriangleLocalBasis<GV,R>;

        public:
        using Element = typename GV::template Codim<0>::Entity;

        /** \brief Export various types related to this LocalFiniteElement
        */
        typedef LocalFiniteElementTraits<ReducedCubicHermiteTriangleLocalBasis<GV,R>,
                                         ReducedCubicHermiteTriangleLocalCoefficients<dim>,
                                         ReducedCubicHermiteTriangleLocalInterpolation<ReducedCubicHermiteTriangleLocalBasis<GV,R>>> Traits;


        /** \brief Constructor with a given ReducedCubicHermiteTriangle basis
         */
        ReducedCubicHermiteTriangleLocalFiniteElement(const ReducedCubicHermiteTrianglePreBasis<GV> &preBasis)
            : localBasis_(preBasis, *this),
              element_(nullptr),
              localInterpolation_(*this){}
        /** \brief Bind LocalFiniteElement to a specific element
         *
         *  Compute the basis transformation matrix on each element for the
         *  transformation of the standard cubic Monomial-Basis to the dual Basis
         *  corresponding to the degrees of freedom: function values and gradients at the nodes.
         */
        void bind(const Element& e)
        {
        element_ = &e;
        auto geometry = element_->geometry();

        /*
        * left-hand side matrix N are the degrees of freedom
        * applied to the cubic monomial basis.
        */
        FieldMatrix<double, 10, 10> N(0);

        for (int i=0; i<geometry.corners(); i++)
        {
            auto c = geometry.corner(i);
            // first three degrees of freedom: point-evaluations at corners
            N[i][0] = 1.0;            //(1)
            N[i][1] = c[0];           //(x) ...
            N[i][2] = c[1];           //(y)
            N[i][3] = c[0]*c[0];      //(x^2)
            N[i][4] = c[0]*c[1];      //(x*y)
            N[i][5] = c[1]*c[1];      //(y^2)
            N[i][6] = c[0]*c[0]*c[0]; //(x^3)
            N[i][7] = c[0]*c[0]*c[1]; //(x^2*y)
            N[i][8] = c[0]*c[1]*c[1]; //(x*y^2)
            N[i][9] = c[1]*c[1]*c[1]; //(y^3)
            // next three degrees of freedom: point-evaluations of partial_x derivative at corners
            N[i+3][0] = 0.0;
            N[i+3][1] = 1.0;
            N[i+3][2] = 0.0;
            N[i+3][3] = 2.0*c[0];
            N[i+3][4] = c[1];
            N[i+3][5] = 0.0;
            N[i+3][6] = 3.0*c[0]*c[0];
            N[i+3][7] = 2.0*c[0]*c[1];
            N[i+3][8] = c[1]*c[1];
            N[i+3][9] = 0.0 ;
            // next three degrees of freedom: point-evaluations of partial_y derivative at corners
            N[i+6][0] = 0.0;
            N[i+6][1] = 0.0;
            N[i+6][2] = 1.0;
            N[i+6][3] = 0.0;
            N[i+6][4] = c[0];
            N[i+6][5] = 2.0*c[1];
            N[i+6][6] = 0.0;
            N[i+6][7] = c[0]*c[0];
            N[i+6][8] = 2.0*c[0]*c[1];
            N[i+6][9] = 3.0*c[1]*c[1];
        }

        /*
        *  The following "kinematic condition" removes degree of freedom in the center of the element.
        */
        auto c0 = geometry.corner(0);
        auto c1 = geometry.corner(1);
        auto c2 = geometry.corner(2);
        auto center = geometry.center();
        N[9][0] = 0.0;
        N[9][1] = 3.0*center[0] - c0[0] - c1[0] - c2[0];
        N[9][2] = 3.0*center[1] - c0[1] - c1[1] - c2[1];
        N[9][3] = 6.0*center[0]*center[0] - (c0[0]+c1[0]+c2[0])*2.0*center[0];
        N[9][4] = 6.0*center[0]*center[1] - c0[1]*center[0] - c0[0]*center[1]-c1[1]*center[0] - c1[0]*center[1] - c2[1]*center[0] - c2[0]*center[1];
        N[9][5] = 6.0*center[1]*center[1] - (c0[1] + c1[1] + c2[1])*2.0*center[1];
        N[9][6] = 6.0*center[0]*center[0]*center[0] + c0[0]*c0[0]*c0[0] + c1[0]*c1[0]*c1[0] + c2[0]*c2[0]*c2[0]
                                            - (c0[0]*c0[0] + c1[0]*c1[0] + c2[0]*c2[0])*3.0*center[0];
        N[9][7] = 6.0*center[0]*center[0]*center[1] - 2.0*center[0]*(c0[0]*c0[1] + c1[0]*c1[1]+ c2[0]*c2[1])
                                            + c0[0]*c0[0]*(c0[1]-center[1]) + c1[0]*c1[0]*(c1[1]-center[1]) + c2[0]*c2[0]*(c2[1]-center[1]);
        N[9][8] = 6.0*center[0]*center[1]*center[1] - 2.0*center[1]*(c0[0]*c0[1]+ c1[0]*c1[1]+ c2[0]*c2[1]) + c0[1]*c0[1]*(c0[0]-center[0])
                                            + c1[1]*c1[1]*(c1[0]-center[0]) + c2[1]*c2[1]*(c2[0]-center[0]);
        N[9][9] = 6.0*center[1]*center[1]*center[1] + c0[1]*c0[1]*c0[1] + c1[1]*c1[1]*c1[1] + c2[1]*c2[1]*c2[1]
                                            - 3.0*center[1]*(c0[1]*c0[1]+c1[1]*c1[1]+c2[1]*c2[1]);

        // Right-Hand side is identity matrix (of size 9x9) with added zero-row
        FieldMatrix<double, 10, 9> b(0);
        for(int i=0; i<size(); i++)
            b[i][i] = 1.0;


        N.invert();// Todo: improve
        //printmatrix(std::cout, N, "inverse of N: ", "--");

        C_ = b.leftmultiply(N);
        }

        /** \brief Hand out a LocalBasis object */
        const ReducedCubicHermiteTriangleLocalBasis<GV,R>& localBasis() const
        {
            return localBasis_;
        }

        /** \brief Hand out a LocalCoefficients object */
        const ReducedCubicHermiteTriangleLocalCoefficients<dim>& localCoefficients() const
        {
            return localCoefficients_;
        }

        /** \brief Hand out a LocalInterpolation object */
        const ReducedCubicHermiteTriangleLocalInterpolation<ReducedCubicHermiteTriangleLocalBasis<GV,R> >& localInterpolation() const
        {
            return localInterpolation_;
        }

        // //! \brief Number of shape functions in this finite element
        std::size_t size() const
        {
            return localCoefficients_.size();
        }

        //! Return current element, throw if unbound
        const Element& element() const
        {
            return *element_;
        }

        /** \brief Return the reference element that the local finite element is defined on (here, a Simplex)
        */
        GeometryType type () const
        {
            return GeometryTypes::simplex(dim);
        }

        private:
            ReducedCubicHermiteTriangleLocalBasis<GV,R> localBasis_;
            ReducedCubicHermiteTriangleLocalCoefficients<dim> localCoefficients_;
            ReducedCubicHermiteTriangleLocalInterpolation<ReducedCubicHermiteTriangleLocalBasis<GV,R> > localInterpolation_;

            FieldMatrix<double,10,9> C_;    // (transposed of) basis transformation matrix

            const Element* element_;
    };





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
    template<typename GV, typename R>
    class ReducedCubicHermiteTrianglePreBasis
    {
    static const int dim = GV::dimension;

    public:
        /** \brief The grid view that the FE space is defined on */
        using GridView = GV;
        using size_type = std::size_t;
        using Node = ReducedCubicHermiteTriangleNode<GV>;

        static constexpr size_type maxMultiIndexSize = 1;
        static constexpr size_type minMultiIndexSize = 1;
        static constexpr size_type multiIndexBufferSize = 1;

        //! Constructor for a given grid view object
        ReducedCubicHermiteTrianglePreBasis(const GridView& gv) : gridView_(gv) , mcmgMap_(gv,Layout())
        {}

        //! Initialize the global indices
        void initializeIndices()
        {}

        //! Obtain the grid view that the basis is defined on
        const GridView& gridView() const
        {
            return gridView_;
        }

        //! Update the stored grid view & `MultipleCodimMultipleGeomTypeMapper`, to be called if the grid has changed
        void update (const GridView& gv)
        {
            gridView_ = gv;
            mcmgMap_.update();
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

        template<typename It>
        It indices(const Node& node, It it) const
        {
            for (size_type i = 0, end = node.finiteElement().size() ; i < end ; ++it, ++i)
            {
                Dune::LocalKey localKey = node.finiteElement().localCoefficients().localKey(i);
                const auto& element = node.element();
                *it = {{ (size_type)(mcmgMap_.subIndex(element,localKey.subEntity(),localKey.codim())) + (size_type)localKey.index() }};
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
        const auto& mcmgMap() const
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
                    if (type.isVertex())
                    return 3;
                    if (type.isLine())
                    return 0;
                    if (type.isTriangle())
                    return 0;
                    assert(type.isTetrahedron());
                    return 0;
                };
        }
    };





    template<typename GV>
    class ReducedCubicHermiteTriangleNode :
        public LeafBasisNode
    {
    static const int dim = GV::dimension;

    public:
        using size_type = std::size_t;
        using Element = typename GV::template Codim<0>::Entity;
        using FiniteElement = ReducedCubicHermiteTriangleLocalFiniteElement<GV,double>;

        ReducedCubicHermiteTriangleNode(const ReducedCubicHermiteTrianglePreBasis<GV>* preBasis) :
            finiteElement_(*preBasis),
            element_(nullptr)
        {}

        //! Return current element, throw if unbound
        //     static const Element& element()
        const Element& element() const
        {
            return *element_;
        }

        //! \brief Return the LocalFiniteElement for the element we are bound to
        const FiniteElement& finiteElement() const
        {
            return finiteElement_;
        }

        //! Bind to element.
        void bind(const Element& e)
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
        const Element* element_;
    };


    namespace BasisFactory {

    /**
     * \brief Create a pre-basis factory that can create a reduced cubic Hermite triangle pre-basis
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     * \tparam R   The range type of the local basis
     */
    template<typename R=double>
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
    template<typename GV, typename R=double>
    using ReducedCubicHermiteTriangleBasis = DefaultGlobalBasis<ReducedCubicHermiteTrianglePreBasis<GV,R>>;


    }   // namespace Functions
}   // namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REDUCEDCUBICHERMITETRIANGLEBASIS_HH
