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
 *
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
    template<typename GV, typename R, typename MI>
    class ReducedCubicHermiteTriangleLocalFiniteElement;

    template<typename GV, class MI, typename R=double>
    class ReducedCubicHermiteTrianglePreBasis;

    /** \brief LocalBasis class in the sense of dune-localfunctions, presenting the restriction
    * of a ReducedCubicHermiteTriangleBasis to an element
    *
    * \ingroup FunctionSpaceBasesImplementations
    *
    * \tparam GV Grid view that the basis is defined on
    * \tparam R Number type used for function values
    * \tparam MI Type to be used for multi-indices
    */
    template<class GV, class R, class MI>
    class ReducedCubicHermiteTriangleLocalBasis
    {
    friend class ReducedCubicHermiteTriangleLocalFiniteElement<GV,R,MI>;

    using Element = typename GV::template Codim<0>::Entity;  //added

    //   typedef typename MI::ctype MultiIndex;  //added;
    typedef typename GV::ctype D;
    enum {dim = GV::dimension};
    public:

    using grid_view = GV;  //added
    using MultiIndex = MI; //added
    using Range = R;       //added

    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,FieldMatrix<R,1,dim> > Traits;  // same as for Lagrangebasis

    //! \brief Constructor with a given preBasis
    ReducedCubicHermiteTriangleLocalBasis(const ReducedCubicHermiteTrianglePreBasis<GV,MI>& preBasis,         // so wird aktuell preBasis nicht gebraucht!!
                        const ReducedCubicHermiteTriangleLocalFiniteElement<GV,R,MI>& lFE)                    // BSpline nutzt preBasis_.evaluateFunction
    : preBasis_(preBasis),
        lFE_(lFE)
    {}

        //! \brief  Evaluate all shape functions
    //   void evaluateFunction (const FieldVector<D,dim>& in,
    //                          std::vector<FieldVector<R,1> >& out) const
    void evaluateFunction(const typename Traits::DomainType& in,
                          std::vector<typename Traits::RangeType>& out) const
    {                                                                                      // Hier braucht man auswertung.. ( an quadraturpunkten ..als input x )


        //  ----- TO DO ------
        // [Konvention] evaluate... nimmt input auf Referenzelement..
        // könnte hier Mapping davorschalten!
        auto element = lFE_.element();
        auto geometry = element.geometry();
//         std::cout << "element corner(0):"<< geometry.corner(0)<< std::endl;
//         std::cout << "element corner(1):"<< geometry.corner(1)<< std::endl;
//         std::cout << "element corner(2):"<< geometry.corner(2)<< std::endl;

        // map input to global coordinates
        auto x = geometry.global(in);


//         out.resize(size());
          out.resize(9);


        auto C = lFE_.C_;   // 9x10                    besser direkt lFE_.C_multiplizieren?


        FieldVector<R,10> tmp;

        // evaluate standard Monomial-Basis at given input
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

        // Multiply Monomial Basis-evaluation with C....

        C.mv(tmp,out);    // or use C.mtv .. with Matrix C that has not been tranposed....    besser direkt lFE_.C_multiplizieren?
        return;
    }




    /** \brief Evaluate Jacobian of all shape functions
        *
        * \param x Point in the reference simplex where to evaluation the Jacobians
        * \param[out] out The Jacobians of all shape functions at the point x
        */
    //   void evaluateJacobian (const FieldVector<D,dim>& in,
    //                          std::vector<FieldMatrix<D,1,dim> >& out) const
    void evaluateJacobian (const typename Traits::DomainType& in,                   // position
                        std::vector<typename Traits::JacobianType>& out) const      // return value
    {

        //  ----- TO DO ------
        // [Konvention] evaluate... nimmt input auf Referenzelement..
        // könnte hier Mapping davorschalten!
        auto element = lFE_.element();
        auto geometry = element.geometry();
//         std::cout << "element corner(0):"<< geometry.corner(0)<< std::endl;
//         std::cout << "element corner(1):"<< geometry.corner(1)<< std::endl;
//         std::cout << "element corner(2):"<< geometry.corner(2)<< std::endl;

        // map input to global coordinates
        auto x = geometry.global(in);




        auto C = lFE_.C_;
//         printmatrix(std::cout, C, "C ", "--");

        out.resize(size());    //9


    //       std::vector<typename Traits::JacobianType> tmp;
    //       tmp.resize(10);

        FieldVector<double,10> tmp_x(0);
        FieldVector<double,10> tmp_y(0);

        tmp_x[0] = 0.0;                   tmp_y[0] = 0.0;       // -1.0 ??
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

        FieldVector<double,9> out_x;
        FieldVector<double,9> out_y;

        C.mv(tmp_x,out_x);
        C.mv(tmp_y,out_y);

    //       tmp[0][0][0] = 0.0;                   tmp[0][0][1] = -1.0;
    //       tmp[1][0][0] = 1.0;                   tmp[1][0][1] = 0.0;
    //       tmp[2][0][0] = 0.0;                   tmp[2][0][1] = 1.0;
    //       tmp[3][0][0] = 2.0*x[0];              tmp[3][0][1] = 0.0;
    //       tmp[4][0][0] = x[1];                  tmp[4][0][1] = x[0];
    //       tmp[5][0][0] = 0.0;                   tmp[5][0][1] = 2.0*x[1];
    //       tmp[6][0][0] = 3.0*x[0]*x[0];         tmp[6][0][1] = 0.0;
    //       tmp[7][0][0] = 2.0*x[0]*x[1];         tmp[7][0][1] = x[0]*x[0];
    //       tmp[8][0][0] = x[1]*x[1];             tmp[8][0][1] = 2.0*x[0]*x[1];
    //       tmp[9][0][0] = 0.0 ;                  tmp[9][0][1] = 3.0*x[1]*x[1];


//         for (int i=0; i<9; i++)                      // ! ! !
//         {
//             if(abs(out_x[i]) < 1e-12)
//                 out_x[i] = 0.0;
//             if(abs(out_y[i]) < 1e-12)
//                 out_y[i] = 0.0;
//         }

        // output in usual form...
        out[0][0][0] = out_x[0];                   out[0][0][1] = out_y[0];
        out[1][0][0] = out_x[1];                   out[1][0][1] = out_y[1];
        out[2][0][0] = out_x[2];                   out[2][0][1] = out_y[2];
        out[3][0][0] = out_x[3];                   out[3][0][1] = out_y[3];
        out[4][0][0] = out_x[4];                   out[4][0][1] = out_y[4];
        out[5][0][0] = out_x[5];                   out[5][0][1] = out_y[5];
        out[6][0][0] = out_x[6];                   out[6][0][1] = out_y[6];
        out[7][0][0] = out_x[7];                   out[7][0][1] = out_y[7];
        out[8][0][0] = out_x[8];                   out[8][0][1] = out_y[8];
//         out[9][0][0] = out_x[9];                   out[9][0][1] = out_y[9];        // !!!! ?!?!?!
        return;
    }



    /** \brief Polynomial order of the shape functions
    */
    unsigned int order () const
    {
        return 3;
    //     return *preBasis_.order_;
    //     Lagrange: preBasis_->order()
    }


    /** \brief Return the number of basis functions on the current knot span
    */
    std::size_t size() const
    {
        return lFE_.size();
    }
    //   constexpr unsigned int size()
    //   {
    //     unsigned int r = 9;
    //     return r;
    //   }

    //! Return current element, throw if unbound
    const Element& element() const
    {
        return *lFE_.element_;
    }



    private:
    const ReducedCubicHermiteTrianglePreBasis<GV,MI,R>& preBasis_;

    const ReducedCubicHermiteTriangleLocalFiniteElement<GV,R,MI>& lFE_;

    // Coordinates in a single knot span differ from coordinates on the B-spline patch
    // by an affine transformation.  This transformation is stored in offset_ and scaling_.
    //   FieldVector<D,dim>    offset_;                                                                 // Trafo sollte bei ReducedCubicHermiteTriangle nicht nötig sein?
    //   DiagonalMatrix<D,dim> scaling_;
    };




    /** \brief Associations of degrees of freedom to subentities of the reference simplex
    *
    * \tparam dim Dimension of the reference simplex
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

            localKeys_[0] = LocalKey(0,2,0);          // function evaluation - vertex dof
            localKeys_[1] = LocalKey(1,2,0);          // function evaluation - vertex dof
            localKeys_[2] = LocalKey(2,2,0);          // function evaluation - vertex dof
            localKeys_[3] = LocalKey(0,2,1);          // partial_x - derivative-vertex dof
            localKeys_[4] = LocalKey(1,2,1);          // partial_x - derivative-vertex dof
            localKeys_[5] = LocalKey(2,2,1);          // partial_x - derivative-vertex dof
            localKeys_[6] = LocalKey(0,2,2);          // partial_y - derivative-vertex dof
            localKeys_[7] = LocalKey(1,2,2);          // partial_y - derivative-vertex dof
            localKeys_[8] = LocalKey(2,2,2);          // partial_y - derivative-vertex dof
            return;

        }

        //! number of coefficients
        static constexpr std::size_t size ()
        {
        constexpr std::size_t r = 9;
        return r;
        }

        //! get i'th index
        const LocalKey& localKey (std::size_t i) const
        {
        return localKeys_[i];
        }

    private:
        std::vector<LocalKey> localKeys_;

    };






        template<typename GV, typename MI>               // Multiindex needed ??  TODO
        class ReducedCubicHermiteTriangleNode;

    //     template<typename GV, class MI>
    //     class ReducedCubicHermiteTriangleNodeIndexSet;
    //


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    /** \brief Evaluate the degrees of freedom of a reduced cubic Hermite triangle basis
    *
    * \tparam LocalBasis The corresponding set of shape functions
    */
    template<class LocalBasis>
    class ReducedCubicHermiteTriangleLocalInterpolation
    {

    public:
        using GV = typename LocalBasis::grid_view;
        using MI = typename LocalBasis::MultiIndex;
        using R  = typename LocalBasis::Range;

        using Element = typename GV::template Codim<0>::Entity;  //added

        ReducedCubicHermiteTriangleLocalInterpolation(const ReducedCubicHermiteTriangleLocalFiniteElement<GV,R,MI>& lFE)      // BSpline nutzt preBasis_.evaluateFunction
        : lFE_(lFE)
        {}

        using D = typename LocalBasis::Traits::DomainFieldType;


//         friend class ReducedCubicHermiteTriangleLocalFiniteElement<GV,R,MI>;       // needed??

        //! \brief Local interpolation of a function
        template <typename F, typename C>
        void interpolate(const F &f, std::vector<C> &out) const
        {
                DUNE_THROW(NotImplemented, "ReducedCubicHermiteTriangleLocalInterpolation::interpolate");
        }
//         const Element* element_;
        const ReducedCubicHermiteTriangleLocalFiniteElement<GV,R,MI>& lFE_;
    };


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /** \brief LocalFiniteElement in the sense of dune-localfunctions, for the Reduced cubic Hermite triangle basis on Simplex grids
    *
    * \ingroup FunctionSpaceBasesImplementations
    *
    * This class ties together the implementation classes ReducedCubicHermiteTriangleLocalBasis, ReducedCubicHermiteTriangleLocalCoefficients, and ReducedCubicHermiteTriangleLocalInterpolation
    *
    * \tparam D Number type used for domain coordinates
    * \tparam R Number type used for spline function values
    * \tparam MI Global multi-index type.  Only here for technical reasons
    */
    template<class GV, class R, class MI>
    class ReducedCubicHermiteTriangleLocalFiniteElement
    {
        typedef typename GV::ctype D;
        enum {dim = GV::dimension};

        friend class ReducedCubicHermiteTriangleLocalBasis<GV,R,MI>;

        public:
        using Element = typename GV::template Codim<0>::Entity;   // (added) ?

        /** \brief Export various types related to this LocalFiniteElement
        */
        typedef LocalFiniteElementTraits<ReducedCubicHermiteTriangleLocalBasis<GV,R,MI>,
                                         ReducedCubicHermiteTriangleLocalCoefficients<dim>,
                                         ReducedCubicHermiteTriangleLocalInterpolation<ReducedCubicHermiteTriangleLocalBasis<GV,R,MI>>> Traits;
        //   ReducedCubicHermiteTriangleLocalInterpolation<dim,ReducedCubicHermiteTriangleLocalBasis<GV,R,MI> > Traits;      //changed

        /** \brief Constructor with a given ReducedCubicHermiteTriangle basis
        */
        //   ReducedCubicHermiteTriangleLocalFiniteElement(const ReducedCubicHermiteTrianglePreBasis<GV,MI>& preBasis)
        //   : preBasis_(preBasis),
        //     localBasis_(preBasis,*this)
        //   {}

        ReducedCubicHermiteTriangleLocalFiniteElement(const ReducedCubicHermiteTrianglePreBasis<GV,MI>& preBasis)
        : preBasis_(preBasis),
            localBasis_(preBasis,*this),
            element_(nullptr),
            localInterpolation_(*this)          // added
        {}

        /** \brief Bind LocalFiniteElement to a specific element
         *
         *  Compute the basis transformation matrix on each element for the
         *  transformation of the  standard Monomial-Basis to the dual Basis
         *  corresponding to the degrees of freedom (function values and gradients at the nodes).
         */
        void bind(const Element& e)
        {

        //       element_ = e;
        element_ = &e;

        //    C_.setSize(10,9);               // bereits in Klassendefinition deklariert    is this needed?

        //    Dune::Matrix<double> N;         //
        //    N.setSize(10,10);               //  later: (LocalFE.size() +1, dim P3 = 10  )
        //    N = 0;
        FieldMatrix<double,10,10> N(0);

        //    Dune::Matrix<double> b;         //
        //    b.setSize(10,9);               // #Functionals + 1. ReducedCondition = 9 + 1  times #Basis-Functions W_h = 9
        //    b = 0;
        FieldMatrix<double,10,9> b(0);

        auto geometry = e.geometry();        //dereferenzieren * ?
        //    std::cout << "element corner(0):"<< geometry.corner(0)<< std::endl;
        //    std::cout << "element corner(1):"<< geometry.corner(1)<< std::endl;
        //    std::cout << "element corner(2):"<< geometry.corner(2)<< std::endl;
    //

    //     auto gv = preBasis_.gridView();
    //     const auto& indexSet = gv.indexSet();
    //     // Global vertex indices within the grid
    //     auto globalV0 = indexSet.subIndex(e,0,dim);
    //     auto globalV1 = indexSet.subIndex(e,1,dim);


        auto c0 = geometry.corner(0);
        auto c1 = geometry.corner(1);
        auto c2 = geometry.corner(2);

        //    std::cout << "element center:" << geometry.center() << std::endl; //barycenter
        auto center = geometry.center();

        for (int i=0; i<geometry.corners(); i++)
        {
                auto c = geometry.corner(i);

                // first three functionals ... point-evaluations at corners
                N[i][0] = 1.0;            //X_1(1)
                N[i][1] = c[0];           //X_1(x) ...
                N[i][2] = c[1];           //X_1(y)
                N[i][3] = c[0]*c[0];      //X_1(x^2)
                N[i][4] = c[0]*c[1];      //X_1(x*y)
                N[i][5] = c[1]*c[1];      //X_1(y^2)
                N[i][6] = c[0]*c[0]*c[0]; //X_1(x^3)
                N[i][7] = c[0]*c[0]*c[1]; //X_1(x^2*y)
                N[i][8] = c[0]*c[1]*c[1]; //X_1(x*y^2)
                N[i][9] = c[1]*c[1]*c[1]; //X_1(y^3)

                // next three functionals ... point-evaluations of partial_x derivative at corners
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

                // next three functionals ... point-evaluations of partial_y derivative at corners
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

        // "kinematic condition" (removes degree of freedom in the center of the element)
        N[9][0] = 0.0;
        N[9][1] = 3.0*center[0] - c0[0] - c1[0] - c2[0];
        N[9][2] = 3.0*center[1] - c0[1] - c1[1] - c2[1];                  // !!
        N[9][3] = 6.0*center[0]*center[0] - (c0[0]+c1[0]+c2[0])*2.0*center[0];
        N[9][4] = 6.0*center[0]*center[1] - c0[1]*center[0]- c0[0]*center[1]-c1[1]*center[0]-c1[0]*center[1]-c2[1]*center[0]-c2[0]*center[1];
        N[9][5] = 6.0*center[1]*center[1] - (c0[1] + c1[1] + c2[1])*2.0*center[1];
        N[9][6] = 6.0*center[0]*center[0]*center[0] + c0[0]*c0[0]*c0[0] + c1[0]*c1[0]*c1[0] + c2[0]*c2[0]*c2[0] - (c0[0]*c0[0] + c1[0]*c1[0] + c2[0]*c2[0])*3.0*center[0];
        N[9][7] = 6.0*center[0]*center[0]*center[1] - 2.0*center[0]*(c0[0]*c0[1] + c1[0]*c1[1]+ c2[0]*c2[1])
                                            + c0[0]*c0[0]*(c0[1]-center[1]) + c1[0]*c1[0]*(c1[1]-center[1]) + c2[0]*c2[0]*(c2[1]-center[1]);
        N[9][8] = 6.0*center[0]*center[1]*center[1] - 2.0*center[1]*(c0[0]*c0[1]+ c1[0]*c1[1]+ c2[0]*c2[1]) + c0[1]*c0[1]*(c0[0]-center[0])
                                            + c1[1]*c1[1]*(c1[0]-center[0]) + c2[1]*c2[1]*(c2[0]-center[0]);
        N[9][9] = 6.0*center[1]*center[1]*center[1] + c0[1]*c0[1]*c0[1] + c1[1]*c1[1]*c1[1] + c2[1]*c2[1]*c2[1]
                                            - 3.0*center[1]*(c0[1]*c0[1]+c1[1]*c1[1]+c2[1]*c2[1]);    // !!

        // Right-Hand side = Identity Matrix with added zero-row
        b[0][0] = 1.0;
        b[1][1] = 1.0;
        b[2][2] = 1.0;
        b[3][3] = 1.0;
        b[4][4] = 1.0;
        b[5][5] = 1.0;
        b[6][6] = 1.0;
        b[7][7] = 1.0;
        b[8][8] = 1.0;

        N.invert();// Todo: improve
    //     printmatrix(std::cout, N, "N inverse: ", "--");
        FieldMatrix<double,10,9> C;
        C = b.leftmultiply(N);
    //     printmatrix(std::cout, C, "C ", "--");

        // Transpose FieldMatrix C
        const int rows = C_.N();
        const int cols = C_.M();

        for(int i=0; i<rows; ++i)
            for(int j=0; j<cols; ++j)
            {
    //             if( abs(C[j][i]) > 1.0e-12)                          // ! ! ! !
                C_[i][j] = C[j][i];
            }
    //        printmatrix(std::cout, C, "C:", "--");
    //        printmatrix(std::cout, C_, "C_ (Transposed of C):", "--");
        }
        // ------------------------------------------------------------------------------------------------------------------
        // ------------------------------------------------------------------------------------------------------------------
        // ------------------------------------------------------------------------------------------------------------------



        /** \brief Hand out a LocalBasis object */
        const ReducedCubicHermiteTriangleLocalBasis<GV,R,MI>& localBasis() const
        {
            return localBasis_;
        }

        /** \brief Hand out a LocalCoefficients object */
        const ReducedCubicHermiteTriangleLocalCoefficients<dim>& localCoefficients() const
        {
            return localCoefficients_;
        }

        /** \brief Hand out a LocalInterpolation object */
        //   const ReducedCubicHermiteTriangleLocalInterpolation<dim,ReducedCubicHermiteTriangleLocalBasis<GV,R,MI> >& localInterpolation() const    //changed
        const ReducedCubicHermiteTriangleLocalInterpolation<ReducedCubicHermiteTriangleLocalBasis<GV,R,MI> >& localInterpolation() const
        {
            return localInterpolation_;
        }

        //! \brief Number of shape functions in this finite element
        unsigned size () const
        {
            std::size_t r = 9;
            return r;
        }
        //   static constexpr std::size_t size ()
        //   {
        //       return ReducedCubicHermiteTriangleLocalBasis<GV,R,MI>::size();
        //   }
        //   static constexpr unsigned int size ()
        //   {
        //         constexpr unsigned int r = 9;
        //         return r;
        //   }



        //! Return current element, throw if unbound
        const Element& element() const
        {
            return *element_;
        }


        //    static Element& element()  //added
        //   {
        //     return element_;
        //   }

        /** \brief Return the reference element that the local finite element is defined on (here, a Simplex)
        */
        GeometryType type () const
        {
            return GeometryTypes::simplex(dim);
        }



        const ReducedCubicHermiteTrianglePreBasis<GV,MI>& preBasis_;                  // never needed here?!

        ReducedCubicHermiteTriangleLocalBasis<GV,R,MI> localBasis_;
        ReducedCubicHermiteTriangleLocalCoefficients<dim> localCoefficients_;
        //   ReducedCubicHermiteTriangleLocalInterpolation<dim,ReducedCubicHermiteTriangleLocalBasis<GV,R,MI> > localInterpolation_;
        ReducedCubicHermiteTriangleLocalInterpolation<ReducedCubicHermiteTriangleLocalBasis<GV,R,MI> > localInterpolation_;


        //   Dune::Matrix<double> C_;    //BasisTrafoMatrix_;     // !!!!
        //   FieldMatrix<double,10,9> C_;
        FieldMatrix<double,9,10> C_;      // Transposed already!

        //   Element element_;
        const Element* element_;

    };








    /** \brief Pre-basis for the Reduced cubic Hermite triangle basis
    *
    * \ingroup FunctionSpaceBasesImplementations
    *
    * \tparam GV The GridView that the space is defined on
    * \tparam MI Type to be used for multi-indices
    * \tparam R   Range type used for shape function values
    *
    * The ReducedCubicHermiteTrianglePreBasis can be used to embed a ReducedCubicHermiteTriangleBasis
    * in a larger basis for the construction of product spaces.
    */

    template<typename GV, class MI, typename R>
    class ReducedCubicHermiteTrianglePreBasis
    {
    static const int dim = GV::dimension;

    public:

    /** \brief The grid view that the FE space is defined on */
    using GridView = GV;
    using size_type = std::size_t;

    using Node = ReducedCubicHermiteTriangleNode<GV, MI>;

    //   using IndexSet = ReducedCubicHermiteTriangleNodeIndexSet<GV, MI>;
    //    //! Type of created tree node index set. \deprecated
    using IndexSet = Impl::DefaultNodeIndexSet<ReducedCubicHermiteTrianglePreBasis>;     // CHANGE ala Hierarchical

    /** \brief Type used for global numbering of the basis vectors */
    using MultiIndex = MI;

    using SizePrefix = Dune::ReservedVector<size_type, 1>;

    // Type used for function values
    //   using R = double;


    //! Constructor for a given grid view object with compile-time order
    //   ReducedCubicHermiteTrianglePreBasis(const GridView& gv)
    //   : ReducedCubicHermiteTrianglePreBasis(gv, std::numeric_limits<unsigned int>::max())
    //   {}
    ReducedCubicHermiteTrianglePreBasis(const GridView& gv) : gridView_(gv) , mcmgMap_(gv,Layout())
    {}



    //! Initialize the global indices
    void initializeIndices()                                   // Offsets for derivatives needed?
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

    /**
    * \brief Create tree node index set
    *
    * Create an index set suitable for the tree node obtained
    * by makeNode().
    */
    IndexSet makeIndexSet() const
    {
        return IndexSet{*this};
    }


    //! \brief Total number of Reduced cubic Hermite triangle basis functions
    //   unsigned int size () const
    //   {
    //       return 3*((size_type)gridView_.size(dim)); // 3 times the number of all Vertices             BETTER WITH MCMG MAPPER
    //
    //   }
    size_type size() const
    {
        return mcmgMap_.size();
    }



    //! Return number of possible values for next position in multi index
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
        return 9;        //TODO
    }

    template<typename It>
    It indices(const Node& node, It it) const
    {
        for (size_type i = 0, end = node.finiteElement().size() ; i < end ; ++it, ++i)
        {
            Dune::LocalKey localKey = node.finiteElement().localCoefficients().localKey(i);
            const auto& element = node.element();

//             std::cout << "Index: "<< (size_type)(mcmgMap_.subIndex(element,localKey.subEntity(),localKey.codim())) + (size_type)localKey.index()  << std::endl;

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









    template<typename GV, typename MI>
    class ReducedCubicHermiteTriangleNode :
    public LeafBasisNode
    {
    static const int dim = GV::dimension;

    public:

    using size_type = std::size_t;
    using Element = typename GV::template Codim<0>::Entity;
    using FiniteElement = ReducedCubicHermiteTriangleLocalFiniteElement<GV,double,MI>;

    ReducedCubicHermiteTriangleNode(const ReducedCubicHermiteTrianglePreBasis<GV, MI>* preBasis) :
        preBasis_(preBasis),                    // preBasis wird hier nie gebraucht? vgl Lagrange!
        finiteElement_(*preBasis),              //already determines LocalfiniteElement ??             // ????
        element_(nullptr)                     // vgl Lagrange
    {}

    //! Return current element, throw if unbound
    //     static const Element& element()
        const Element& element() const
    {
        return *element_;
    }

    //    Element& element()
    //   {
    //     return element_;
    //   }

    //! \brief Return the LocalFiniteElement for the element we are bound to
    const FiniteElement& finiteElement() const
    {
        return finiteElement_;
    }

    //! Bind to element.
    void bind(const Element& e)                                                    ///////// !!!!!!!!!! TODO
    {
        element_ = &e;                                                                 //changed
    //     element_ = e;                     //eher Referenz?!
    //     auto elementIndex = preBasis_->gridView().indexSet().index(e);                 //global Index of element


        // Hier wird für ein konkretes LocalFiniteElement ein currentKnotSpan bestimmt.
        // Stattdessen bestimmen wir für das FiniteElement die BasisTrafoMatrix C

        if (e.type() != finiteElement_.type())
            DUNE_THROW(Dune::Exception,
                        "Reduced cubic Hermite-elements do not exist for elements of type " << e.type());

        finiteElement_.bind(*element_);

    //     finiteElement_.bind(preBasis_->getIJK(elementIndex,preBasis_->elements_));       // bind@LocalFE takes element index
    //     finiteElement_ = new HierarchicalP2LocalFiniteElement<typename GV::ctype,R,dim>;
    //     finiteElement_.bind(preBasis_->getIJK(elementIndex,preBasis_->elements_));       // not needed ????
        this->setSize(finiteElement_.size());
    }

    protected:

    unsigned int order() const
    {
        return 3;
    }

    const ReducedCubicHermiteTrianglePreBasis<GV,MI>* preBasis_;      // NEEDED ?!?!?!

    FiniteElement finiteElement_;
    const Element* element_;   // (vgl. Lagrangebasis)
    //   Element element_;
    };




    namespace BasisFactory {

    namespace Impl {                                                //TODO REMOVE!

    // template<int k, typename R=double>
    template<typename R=double>
    class ReducedCubicHermiteTrianglePreBasisFactory
    {
    public:
    static const std::size_t requiredMultiIndexSize=1;

    //   ReducedCubicHermiteTrianglePreBasisFactory(const std::vector<double>& knotVector,
    //                          unsigned int order,
    //                          bool makeOpen = true)
    //   : knotVector_(knotVector),
    //     order_(order),
    //     makeOpen_(makeOpen)
    //   {}

    // \brief Constructor for factory with compile-time order
    //   ReducedCubicHermiteTrianglePreBasisFactory() : order_(0)                                 // Why needed ??
    //   {}

    //   // \brief Constructor for factory with run-time order (template argument k is disregarded)
    //   ReducedCubicHermiteTrianglePreBasisFactory(unsigned int order)
    //   : order_(order)
    //   {}

    template<class MultiIndex, class GridView>
    auto makePreBasis(const GridView& gridView) const
    {
    //     return ReducedCubicHermiteTrianglePreBasis<GridView, MultiIndex>(gridView, order);
        return ReducedCubicHermiteTrianglePreBasis<GridView, MultiIndex, R>(gridView);
    }

    // private:
    //   unsigned int order_;

    };

    } // end namespace BasisFactory::Impl

    /**
     * \brief Create a pre-basis factory that can create a reduced cubic Hermite triangle pre-basis
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     * \tparam R   The range type of the local basis
     */

    // template<std::size_t k, typename R=double>
    template<typename R=double>
    auto reducedCubicHermiteTriangle()

    {
    return Impl::ReducedCubicHermiteTrianglePreBasisFactory<R>();
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
    using ReducedCubicHermiteTriangleBasis = DefaultGlobalBasis<ReducedCubicHermiteTrianglePreBasis<GV, FlatMultiIndex<std::size_t>, R> >;


    }   // namespace Functions

}   // namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_REDUCEDCUBICHERMITETRIANGLEBASIS_HH
