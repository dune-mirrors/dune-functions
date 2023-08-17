// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ARGYRISBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ARGYRISBASIS_HH

#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dune/localfunctions/argyris.hh>

#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lineartransformedlocalfiniteelement.hh>

#include <dune/functions/functionspacebases/morleybasis.hh>
#include <dune/functions/functionspacebases/vectorfloatcmp.hh>

#include <dune/c1elements/tensormatvec.hpp>

namespace Dune{
  namespace Functions{
    namespace Impl{
      /** \brief Class implementing a transformation from pullbacks of reference basis functions to
       * global basisfunctions for Argyris elements. For more information, see
       * c1element/dune/functions/functionspacebases/lineartransformedlocalfinitelement.hh.
       * \tparam R RangeFieldType of the finite element
       */
      template <class R>
      class ArgyrisTransformator
      {
      public:
        /** \brief ElementInformation object
         *  Stores information that specify the directions of derivative dofs of an element
         *  Also stores, if those directions are modified (i.e. not the coordinate axes) and if they
         * are to be used for dirichlet or clamped boundary conditions
         * */

        // TODO include a check, whether strong enforcement of Dirichlet conditions is possible and
        // throw Error on call to isDirichlet if not
        template <class Element>
        class ElementInformation
        {
          using GlobalCoordinate = typename Element::Geometry::GlobalCoordinate;
          using ctype = typename Element::Geometry::ctype;
          static constexpr int dim = Element::Geometry::coorddimension;
          static_assert(dim == 2);
          static_assert(std::is_same<GlobalCoordinate, FieldVector<ctype, dim>>::
                            value); // not sure my code works for other types of global coordinates
        public:
          ElementInformation(): ElementInformation((unsigned short) 0) {}
          ElementInformation(std::bitset<3> orientationBitset)
          {
            for (std::size_t i = 0; i < dim + 1; ++i)
            {
              edgeOrientation_[i] = orientationBitset[i] ? -1. : 1.;
              // default directions are global coordinates
              derivativeDirections_[i] = 0;
              for (std::size_t j = 0; j < dim; ++j)
                derivativeDirections_[i][j][j] = 1.;
            }
          }

          ElementInformation(
              std::bitset<3> orientationBitset,
              std::array<FieldMatrix<ctype, dim, dim>, dim + 1> const &directionsPerVertex,
              std::array<std::bitset<3 * dim>, dim + 1> setTangential)
              : derivativeDirections_(directionsPerVertex), booleanInformation_(setTangential)
          {
            for (std::size_t i = 0; i < dim + 1; ++i)
            {
              edgeOrientation_[i] = orientationBitset[i] ? -1. : 1.;
            }
          }

          std::array<R, dim + 1> const &getEdgeOrientation() const { return edgeOrientation_; }

          /**
           * \brief return whether the some direction of vertex"th vertex was modified
           * \param vertex
           * \return bool
           * */
          bool isModified(std::size_t vertex) const
          {
            bool ret = false;
            for (std::size_t i = 0; i < dim; ++i)
              ret |= booleanInformation_[vertex][i];
            return ret;
          }
          /**
           * \brief return whether the "direction"th direction of "vertex"th vertex is a Dirichlet
           * dof (i.e. is a tangential to this element)
           * \param vertex
           * \param direction
           * \return bool
           */
          bool isTangential(std::size_t vertex, std::size_t direction) const
          {
            return booleanInformation_[vertex][dim + direction]
                && booleanInformation_[vertex][2 * dim + direction];
          }

          /**
           * \brief return whether the "direction"th direction of "vertex"th vertex is a normal to
           * this element \param vertex \param direction \return bool
           */
          bool isNormal(std::size_t vertex, std::size_t direction) const
          {
            return !booleanInformation_[vertex][dim + direction]
                && booleanInformation_[vertex][2 * dim + direction];
          }

          /**
           * \brief return the dof corresponding to localKey is a Dirichlet
           * dof (i.e. is a tangential to this very element)
           * \param localKey
           *
           * \return bool
           */
          bool isDirichlet(LocalKey localKey) const
          {
            if (localKey.codim() == 1)
              return false;
            else if (localKey.codim() == 2)
            {
              if (localKey.index() == 0) // value -> always Dirichlet
                return true;
              else if (localKey.index() == 1
                || localKey.index() == 3) // first direction, first or second derivative
                return isTangential(localKey.subEntity(), 0);
              else if (localKey.index() == 2
                || localKey.index() == 5) // second direction, first or second derivative
                return isTangential(localKey.subEntity(), 1);
              else if (localKey.index() == 4) // cross derivative, never Dirichlet
                return false;
            }

            DUNE_THROW(NotImplemented, "Invalid LocalKey");
          }

          /**
           * \brief return the dof corresponding to localKey is a Clamped Dof
           * \param localKey
           *
           * \return bool
           */
          bool isClamped(LocalKey localKey) const
          {
            // normal derivate at edge midpoint
            if (localKey.codim() == 1)
              return true;
            else if (localKey.codim() == 2)
            {
              if (localKey.index() == 0 || localKey.index() == 1 || localKey.index() == 2
                  || localKey.index() == 4) // value  / gradient / cross derivative-> always clamped
                return true;
              else if (localKey.index() == 3) // first direction, second derivative, clamped if
                                              // tangential to this element
                return isTangential(localKey.subEntity(), 0);
              else if (localKey.index() == 5) // second direction, second derivative, clamped if
                                              // tangential to this element
                return isTangential(localKey.subEntity(), 1);
              else
                DUNE_THROW(NotImplemented, "Invalid LocalKey");
            }

            DUNE_THROW(NotImplemented, "Invalid LocalKey");
          }

          /**
           * \brief Get the Derivative Directions object, and array storing the directions of
           * derivative dofs as rows of a fieldmatrix for each vertex
           *
           * \return std::array<FieldMatrix<ctype, dim, dim>, dim + 1> const&
           */
          std::array<FieldMatrix<ctype, dim, dim>, dim + 1> const &getDerivativeDirections() const
          {
            return derivativeDirections_;
          }

        private:
          std::array<R, dim + 1> edgeOrientation_;
          // For every corner of Element, dim directions as rows(!) of a fieldmatrix
          std::array<FieldMatrix<ctype, dim, dim>, dim + 1> derivativeDirections_;
          std::array<std::bitset<3 * dim>, dim + 1>
              booleanInformation_; // first dim bits encode whether direction was modified,
                                   // second dim bits encode whether the direction is tangential to
                                   // the grid, third dim bits encode whether this direction is
                                   // tangential/normal to this very element
        };

        ArgyrisTransformator(): mat_(21, 21, BCRSMatrix<R>::random) { setupMatrix(); }
        /**
         * \brief Binding Method. This Method fills the Transformation matrix, including direction
         * of derivative dofs and orientation of normals
         *
         * \tparam Element Type of Element we are bound to
         * \param element
         * \param elementInfo ElementInformation object
         */
        template <class Element>
        void bind(Element const &element, ElementInformation<Element> const &elementInfo)
        {

          if (element.geometry().affine())
            fillMatrix(element.geometry(), elementInfo);
          else
            DUNE_THROW(Dune::NotImplemented,
                       "Argyris Element is only implemented for affine transformations");
        }

        /**
         * \brief Applies the transformation onto a vector. In constrast to the transformations used
         * by GlobalValuedLocalFiniteElement, this method can be called for shapevalues, Jacobians
         * and Hessians
         *
         * \tparam Values Vector type
         * \param values
         */
        template <class Values>
        void apply(Values &values) const
        {
          Values tmp = values; // needs to be deep copy
          mat_.mv(tmp, values);
        }

      private:
        void setupMatrix()
        { // created with sympy, see argyrisTransformationMatrix.py in module directory
          mat_.setrowsize(0, 3);
          mat_.setrowsize(1, 4);
          mat_.setrowsize(2, 4);
          mat_.setrowsize(3, 5);
          mat_.setrowsize(4, 5);
          mat_.setrowsize(5, 5);
          mat_.setrowsize(6, 3);
          mat_.setrowsize(7, 4);
          mat_.setrowsize(8, 4);
          mat_.setrowsize(9, 5);
          mat_.setrowsize(10, 5);
          mat_.setrowsize(11, 5);
          mat_.setrowsize(12, 3);
          mat_.setrowsize(13, 4);
          mat_.setrowsize(14, 4);
          mat_.setrowsize(15, 5);
          mat_.setrowsize(16, 5);
          mat_.setrowsize(17, 5);
          mat_.setrowsize(18, 1);
          mat_.setrowsize(19, 1);
          mat_.setrowsize(20, 1);
          mat_.endrowsizes();
          mat_.addindex(0, 0);
          mat_.addindex(0, 18);
          mat_.addindex(0, 19);
          mat_.addindex(1, 1);
          mat_.addindex(1, 2);
          mat_.addindex(1, 18);
          mat_.addindex(1, 19);
          mat_.addindex(2, 1);
          mat_.addindex(2, 2);
          mat_.addindex(2, 18);
          mat_.addindex(2, 19);
          mat_.addindex(3, 3);
          mat_.addindex(3, 4);
          mat_.addindex(3, 5);
          mat_.addindex(3, 18);
          mat_.addindex(3, 19);
          mat_.addindex(4, 3);
          mat_.addindex(4, 4);
          mat_.addindex(4, 5);
          mat_.addindex(4, 18);
          mat_.addindex(4, 19);
          mat_.addindex(5, 3);
          mat_.addindex(5, 4);
          mat_.addindex(5, 5);
          mat_.addindex(5, 18);
          mat_.addindex(5, 19);
          mat_.addindex(6, 6);
          mat_.addindex(6, 18);
          mat_.addindex(6, 20);
          mat_.addindex(7, 7);
          mat_.addindex(7, 8);
          mat_.addindex(7, 18);
          mat_.addindex(7, 20);
          mat_.addindex(8, 7);
          mat_.addindex(8, 8);
          mat_.addindex(8, 18);
          mat_.addindex(8, 20);
          mat_.addindex(9, 9);
          mat_.addindex(9, 10);
          mat_.addindex(9, 11);
          mat_.addindex(9, 18);
          mat_.addindex(9, 20);
          mat_.addindex(10, 9);
          mat_.addindex(10, 10);
          mat_.addindex(10, 11);
          mat_.addindex(10, 18);
          mat_.addindex(10, 20);
          mat_.addindex(11, 9);
          mat_.addindex(11, 10);
          mat_.addindex(11, 11);
          mat_.addindex(11, 18);
          mat_.addindex(11, 20);
          mat_.addindex(12, 12);
          mat_.addindex(12, 19);
          mat_.addindex(12, 20);
          mat_.addindex(13, 13);
          mat_.addindex(13, 14);
          mat_.addindex(13, 19);
          mat_.addindex(13, 20);
          mat_.addindex(14, 13);
          mat_.addindex(14, 14);
          mat_.addindex(14, 19);
          mat_.addindex(14, 20);
          mat_.addindex(15, 15);
          mat_.addindex(15, 16);
          mat_.addindex(15, 17);
          mat_.addindex(15, 19);
          mat_.addindex(15, 20);
          mat_.addindex(16, 15);
          mat_.addindex(16, 16);
          mat_.addindex(16, 17);
          mat_.addindex(16, 19);
          mat_.addindex(16, 20);
          mat_.addindex(17, 15);
          mat_.addindex(17, 16);
          mat_.addindex(17, 17);
          mat_.addindex(17, 19);
          mat_.addindex(17, 20);
          mat_.addindex(18, 18);
          mat_.addindex(19, 19);
          mat_.addindex(20, 20);
          mat_.endindices();
        }

        template <class Element>
        void fillMatrix(typename Element::Geometry const &geometry,
                        ElementInformation<Element> const &elementInfo)
        {
          std::array<R, 3> const &edgeOrientation = elementInfo.getEdgeOrientation();
          // get geometrical information
          FieldVector<R, 3> l;

          std::array<Dune::FieldVector<R, 2>, 3> referenceTangents;
          std::array<Dune::FieldVector<R, 2>, 3> globalTangents;

          // By default, edges point from the vertex with the smaller index
          // to the vertex with the larger index. Note that the B matrix is invariant of
          // orientation, since the -1 s of the reference normals/tangents cancel the -1 s
          // from the global normals/tangents normalize

          // get local and global Tangents
          auto refElement
              = Dune::referenceElement<typename Element::Geometry::ctype, 2>(geometry.type());
          for (std::size_t i = 0; i < 3; ++i)
          {
            std::size_t lower = (i == 2) ? 1 : 0;
            std::size_t upper = (i == 0) ? 1 : 2;
            auto edge = refElement.position(upper, 2) - refElement.position(lower, 2);

            referenceTangents[i] = edge / edge.two_norm();

            auto globalEdge = geometry.global(refElement.position(upper, 2))
                            - geometry.global(refElement.position(lower, 2));

            l[i] = globalEdge.two_norm();
            globalTangents[i] = globalEdge / l[i];
          }

          // fill G matrices and multiply to B matrix
          std::array<FieldMatrix<R, 2, 2>, 3> referenceG;
          std::array<FieldMatrix<R, 2, 2>, 3> globalG;
          std::array<Dune::FieldVector<R, 3>, 3> tau;

          std::array<Dune::FieldMatrix<R, 2, 2>, 3> b;

          for (std::size_t i = 0; i < 3; ++i)
          {
            referenceG[i][0] = FieldVector<R, 2>{-referenceTangents[i][1], referenceTangents[i][0]};
            referenceG[i][1] = referenceTangents[i];
            globalG[i][0] = FieldVector<R, 2>{-globalTangents[i][1], globalTangents[i][0]};
            globalG[i][1] = globalTangents[i];
            tau[i] = FieldVector<R, 3>{globalTangents[i][0] * globalTangents[i][0],
                                       2. * globalTangents[i][0] * globalTangents[i][1],
                                       globalTangents[i][1] * globalTangents[i][1]};
          }

          // Multiply jacobian with directional matrix
          auto const &directions = elementInfo.getDerivativeDirections();

          // Create directional, Jacobian and Theta Matrices
          std::array<FieldMatrix<R, 3, 3>, 3> theta;
          std::array<FieldMatrix<R, 3, 3>, 3> thetaDir;
          std::array<FieldMatrix<R, 2, 2>, 3> dir;
          std::array<FieldMatrix<R, 2, 2>, 3> jacobian;
          for (int i = 0; i < geometry.corners(); ++i)
          {
            b[i] = globalG[i] * transpose(geometry.jacobianTransposed(refElement.position(i, 2)))
                 * referenceG[i];

            dir[i] = directions[i];
            dir[i].invert();
            jacobian[i] = geometry.jacobianTransposed(refElement.position(i, 2));
            theta[i][0][0] = jacobian[i][0][0] * jacobian[i][0][0];
            theta[i][0][1] = jacobian[i][0][0] * jacobian[i][1][0];
            theta[i][0][2] = jacobian[i][1][0] * jacobian[i][1][0];

            theta[i][1][0] = 2. * jacobian[i][0][0] * jacobian[i][0][1];
            theta[i][1][1]
                = jacobian[i][0][0] * jacobian[i][1][1] + jacobian[i][1][0] * jacobian[i][0][1];
            theta[i][1][2] = 2. * jacobian[i][1][0] * jacobian[i][1][1];

            theta[i][2][0] = jacobian[i][0][1] * jacobian[i][0][1];
            theta[i][2][1] = jacobian[i][0][1] * jacobian[i][1][1];
            theta[i][2][2] = jacobian[i][1][1] * jacobian[i][1][1];

            // directional matrix for hessian
            thetaDir[i][0][0] = dir[i][0][0] * dir[i][0][0];
            thetaDir[i][0][1] = dir[i][0][0] * dir[i][1][0];
            thetaDir[i][0][2] = dir[i][1][0] * dir[i][1][0];

            thetaDir[i][1][0] = 2. * dir[i][0][0] * dir[i][0][1];
            thetaDir[i][1][1] = dir[i][0][0] * dir[i][1][1] + dir[i][1][0] * dir[i][0][1];
            thetaDir[i][1][2] = 2. * dir[i][1][0] * dir[i][1][1];

            thetaDir[i][2][0] = dir[i][0][1] * dir[i][0][1];
            thetaDir[i][2][1] = dir[i][0][1] * dir[i][1][1];
            thetaDir[i][2][2] = dir[i][1][1] * dir[i][1][1];
          }

          auto &[b_0, b_1, b_2] = b; // compability with sympy code below
          auto &[dir_0, dir_1, dir_2] = dir;
          auto &[J_0, J_1, J_2] = jacobian;
          auto &[theta_0, theta_1, theta_2] = theta;
          auto &[thetaDir_0, thetaDir_1, thetaDir_2] = thetaDir;

          std::array<Dune::FieldVector<R, 2>, 3> const &t = globalTangents;

          // created with sympy, see argyrisTransformationMatrix.py in module directory
          // note that J_i that is the transposed Jacobian, but theta_i is already incorporating
          // this and corresponds to the nontransposed Jacobian applied twice, same for direction
          // matrices
          // TODO rewrite all of this once there is a geometry.jacobian(x) method

          mat_[0][0] = 1;
          mat_[0][18] = -15.0 / 8.0 * b_0[1][0] / l[0];
          mat_[0][19] = -15.0 / 8.0 * b_1[1][0] / l[1];
          mat_[1][1] = J_0[0][0] * dir_0[0][0] + J_0[0][1] * dir_0[1][0];
          mat_[1][2] = J_0[1][0] * dir_0[0][0] + J_0[1][1] * dir_0[1][0];
          mat_[1][18]
              = (-0.4375 * dir_0[0][0] * t[0][0] - 0.4375 * dir_0[1][0] * t[0][1]) * b_0[1][0];
          mat_[1][19]
              = (-0.4375 * dir_0[0][0] * t[1][0] - 0.4375 * dir_0[1][0] * t[1][1]) * b_1[1][0];
          mat_[2][1] = J_0[0][0] * dir_0[0][1] + J_0[0][1] * dir_0[1][1];
          mat_[2][2] = J_0[1][0] * dir_0[0][1] + J_0[1][1] * dir_0[1][1];
          mat_[2][18]
              = (-0.4375 * dir_0[0][1] * t[0][0] - 0.4375 * dir_0[1][1] * t[0][1]) * b_0[1][0];
          mat_[2][19]
              = (-0.4375 * dir_0[0][1] * t[1][0] - 0.4375 * dir_0[1][1] * t[1][1]) * b_1[1][0];
          mat_[3][3] = thetaDir_0[0][0] * theta_0[0][0] + thetaDir_0[0][1] * theta_0[1][0]
                     + thetaDir_0[0][2] * theta_0[2][0];
          mat_[3][4] = thetaDir_0[0][0] * theta_0[0][1] + thetaDir_0[0][1] * theta_0[1][1]
                     + thetaDir_0[0][2] * theta_0[2][1];
          mat_[3][5] = thetaDir_0[0][0] * theta_0[0][2] + thetaDir_0[0][1] * theta_0[1][2]
                     + thetaDir_0[0][2] * theta_0[2][2];
          mat_[3][18] = (-1.0 / 32.0 * l[0] * tau[0][0] * thetaDir_0[0][0]
                         - 1.0 / 32.0 * l[0] * tau[0][1] * thetaDir_0[0][1]
                         - 1.0 / 32.0 * l[0] * tau[0][2] * thetaDir_0[0][2])
                      * b_0[1][0];
          mat_[3][19] = (-1.0 / 32.0 * l[1] * tau[1][0] * thetaDir_0[0][0]
                         - 1.0 / 32.0 * l[1] * tau[1][1] * thetaDir_0[0][1]
                         - 1.0 / 32.0 * l[1] * tau[1][2] * thetaDir_0[0][2])
                      * b_1[1][0];
          mat_[4][3] = thetaDir_0[1][0] * theta_0[0][0] + thetaDir_0[1][1] * theta_0[1][0]
                     + thetaDir_0[1][2] * theta_0[2][0];
          mat_[4][4] = thetaDir_0[1][0] * theta_0[0][1] + thetaDir_0[1][1] * theta_0[1][1]
                     + thetaDir_0[1][2] * theta_0[2][1];
          mat_[4][5] = thetaDir_0[1][0] * theta_0[0][2] + thetaDir_0[1][1] * theta_0[1][2]
                     + thetaDir_0[1][2] * theta_0[2][2];
          mat_[4][18] = (-1.0 / 32.0 * l[0] * tau[0][0] * thetaDir_0[1][0]
                         - 1.0 / 32.0 * l[0] * tau[0][1] * thetaDir_0[1][1]
                         - 1.0 / 32.0 * l[0] * tau[0][2] * thetaDir_0[1][2])
                      * b_0[1][0];
          mat_[4][19] = (-1.0 / 32.0 * l[1] * tau[1][0] * thetaDir_0[1][0]
                         - 1.0 / 32.0 * l[1] * tau[1][1] * thetaDir_0[1][1]
                         - 1.0 / 32.0 * l[1] * tau[1][2] * thetaDir_0[1][2])
                      * b_1[1][0];
          mat_[5][3] = thetaDir_0[2][0] * theta_0[0][0] + thetaDir_0[2][1] * theta_0[1][0]
                     + thetaDir_0[2][2] * theta_0[2][0];
          mat_[5][4] = thetaDir_0[2][0] * theta_0[0][1] + thetaDir_0[2][1] * theta_0[1][1]
                     + thetaDir_0[2][2] * theta_0[2][1];
          mat_[5][5] = thetaDir_0[2][0] * theta_0[0][2] + thetaDir_0[2][1] * theta_0[1][2]
                     + thetaDir_0[2][2] * theta_0[2][2];
          mat_[5][18] = (-1.0 / 32.0 * l[0] * tau[0][0] * thetaDir_0[2][0]
                         - 1.0 / 32.0 * l[0] * tau[0][1] * thetaDir_0[2][1]
                         - 1.0 / 32.0 * l[0] * tau[0][2] * thetaDir_0[2][2])
                      * b_0[1][0];
          mat_[5][19] = (-1.0 / 32.0 * l[1] * tau[1][0] * thetaDir_0[2][0]
                         - 1.0 / 32.0 * l[1] * tau[1][1] * thetaDir_0[2][1]
                         - 1.0 / 32.0 * l[1] * tau[1][2] * thetaDir_0[2][2])
                      * b_1[1][0];
          mat_[6][6] = 1;
          mat_[6][18] = (15.0 / 8.0) * b_0[1][0] / l[0];
          mat_[6][20] = -15.0 / 8.0 * b_2[1][0] / l[2];
          mat_[7][7] = J_1[0][0] * dir_1[0][0] + J_1[0][1] * dir_1[1][0];
          mat_[7][8] = J_1[1][0] * dir_1[0][0] + J_1[1][1] * dir_1[1][0];
          mat_[7][18]
              = (-0.4375 * dir_1[0][0] * t[0][0] - 0.4375 * dir_1[1][0] * t[0][1]) * b_0[1][0];
          mat_[7][20]
              = (-0.4375 * dir_1[0][0] * t[2][0] - 0.4375 * dir_1[1][0] * t[2][1]) * b_2[1][0];
          mat_[8][7] = J_1[0][0] * dir_1[0][1] + J_1[0][1] * dir_1[1][1];
          mat_[8][8] = J_1[1][0] * dir_1[0][1] + J_1[1][1] * dir_1[1][1];
          mat_[8][18]
              = (-0.4375 * dir_1[0][1] * t[0][0] - 0.4375 * dir_1[1][1] * t[0][1]) * b_0[1][0];
          mat_[8][20]
              = (-0.4375 * dir_1[0][1] * t[2][0] - 0.4375 * dir_1[1][1] * t[2][1]) * b_2[1][0];
          mat_[9][9] = thetaDir_1[0][0] * theta_1[0][0] + thetaDir_1[0][1] * theta_1[1][0]
                     + thetaDir_1[0][2] * theta_1[2][0];
          mat_[9][10] = thetaDir_1[0][0] * theta_1[0][1] + thetaDir_1[0][1] * theta_1[1][1]
                      + thetaDir_1[0][2] * theta_1[2][1];
          mat_[9][11] = thetaDir_1[0][0] * theta_1[0][2] + thetaDir_1[0][1] * theta_1[1][2]
                      + thetaDir_1[0][2] * theta_1[2][2];
          mat_[9][18] = ((1.0 / 32.0) * l[0] * tau[0][0] * thetaDir_1[0][0]
                         + (1.0 / 32.0) * l[0] * tau[0][1] * thetaDir_1[0][1]
                         + (1.0 / 32.0) * l[0] * tau[0][2] * thetaDir_1[0][2])
                      * b_0[1][0];
          mat_[9][20] = (-1.0 / 32.0 * l[2] * tau[2][0] * thetaDir_1[0][0]
                         - 1.0 / 32.0 * l[2] * tau[2][1] * thetaDir_1[0][1]
                         - 1.0 / 32.0 * l[2] * tau[2][2] * thetaDir_1[0][2])
                      * b_2[1][0];
          mat_[10][9] = thetaDir_1[1][0] * theta_1[0][0] + thetaDir_1[1][1] * theta_1[1][0]
                      + thetaDir_1[1][2] * theta_1[2][0];
          mat_[10][10] = thetaDir_1[1][0] * theta_1[0][1] + thetaDir_1[1][1] * theta_1[1][1]
                       + thetaDir_1[1][2] * theta_1[2][1];
          mat_[10][11] = thetaDir_1[1][0] * theta_1[0][2] + thetaDir_1[1][1] * theta_1[1][2]
                       + thetaDir_1[1][2] * theta_1[2][2];
          mat_[10][18] = ((1.0 / 32.0) * l[0] * tau[0][0] * thetaDir_1[1][0]
                          + (1.0 / 32.0) * l[0] * tau[0][1] * thetaDir_1[1][1]
                          + (1.0 / 32.0) * l[0] * tau[0][2] * thetaDir_1[1][2])
                       * b_0[1][0];
          mat_[10][20] = (-1.0 / 32.0 * l[2] * tau[2][0] * thetaDir_1[1][0]
                          - 1.0 / 32.0 * l[2] * tau[2][1] * thetaDir_1[1][1]
                          - 1.0 / 32.0 * l[2] * tau[2][2] * thetaDir_1[1][2])
                       * b_2[1][0];
          mat_[11][9] = thetaDir_1[2][0] * theta_1[0][0] + thetaDir_1[2][1] * theta_1[1][0]
                      + thetaDir_1[2][2] * theta_1[2][0];
          mat_[11][10] = thetaDir_1[2][0] * theta_1[0][1] + thetaDir_1[2][1] * theta_1[1][1]
                       + thetaDir_1[2][2] * theta_1[2][1];
          mat_[11][11] = thetaDir_1[2][0] * theta_1[0][2] + thetaDir_1[2][1] * theta_1[1][2]
                       + thetaDir_1[2][2] * theta_1[2][2];
          mat_[11][18] = ((1.0 / 32.0) * l[0] * tau[0][0] * thetaDir_1[2][0]
                          + (1.0 / 32.0) * l[0] * tau[0][1] * thetaDir_1[2][1]
                          + (1.0 / 32.0) * l[0] * tau[0][2] * thetaDir_1[2][2])
                       * b_0[1][0];
          mat_[11][20] = (-1.0 / 32.0 * l[2] * tau[2][0] * thetaDir_1[2][0]
                          - 1.0 / 32.0 * l[2] * tau[2][1] * thetaDir_1[2][1]
                          - 1.0 / 32.0 * l[2] * tau[2][2] * thetaDir_1[2][2])
                       * b_2[1][0];
          mat_[12][12] = 1;
          mat_[12][19] = (15.0 / 8.0) * b_1[1][0] / l[1];
          mat_[12][20] = (15.0 / 8.0) * b_2[1][0] / l[2];
          mat_[13][13] = J_2[0][0] * dir_2[0][0] + J_2[0][1] * dir_2[1][0];
          mat_[13][14] = J_2[1][0] * dir_2[0][0] + J_2[1][1] * dir_2[1][0];
          mat_[13][19]
              = (-0.4375 * dir_2[0][0] * t[1][0] - 0.4375 * dir_2[1][0] * t[1][1]) * b_1[1][0];
          mat_[13][20]
              = (-0.4375 * dir_2[0][0] * t[2][0] - 0.4375 * dir_2[1][0] * t[2][1]) * b_2[1][0];
          mat_[14][13] = J_2[0][0] * dir_2[0][1] + J_2[0][1] * dir_2[1][1];
          mat_[14][14] = J_2[1][0] * dir_2[0][1] + J_2[1][1] * dir_2[1][1];
          mat_[14][19]
              = (-0.4375 * dir_2[0][1] * t[1][0] - 0.4375 * dir_2[1][1] * t[1][1]) * b_1[1][0];
          mat_[14][20]
              = (-0.4375 * dir_2[0][1] * t[2][0] - 0.4375 * dir_2[1][1] * t[2][1]) * b_2[1][0];
          mat_[15][15] = thetaDir_2[0][0] * theta_2[0][0] + thetaDir_2[0][1] * theta_2[1][0]
                       + thetaDir_2[0][2] * theta_2[2][0];
          mat_[15][16] = thetaDir_2[0][0] * theta_2[0][1] + thetaDir_2[0][1] * theta_2[1][1]
                       + thetaDir_2[0][2] * theta_2[2][1];
          mat_[15][17] = thetaDir_2[0][0] * theta_2[0][2] + thetaDir_2[0][1] * theta_2[1][2]
                       + thetaDir_2[0][2] * theta_2[2][2];
          mat_[15][19] = ((1.0 / 32.0) * l[1] * tau[1][0] * thetaDir_2[0][0]
                          + (1.0 / 32.0) * l[1] * tau[1][1] * thetaDir_2[0][1]
                          + (1.0 / 32.0) * l[1] * tau[1][2] * thetaDir_2[0][2])
                       * b_1[1][0];
          mat_[15][20] = ((1.0 / 32.0) * l[2] * tau[2][0] * thetaDir_2[0][0]
                          + (1.0 / 32.0) * l[2] * tau[2][1] * thetaDir_2[0][1]
                          + (1.0 / 32.0) * l[2] * tau[2][2] * thetaDir_2[0][2])
                       * b_2[1][0];
          mat_[16][15] = thetaDir_2[1][0] * theta_2[0][0] + thetaDir_2[1][1] * theta_2[1][0]
                       + thetaDir_2[1][2] * theta_2[2][0];
          mat_[16][16] = thetaDir_2[1][0] * theta_2[0][1] + thetaDir_2[1][1] * theta_2[1][1]
                       + thetaDir_2[1][2] * theta_2[2][1];
          mat_[16][17] = thetaDir_2[1][0] * theta_2[0][2] + thetaDir_2[1][1] * theta_2[1][2]
                       + thetaDir_2[1][2] * theta_2[2][2];
          mat_[16][19] = ((1.0 / 32.0) * l[1] * tau[1][0] * thetaDir_2[1][0]
                          + (1.0 / 32.0) * l[1] * tau[1][1] * thetaDir_2[1][1]
                          + (1.0 / 32.0) * l[1] * tau[1][2] * thetaDir_2[1][2])
                       * b_1[1][0];
          mat_[16][20] = ((1.0 / 32.0) * l[2] * tau[2][0] * thetaDir_2[1][0]
                          + (1.0 / 32.0) * l[2] * tau[2][1] * thetaDir_2[1][1]
                          + (1.0 / 32.0) * l[2] * tau[2][2] * thetaDir_2[1][2])
                       * b_2[1][0];
          mat_[17][15] = thetaDir_2[2][0] * theta_2[0][0] + thetaDir_2[2][1] * theta_2[1][0]
                       + thetaDir_2[2][2] * theta_2[2][0];
          mat_[17][16] = thetaDir_2[2][0] * theta_2[0][1] + thetaDir_2[2][1] * theta_2[1][1]
                       + thetaDir_2[2][2] * theta_2[2][1];
          mat_[17][17] = thetaDir_2[2][0] * theta_2[0][2] + thetaDir_2[2][1] * theta_2[1][2]
                       + thetaDir_2[2][2] * theta_2[2][2];
          mat_[17][19] = ((1.0 / 32.0) * l[1] * tau[1][0] * thetaDir_2[2][0]
                          + (1.0 / 32.0) * l[1] * tau[1][1] * thetaDir_2[2][1]
                          + (1.0 / 32.0) * l[1] * tau[1][2] * thetaDir_2[2][2])
                       * b_1[1][0];
          mat_[17][20] = ((1.0 / 32.0) * l[2] * tau[2][0] * thetaDir_2[2][0]
                          + (1.0 / 32.0) * l[2] * tau[2][1] * thetaDir_2[2][1]
                          + (1.0 / 32.0) * l[2] * tau[2][2] * thetaDir_2[2][2])
                       * b_2[1][0];
          mat_[18][18] = b_0[0][0] * edgeOrientation[0];
          mat_[19][19] = b_1[0][0] * edgeOrientation[1];
          mat_[20][20] = b_2[0][0] * edgeOrientation[2];
        }
        BCRSMatrix<R> mat_;

      public:
        /**
         * \brief Class that evaluates the push forwards of the global nodes of a LocalFunction.
         *        It stretches the LocalInterpolation interface, because we evaluate the derivatives
         * of f
         * \tparam LocalBasis Basis of the reference element, to get the Domain and RangeFieldType
         * \tparam Element
         */
        template <class LocalBasis, class Element>
        class GlobalValuedInterpolation
        {
          using size_type = std::size_t;
          using LocalCoordinate = typename LocalBasis::Traits::DomainType;
          using RFT = typename LocalBasis::Traits::RangeFieldType;
          static constexpr unsigned int numberOfEdges = 3;

        public:
          GlobalValuedInterpolation() {}

          /**
           * \brief Binding routine. Collects tangentials and normals
           * \param element
           * \param elementInfo
           */
          void bind(const Element &element, const ElementInformation<Element> &elementInfo)
          {
            elementInfo_ = &elementInfo;
            auto geometry = element.geometry();

            // get global Normals and midpoints
            auto refElement = Dune::referenceElement<double, 2>(geometry.type());
            for (std::size_t i = 0; i < 3; ++i)
            {
              localVertices_[i] = refElement.position(i, 2);

              localMidpoints_[i] = refElement.position(i, 1);
              std::size_t lower = (i == 2) ? 1 : 0;
              std::size_t upper = (i == 0) ? 1 : 2;

              auto edge = geometry.global(refElement.position(upper, 2))
                        - geometry.global(refElement.position(lower, 2));
              // normalize and orient
              edge /= edge.two_norm() * elementInfo_->getEdgeOrientation()[i];
              // Rotation by pi/2. Note that Kirby rotates by 3*pi/2
              globalNormals_[i] = {-edge[1], edge[0]};
            }
          }

          /** \brief Evaluate a given function at the Lagrange nodes
           *
           * \tparam F Type of function to evaluate
           * \tparam C Type used for the values of the function
           * \param[in] ff Function to evaluate
           * \param[out] out Array of function values
           */
          template <typename F, typename C>
          void interpolate(const F &ff, std::vector<C> &out) const
          {
            // constexpr auto dim = LocalBasis::Traits::dimDomain;

            auto &&f
                = Dune::Impl::makeFunctionWithCallOperator<typename LocalBasis::Traits::DomainType>(
                    ff);

            auto df = derivative(ff);
            auto hessf = derivative(df);

            out.resize(LocalBasis::size());
            FieldVector<RFT, 2> derivativeValue;
            // Iterate over vertices, 6 dofs per vertex
            for (unsigned int i = 0; i < 3; ++i)
            {
              // matrix storing directions as rows (!)
              auto const &directionMatrix = elementInfo_->getDerivativeDirections()[i];

              directionMatrix.mv(matrixToVector(df(localVertices_[i])), derivativeValue);
              auto const &hessianValue = directionMatrix * tensorToMatrix(hessf(localVertices_[i]))
                                       * transpose(directionMatrix); //

              out[i * 6 + 0] = f(localVertices_[i]);
              out[i * 6 + 1] = matrixToVector(derivativeValue)[0];
              out[i * 6 + 2] = matrixToVector(derivativeValue)[1];
              out[i * 6 + 3]
                  = matrixToVector(hessianValue[0])[0]; // matrixToVector probably unneccesary here
              out[i * 6 + 4] = matrixToVector(hessianValue[0])[1];
              out[i * 6 + 5] = matrixToVector(hessianValue[1])[1];
            }

            // iterate over edges, one dof per edge
            for (unsigned int i = 0; i < numberOfEdges; ++i)
            {
              out[18 + i] = normalDerivative(ff, i);
            }
          }

        protected:
          ElementInformation<Element> const *elementInfo_;
          std::array<Dune::FieldVector<R, 2>, 3> globalNormals_;
          std::array<Dune::FieldVector<R, 2>, 3> localMidpoints_;
          std::array<Dune::FieldVector<R, 2>, 3> localVertices_;

          // Infrastructure for normal Derivative that allows evaluation of default oriented global
          // normal derivative, if f has this method
          template <class F>
          auto normalDerivative(F const &f, size_type i) const
          {
            return normalDerivativeImpl(f, i, PriorityTag<42>{});
          }

          template <class F,
                    decltype((std::declval<F>().normalDerivative(std::declval<size_type>()), true))
                    = true>
          auto normalDerivativeImpl(F const &f, size_type i, PriorityTag<4>) const
          {
            return f.normalDerivative(i);
          }

          template <class F, decltype((derivative(std::declval<F>()), true)) = true>
          auto normalDerivativeImpl(F const &f, size_type i, PriorityTag<3>) const
          {
            return matrixToVector(derivative(f)(localMidpoints_[i])).dot(globalNormals_[i]);
          }

          template <class F>
          auto normalDerivativeImpl(F const &f, size_type i, PriorityTag<1>) const
          {
            DUNE_THROW(Dune::NotImplemented,
                       Dune::className(f)
                           + " supports neither derivative(f) nor f.normalDerivative(i)!");
            return 0;
          }
        };
      };

      /**
       * \brief Class that creates a Mapping from Elements of a GridView to their ElementInformation
       *
       * \tparam GV GridView
       * \tparam R Rangetype of the finite element
       */
      template <typename GV, typename R>
      class ArgyrisElementInformationMap
      {
        using D = typename GV::ctype;
        using Element = typename GV::template Codim<0>::Entity;
        using GlobalCoordinate = typename Element::Geometry::GlobalCoordinate;
        using LocalCoordinate = typename Element::Geometry::LocalCoordinate;

        using IndexSet = typename GV::IndexSet;
        using IndexType = typename IndexSet::IndexType;
        static constexpr unsigned int dim = 2;
        static_assert(GV::dimension == dim);
        using ElementInformation =
            typename ArgyrisTransformator<R>::template ElementInformation<Element>;
        using VertexDataVector
            = std::vector<std::tuple<FieldMatrix<D, dim, dim>, std::bitset<2 * dim>,
                                     std::array<std::vector<IndexType>, dim>>>;

      public:
        ArgyrisElementInformationMap(GV const &gv, bool useTangentials = true)
            : elementMapper_(gv, mcmgElementLayout()), elementInformation_(gv.size(0)),
              useTangentials_(useTangentials)
        {
          fill(gv);
        }

        template <class Function>
        ArgyrisElementInformationMap(GV const &gv, Function const &f)
            : elementMapper_(gv, mcmgElementLayout()), elementInformation_(gv.size(0)),
              tangentialMap_(f), useTangentials_(true)
        {
          fill(gv);
        }

        void update(GV const &gv)
        {
          elementInformation_.resize(gv.size(0));
          elementMapper_.update(gv);
          fill(gv);
        }

        template <class Element>
        const auto &find(const Element &element) const
        {
          return elementInformation_[elementMapper_.index(element)];
        }

      private:
        class VertexDataHandle:
            public CommDataHandleIF<VertexDataHandle,
                                    std::tuple<FieldMatrix<D, dim, dim>, bool, bool>>
        {
          using VertexData = std::tuple<FieldMatrix<D, dim, dim>, bool, bool>;

        public:
          VertexDataHandle(IndexSet const &indexSet, VertexDataVector &data)
              : indexSet_(indexSet), vertexData_(data)
          {}

          bool contains(int dim, int codim) const { return (codim == dim); }

          template <class Entity>
          std::size_t size(Entity &entity) const
          {
            return 1;
          }

          bool fixedSize(int dim, int codim) const { return true; }

          template <class MessageBuffer, class Entity>
          void gather(MessageBuffer &messageBuffer, const Entity &entity) const
          {
            assert(Entity::mydimension == 0);
            auto ii = indexSet_.index(entity);

            auto const &[dir, bools, indices] = vertexData_[ii];

            messageBuffer.write(VertexData{dir, bools[2], bools[3]});
          }

          template <class MessageBuffer, class Entity>
          void scatter(MessageBuffer &messageBuffer, const Entity &entity, size_t n)
          {
            assert(n == 1);
            VertexData tmp;
            messageBuffer.read(tmp);
            IndexType const &index = indexSet_.index(entity);
            mergeData(tmp, vertexData_[index], index);
          }

        private:
          void mergeData(VertexData const &newData, typename VertexDataVector::value_type &oldData,
                         IndexType const &index)
          {
            auto const &[newDir, bool_1, bool_2] = newData;

            if (bool_1)
              setVertexData<false>(oldData, newDir[0], index);
            if (bool_2)
              setVertexData<false>(oldData, newDir[1], index);
          }
          IndexSet const &indexSet_;
          VertexDataVector &vertexData_;
        };

        void fill(const GV &gv)
        {
          bool sequentialSetup = (gv.comm().size() == 1);

          // compute orientation for all elements
          const auto &indexSet = gv.indexSet();
          unsigned short orientation;
          // vector with directions and bools to indicate whether this dir was already set
          VertexDataVector directionPerVertex(indexSet.size(dim),
            {{{1., 0.}, {0., 1.}}, (unsigned long) 0, {}});
          // interate over intersections
          if (useTangentials_)
          {
            for (const auto &element : elements(gv))
            {
              if (element.hasBoundaryIntersections())
              {
              auto elementIndex = elementMapper_.index(element);
              const auto &refElement = referenceElement(element);
              for (const auto &intersection : intersections(gv, element))
              {
                // fill vertex Data with linear independent tangentials
                if (intersection.boundary())
                {
                  auto vertexIterator = refElement.subEntities(intersection.indexInInside(), 1, 2).begin();
                  auto startIndex = indexSet.subIndex(element, *vertexIterator, 2);
                  auto endIndex = indexSet.subIndex(element, *++vertexIterator, 2);

                  auto tangential
                    = intersection.geometry().corner(1) - intersection.geometry().corner(0);
                  setVertexData(directionPerVertex[startIndex], tangential, elementIndex);
                  setVertexData(directionPerVertex[endIndex], tangential, elementIndex);
                }
                // inner dofs, keep global axes
              }
            }
            }

            if (!sequentialSetup)
            {
            // communicate vertex information Information
            VertexDataHandle dataHandle(indexSet, directionPerVertex);
            gv.communicate(dataHandle, All_All_Interface, ForwardCommunication);
            }

            // iterate over vector and set the remaining directions to normals
            for (auto &[dir, setTangential, indices] : directionPerVertex)
            {
              if (setTangential[2]) // if first is tangential
              {
                if (setTangential[3])
                // both are tangetials, i.e. vertex is corner
                {
                  // check if orthogonal
                  if (dir[0].dot(dir[1]) < 1e-14)
                  {
                    // direction is normal to all elements that other direction is tangential to
                    auto lenInd_1 = indices[1].size();
                    for (auto const &index : indices[0])
                      indices[1].push_back(index);
                    for (std::size_t i = 0; i < lenInd_1; ++i)
                      indices[0].push_back(indices[1][i]);
                  }
                }
                else // straight boundary segment
                {
                  // set second direction normal to first
                  dir[1][0] = -dir[0][1];
                  dir[1][1] = dir[0][0];
                  setTangential[1] = true;  // modified
                  setTangential[3] = false; // not tangential
                  // direction is normal to all elements that first direction is tangential to
                  for (auto const &index : indices[0])
                    indices[1].push_back(index);
                }
              } // change nothing if non was set
              // sanity check
              // TODO this is somewhat arbitrary number, more or less set
              // such that the unit tests trigger this exception rather than
              // failing the test
              // TODO handle this case
              if (!linearIndependent(dir[0], dir[1], 1e-3))
              {
                // Check wheter there are to similar tangentials
                DUNE_THROW(Dune::NotImplemented, "Almost Linear Dependent tangentials!");
              }
            }   // end for loop

          } // end if useTangentials_
            // note that there is no communication if tangentials are not used


          // compute orientation for element facets and fill map
          for (const auto &element : elements(gv))
          {
            const auto &refElement = referenceElement(element);
            auto elementIndex = elementMapper_.index(element);
            orientation = 0;
            std::array<Dune::FieldMatrix<D, dim, dim>, dim + 1> directions;
            // First loop to compute orientations
            for (std::size_t i = 0; i < element.subEntities(dim - 1); i++)
            {
              // Local vertex indices within the element
              auto localV0 = refElement.subEntity(i, dim - 1, 0, dim);
              auto localV1 = refElement.subEntity(i, dim - 1, 1, dim);

              if (sequentialSetup)
              {
                // Global vertex indices within the grid
                auto globalV0 = indexSet.subIndex(element, localV0, dim);
                auto globalV1 = indexSet.subIndex(element, localV1, dim);

                if ((localV0 < localV1 && globalV0 > globalV1)
                    || (localV0 > localV1 && globalV0 < globalV1))
                  orientation |= (1 << i);
              }
              else
              {
                // sort lexicographically by coordinate
                // this ensures consistent orientation also for distributed grids
                auto globalV0 = element.template subEntity<dim>(localV0).geometry().corner(0);
                auto globalV1 = element.template subEntity<dim>(localV1).geometry().corner(0);

                if ((localV0 < localV1 && vectorGreater(globalV0, globalV1))
                    || (localV0 > localV1 && vectorLess(globalV0, globalV1)))
                  orientation |= (1 << i);
              }
            }
            // Second loop to ensure vertices in directions are correctly
            // ordered
            std::array<std::bitset<3 * dim>, 3> booleans;

            for (std::size_t i = 0; i < element.subEntities(dim); ++i)
            {
              auto vertexIndex = indexSet.subIndex(element, i, dim);
              auto const &[dir, b, indices] = directionPerVertex[vertexIndex];
              directions[i] = dir;
              // modification and tangential information are taken directly
              for (std::size_t j = 0; j < b.size(); ++j)
                booleans[i][j] = b[j];
              // last two bits encode if the directions are tangential/normal to this very Element
              for (std::size_t j = 0; j < 2; ++j)
              {
                if (std::find(indices[j].begin(), indices[j].end(), elementIndex)
                    != indices[j].end())
                {
                  booleans[i][4 + j] = true;
                }
              }
            }

            // store information in map
            elementInformation_[elementIndex]
                = ElementInformation(orientation, directions, booleans);
          }
        }

        // properties

        Dune::MultipleCodimMultipleGeomTypeMapper<GV> elementMapper_;
        std::vector<ElementInformation> elementInformation_;
        // Map that assigns a boundary point a direction (usually tangential)
        // used to determine derivative directions at that point as tangetial and normal
        std::function<GlobalCoordinate(GlobalCoordinate)> tangentialMap_;
        bool useTangentials_ = true;

        // helper methods
        static bool linearIndependent(Dune::FieldVector<D, dim> a, Dune::FieldVector<D, dim> b,
                                      double tol = 1e-12)
        {
          return std::abs(a.dot(Dune::FieldVector<D, dim>{-b[1], b[0]})) > tol;
        }

        template <bool setIndices = true>
        static void
        setVertexData(std::tuple<FieldMatrix<D, dim, dim>, std::bitset<2 * dim>,
                                 std::array<std::vector<IndexType>, dim>> &directionPerVertex,
                      GlobalCoordinate const &t, IndexType const &index)
        {
          auto &[dir, setTangential, indices] = directionPerVertex;
          auto tangential = t / t.two_norm();
          // linear dependent to first, only set tangential bit and index
          if (!linearIndependent(dir[0], tangential))
          {
            setTangential[2] = true;
            if constexpr (setIndices)
              indices[0].push_back(index);
            return;
          }
          // linear dependent to second, only set tangential bit
          if (!linearIndependent(dir[1], tangential))
          {
            setTangential[3] = true;
            if constexpr (setIndices)
              indices[1].push_back(index);
            return;
          }
          // new tangential is not yet contained
            // first entry is not tangential, set first to new
            if (!setTangential[2])
            {
              dir[0] = tangential;
              setTangential[0] = true; // is modified
              setTangential[2] = true; // is tangential
              if constexpr (setIndices)
                indices[0].push_back(index);
              return;
            }
          // first is already tangential

              // second is not tangential
              if (!setTangential[3])
              {
                dir[1] = tangential;
                setTangential[1] = true; // is modified
                setTangential[3] = true; // is tangential
                if constexpr (setIndices)
                  indices[1].push_back(index);
              }
          // second is tangential too
              else
              {
                std::cout << "Old direction: " << dir[0] << "\t" << dir[1] << std::endl
                          << "New tangential: " << tangential << std::endl
                          << "Bools: " << setTangential << std::endl;
                DUNE_THROW(Dune::NotImplemented,
                           "More than two pairwise linear independent tangentials at vertex");
              }
        }
      };
    } // namespace Impl

    // forward declaration
    template <class GV, class R>
    class ArgyrisNode;
    /**
     * \brief PreBasis class for the Argyris Element
     *
     * \tparam GV GridView that we are bound to
     * \tparam R RangeFieldType of the finite Element
     */
    template <class GV, typename R>
    class ArgyrisPreBasis
    {
      static const int dim = GV::dimension;
      static_assert(dim == 2,
                    "Argyris PreBasis only implemented for 2d simplices");
      // using FiniteElementMap = Impl::ArgyrisLocalFiniteElementMap<GV, R>;
      using ElementInformationMap = Impl::ArgyrisElementInformationMap<GV, R>;

    public:
      using GridView = GV;
      using Range = R;
      using size_type = std::size_t;

      using Node = ArgyrisNode<GridView, Range>;

      static constexpr size_type maxMultiIndexSize = 1;
      static constexpr size_type minMultiIndexSize = 1;
      static constexpr size_type multiIndexBufferSize = 1;

      //! Constructor for a given gridView object
      ArgyrisPreBasis(const GridView &gv, bool useTangentials): gridView_(gv), elementInformationMap_(gv, useTangentials)
      {
        if (dim != 2)
          DUNE_THROW(Dune::NotImplemented, "ArgyrisPreBasis only implemented for dim == 2");
      }

      //! Initialize the global indices
      void initializeIndices() {}

      //! Obtain the grid view that the basis is defined on
      const GridView &gridView() const { return gridView_; }

      //! Update the stored grid view, to be called if the grid has changed
      void update(const GridView &gv)
      {
        gridView_ = gv;
        elementInformationMap_.update(gv);
      }

      /**
       * \brief Create tree node
       */

      Node makeNode() const { return Node{&elementInformationMap_}; }

      //! Same as size(prefix) with empty prefix
      size_type size() const
      {
        // 6 dofs per vertex, one dof per edge
        return 6 * gridView_.size(2) + gridView_.size(1);
      }

      //! Return number of possible values for next position in multi index
      template <class SizePrefix>
      size_type size(const SizePrefix prefix) const
      {
        assert(prefix.size() == 0 || prefix.size() == 1);
        return (prefix.size() == 0) ? size() : 0;
      }

      //! Get the total dimension of the space spanned by this basis
      size_type dimension() const { return size(); }

      //! Get the maximal number of DOFs associated to node for any element
      size_type maxNodeSize() const { return 21; }

      template <typename It>
      It indices(const Node &node, It it) const
      {
        const auto &gridIndexSet = gridView().indexSet();
        const auto &element = node.element();

        // throw if Element is not simplex
        if (not(element.type().isSimplex()))
          DUNE_THROW(Dune::NotImplemented, "Hermite Basis only implemented for simplex elements");
        for (size_type i = 0, end = node.finiteElement().size(); i < end; ++it, ++i)
        {
          Dune::LocalKey localKey = node.finiteElement().localCoefficients().localKey(i);

          *it = {{(size_type) ((localKey.codim() == dim)
                                  ? (6 * gridIndexSet.subIndex(element, localKey.subEntity(), dim)
                                     + localKey.index())
                                  : (6 * gridIndexSet.size(dim)
                                     + gridIndexSet.subIndex(element, localKey.subEntity(),
                                                             localKey.codim())))}};
        }
        return it;
      }

    protected:
      GridView gridView_;
      ElementInformationMap elementInformationMap_;
    };

    /**
     * \brief Node class for the Argyris element. This class holds a pointer to the
     * ElementInformation map, and add the respective ElementInformation object to the binding
     * routine.
     *
     * \tparam GV GridView that we are defined on
     * \tparam R RangeFieldType of the finite element
     */
    template <class GV, class R>
    class ArgyrisNode: public LeafBasisNode
    {
      static constexpr unsigned int dim = GV::dimension;
      using LocalValuedFE = ArgyrisLocalFiniteElement<typename GV::ctype, R>;
      using ElementInformationMap = Impl::ArgyrisElementInformationMap<GV, R>;

    public:
      using Element = typename GV::template Codim<0>::Entity;
      using size_type = std::size_t;
      using FiniteElement = Impl::LinearTransformedLocalFiniteElement<Impl::ArgyrisTransformator<R>, LocalValuedFE, Element>;

      ArgyrisNode(const ElementInformationMap *elementInformationMap)
          : localValuedFiniteElement(std::make_shared<LocalValuedFE>()), // default oriented!
            elementInformationMap_(elementInformationMap),
            finiteElement_(std::make_shared<FiniteElement>())
      {
        this->setSize(finiteElement_->size());
      }
      //! Return current element, throw if unbound
      const Element &element() const { return element_; }

      /** \brief Return the LocalFiniteElement for the element we are bound to
       *
       * The LocalFiniteElement implements the corresponding interfaces of the dune-localfunctions
       * module
       */
      const FiniteElement &finiteElement() const { return *finiteElement_; }

      //! Bind to element.
      void bind(const Element &e)
      {
        if (not e.type().isSimplex())
          DUNE_THROW(Dune::NotImplemented, "ArgyrisBasis can only be bound to simplex elements");
        element_ = e;

        finiteElement_->bind(*localValuedFiniteElement, element_,
                             elementInformationMap_->find(element_));
        this->setSize(finiteElement_->size());
      }

      unsigned int order() const { return 5; }

    private:
      std::shared_ptr<LocalValuedFE> localValuedFiniteElement;
      const ElementInformationMap *elementInformationMap_;
      std::shared_ptr<FiniteElement> finiteElement_;
      Element element_;
    };

    namespace BasisFactory{
      /**
       * \brief Create a pre-basis factory that can create Argyris pre-basis
       *
       * \ingroup FunctionSpaceBasesImplementations
       *
       * \tparam Range Numbertype used for shape function values
       *  \param useTangentials Whether to use the standart Argyris element with coordinate axis
       * oriented DOFs (false) or to orient the boundary DOFs tangential (and normal) to the
       * boundary
       *
       */
      template <typename Range = double>
      auto argyris(bool useTangentials = true)
      {
        return [=](auto const &gridView)
        { return ArgyrisPreBasis<std::decay_t<decltype(gridView)>, Range>(gridView, useTangentials); };
      }
    } // namespace BasisFactory
  }   // namespace Functions
} // namespace Dune

#endif