// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MORLEYBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_MORLEYBASIS_HH

#include <dune/common/exceptions.hh>
#include <type_traits>

#include <dune/c1elements/tensormatvec.hpp>
#include <dune/functions/functionspacebases/lineartransformedlocalfiniteelement.hh>
#include <dune/functions/functionspacebases/nodes.hh>
#include <dune/functions/functionspacebases/vectorfloatcmp.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/localfunctions/morley.hh>
// TODO Rework this to incorporate some stuff for nonaffine mappings

namespace Dune
{
  namespace Functions
  {
    namespace Impl
    {
      /** \brief Linear transformation that maps the reference basis onto the pull-backs of
       * physical nodal basis for the Morley element
       * \tparam R RangeFieldType of finite element
       */
      template <class R>
      class MorleyTransformator
      {
      public:
        /** \brief class holding the orientation of normal derivatives
         * \tparam Element Element type for compability
         */
        template <class Element>
        class ElementInformation
        {
          using GlobalCoordinate = typename Element::Geometry::GlobalCoordinate;

        public:
          ElementInformation() : edgeOrientation_({1., 1., 1.}) {}
          ElementInformation(std::bitset<3> orientationBitset)
          {
            for (std::size_t i = 0; i < 3; ++i)
              edgeOrientation_[i] = orientationBitset[i] ? -1. : 1.;
          }
          std::array<R, 3> const &getEdgeOrientation() const { return edgeOrientation_; }
          /**
           * \brief return the dof corresponding to localKey is a Dirichlet
           * dof (i.e. is a tangential)
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
              if (localKey.index() == 0) // value ->Dirichlet
                return true;
            }
            DUNE_THROW(NotImplemented, "Invalid LocalKey");
          }

        private:
          std::array<R, 3> edgeOrientation_;
        };

        MorleyTransformator() : mat_(6, 6, BCRSMatrix<R>::random) { setupMatrix(); }

        /**
         * \brief binds the transformation to an element and its elementinformation
         *        Fills the transformation Matrix.
         *
         * \tparam Element
         * \param element
         * \param elementInformation
         */
        template <class Element>
        void bind(Element const &element, ElementInformation<Element> elementInformation)
        {
          fillMatrix(
              Dune::referenceElement<double, 2>(GeometryTypes::simplex(2)).position(0, 0),
              element.geometry(),
              elementInformation.getEdgeOrientation()); // barycenter, because we need some value.
        }

        /**
         * \brief Apply the transformation to some Vector of Shapevalues, Jacobians or Hessians
         *
         * \tparam Values Vector
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
        {
          for (std::size_t i = 0; i < 3; ++i)
          {
            mat_.setrowsize(i, 3);
            mat_.setrowsize(3 + i, 1);
          }
          mat_.endrowsizes();
          for (std::size_t i = 0; i < 3; ++i)
          {
            mat_.addindex(i, i);
            for (std::size_t j = 0; j < 3; ++j)
              if (j != (2 - i))
                mat_.addindex(j, 3 + i);
            mat_.addindex(3 + i, 3 + i);
          }
          mat_.endindices();
        }

        template <class LocalCoordinate, class Geometry>
        void fillMatrix(LocalCoordinate const &x, Geometry const &geometry, std::array<R,3> edgeOrientation)
        {
          std::array<R, 3> B_11;
          std::array<R, 3> B_12;
          std::array<R, 3> l_inv;

          std::array<Dune::FieldVector<R, 2>, 3> referenceTangents;
          std::array<Dune::FieldVector<R, 2>, 3> globalTangents;

          // By default, edges point from the vertex with the smaller index
          // to the vertex with the larger index. Note that the B matrix is invariant of
          // orientation, since the -1 s of the reference normals/tangents cancel the -1 s
          // from the global normals/tangents normalize

          // get local and global Tangents
          auto refElement = Dune::referenceElement<double, 2>(geometry.type());
          for (std::size_t i = 0; i < 3; ++i)
          {
            std::size_t lower = (i == 2) ? 1 : 0;
            std::size_t upper = (i == 0) ? 1 : 2;
            auto edge = refElement.position(upper, 2) - refElement.position(lower, 2);

            referenceTangents[i] = edge / edge.two_norm();

            auto globalEdge = geometry.global(refElement.position(upper, 2))
                            - geometry.global(refElement.position(lower, 2));

            l_inv[i] = 1. / globalEdge.two_norm();
            globalTangents[i] = globalEdge * l_inv[i];
          }

          auto jacobianTransposed = geometry.jacobianTransposed(x);
          for (std::size_t i = 0; i < 3; ++i)
          {
            B_11[i] = -referenceTangents[i][1]
                        * (-globalTangents[i][1] * jacobianTransposed[0][0]
                           + globalTangents[i][0] * jacobianTransposed[0][1])
                    + referenceTangents[i][0]
                          * (-globalTangents[i][1] * jacobianTransposed[1][0]
                             + globalTangents[i][0] * jacobianTransposed[1][1]);
            B_12[i] = -referenceTangents[i][1]
                        * (globalTangents[i][0] * jacobianTransposed[0][0]
                           + globalTangents[i][1] * jacobianTransposed[0][1])
                    + referenceTangents[i][0]
                          * (globalTangents[i][0] * jacobianTransposed[1][0]
                             + globalTangents[i][1] * jacobianTransposed[1][1]);
          }

          // Actually setup matrix
          int sign = -1;
          for (std::size_t i = 0; i < 3; ++i)
          {
            mat_[i][i] = 1.;
            for (std::size_t j = 0; j < 3; ++j)
              if (j != (2 - i)) // dune specific edge order
              {
                mat_[j][3 + i] = sign * B_12[i] * l_inv[i];
                sign *= -1;
              }
            mat_[3 + i][3 + i] = edgeOrientation[i] * B_11[i];
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
        class GlobalValuedInterpolation
        {
          using size_type = std::size_t;
          using LocalCoordinate = typename LocalBasis::Traits::DomainType;
          static constexpr unsigned int numberOfEdges = 3;

        public:
          GlobalValuedInterpolation() {}

          void bind(const Element &element, const ElementInformation<Element> &elementInfo)
          {
            elementInfo_ = elementInfo;
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
              edge /= edge.two_norm() * elementInfo_.getEdgeOrientation()[i];
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
            // using D = typename LocalBasis::Traits::DomainFieldType;

            auto &&f = Dune::Impl::makeFunctionWithCallOperator<LocalCoordinate>(ff);

            out.resize(LocalBasis::size());

            for (size_type i = 0; i < 3; ++i)
            {
              out[i] = f(localVertices_[i]);
              out[3 + i] = normalDerivative(ff, i);
            }
          }

        protected:
          ElementInformation<Element> elementInfo_;
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

      };   // Transformator

      template <typename GV, typename R>
      class MorleyElementInformationMap
      {
        using D = typename GV::ctype;
        using Element = typename GV::template Codim<0>::Entity;
        static constexpr unsigned int dim = 2;
        static_assert(GV::dimension == dim);
        using ElementInformation =
            typename MorleyTransformator<R>::template ElementInformation<Element>;

      public:
        MorleyElementInformationMap(GV const &gv)
            : elementMapper_(gv, mcmgElementLayout()), elementInformation_(gv.size(0))
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
        void fill(const GV &gv)
        {

          // compute orientation for all elements
          unsigned short orientation;
          bool sequentialSetup = (gv.comm().size() == 1);
          auto const &indexSet = gv.indexSet();
          for (const auto &element : elements(gv))
          {
            const auto &refElement = referenceElement(element);
            auto elementIndex = elementMapper_.index(element);

            orientation = 0;

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
            elementInformation_[elementIndex] = ElementInformation(orientation);
          }
        }
        Dune::MultipleCodimMultipleGeomTypeMapper<GV> elementMapper_;
        std::vector<ElementInformation> elementInformation_;
      };

    } // namespace Impl

    template <class GV, class R>
    class MorleyNode;

    template <class GV, typename R>
    class MorleyPreBasis
    {
      static const int dim = GV::dimension;
      static_assert(dim == 2, "Morley PreBasis only implemented for 2d simplices");
      using ElementInformationMap = Impl::MorleyElementInformationMap<GV, R>;

    public:
      //! The grid view that the FE basis is defined on
      using GridView = GV;

      //! Type used for indices and size information
      using size_type = std::size_t;

      //! Template mapping root tree path to type of created tree node
      using Node = MorleyNode<GridView, R>;

      static constexpr size_type maxMultiIndexSize = 1;
      static constexpr size_type minMultiIndexSize = 1;
      static constexpr size_type multiIndexBufferSize = 1;

      //! Constructor for a given grid view object
      MorleyPreBasis(const GV &gv) : gridView_(gv), elementInformationMap_(gv) {}

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
      size_type size() const { return gridView_.size(2) + gridView_.size(1); }

      //! Return number of possible values for next position in multi index
      template <class SizePrefix>
      size_type size(const SizePrefix prefix) const
      {
        // this basically means this is a leaf node with dimrange 1 right?
        assert(prefix.size() == 0 || prefix.size() == 1);
        return (prefix.size() == 0) ? size() : 0;
      }

      //! Get the total dimension of the space spanned by this basis
      size_type dimension() const { return size(); }

      //! Get the maximal number of DOFs associated to node for any element
      size_type maxNodeSize() const { return 6; }

      template <typename It>
      It indices(const Node &node, It it) const
      {
        const auto &gridIndexSet = gridView().indexSet();
        const auto &element = node.element();

        // throw if Element is not simplex
        if (not(element.type().isSimplex()))
          DUNE_THROW(Dune::NotImplemented, "Morley Basis only implemented for simplex elements");
        for (size_type i = 0, end = node.finiteElement().size(); i < end; ++it, ++i)
        {
          Dune::LocalKey localKey = node.finiteElement().localCoefficients().localKey(i);
          // TODO probably unnecessary
          if (!gridIndexSet.contains(element))
            DUNE_THROW(Dune::RangeError, "Element is not in gridIndexSet!");

          *it = {
              {(size_type)((localKey.codim() == dim)
                               ? gridIndexSet.subIndex(element, localKey.subEntity(), dim)
                               : gridIndexSet.size(dim)
                                     + gridIndexSet.subIndex(element, localKey.subEntity(), 1))}};
        }
        return it;
      }

    protected:
      GridView gridView_;
      ElementInformationMap elementInformationMap_;
    };

    template <class GV, class R>
    class MorleyNode: public LeafBasisNode
    {
      static constexpr unsigned int dim = GV::dimension;
      using LocalValuedFE = MorleyLocalFiniteElement<typename GV::ctype, R>;
      using ElementInformationMap = Impl::MorleyElementInformationMap<GV, R>;

    public:
      using size_type = std::size_t;
      using Element = typename GV::template Codim<0>::Entity;

      using FiniteElement = Impl::LinearTransformedLocalFiniteElement<Impl::MorleyTransformator<R>,
                                                                      LocalValuedFE, Element>;

      MorleyNode(const ElementInformationMap *elementInfoMap)
          : localValuedFiniteElement_(std::make_shared<LocalValuedFE>()),
            finiteElement_(std::make_shared<FiniteElement>()),
            elementInformationMap_(elementInfoMap)
      {
        this->setSize(finiteElement_->size());
      }

      ~MorleyNode() {}

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
          DUNE_THROW(Dune::NotImplemented, "MorleyBasis can only be bound to simplex elements");
        element_ = e;

        finiteElement_->bind(*localValuedFiniteElement_, element_,
                             elementInformationMap_->find(element_));
      }

      unsigned int order() const { return 2; }

    protected:
      std::shared_ptr<LocalValuedFE> localValuedFiniteElement_;
      std::shared_ptr<FiniteElement> finiteElement_;
      const ElementInformationMap *elementInformationMap_;
      Element element_;
    };

    namespace BasisFactory
    {
      /**
       * \brief Create a pre-basis factory that can create Morley pre-basis
       *
       * \ingroup FunctionSpaceBasesImplementations
       *
       * \tparam Range Numbertype used for shape function values
       *
       *
       */

      template <typename Range = double>
      auto morley()
      {
        return [](auto const &gridView)
        { return MorleyPreBasis<std::decay_t<decltype(gridView)>, Range>(gridView); };
      }
    } // namespace BasisFactory

  } // namespace Functions
} // namespace Dune

#endif
