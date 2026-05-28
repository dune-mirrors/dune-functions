// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PULLBACKTRANSFORMEDLOCALFINITEELEMENT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PULLBACKTRANSFORMEDLOCALFINITEELEMENT_HH

#include <array>
#include <cassert>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

namespace Dune::Functions {

/**
 * \brief Tag types for selecting values or derivatives of local shape functions.
 *
 * These tags form the experimental dispatch interface used by transformed
 * local finite elements.  A tag object is passed to evaluate(), precompute(),
 * or finalize() to select the quantity to compute without overloading on the
 * output container type.
 *
 * The tags are intentionally lightweight value types.  Tags without data denote
 * a whole derivative quantity, while Partial stores the selected coordinate
 * direction.
 */
namespace Derivatives
{
  //! Select evaluation of shape-function values.
  struct Value {};

  //! Select evaluation of the Jacobian/first derivative of shape functions.
  struct Jacobian {};

  //! Select evaluation of gradients.
  struct Gradient {};

  //! Select evaluation of divergences.
  struct Divergence {};

  //! Select evaluation of Hessians.
  struct Hessian {};

  //! Select evaluation of Laplacians.
  struct Laplacian {};

  //! Select evaluation of a partial derivative in coordinate direction i.
  struct Partial { int i; };

} // end namespace Derivatives

namespace Impl {

template <class LocalBasis, int domainDimension>
using HessianMatrix = FieldMatrix<typename LocalBasis::Traits::RangeFieldType,domainDimension,domainDimension>;

template <class LocalBasis, int domainDimension>
using HessianTensor = std::array<HessianMatrix<LocalBasis,domainDimension>,LocalBasis::Traits::dimRange>;

template <class LocalBasis, class Geometry, class Derivative>
struct StandardDerivativeRange;

template <class LocalBasis, class Geometry>
struct StandardDerivativeRange<LocalBasis,Geometry,Derivatives::Value>
{
  using type = typename LocalBasis::Traits::RangeType;
};

template <class LocalBasis, class Geometry>
struct StandardDerivativeRange<LocalBasis,Geometry,Derivatives::Jacobian>
{
  using K = typename LocalBasis::Traits::RangeFieldType;
  using type = FieldMatrix<K,LocalBasis::Traits::dimRange,Geometry::coorddimension>;
};

template <class LocalBasis, class Geometry>
struct StandardDerivativeRange<LocalBasis,Geometry,Derivatives::Partial>
{
  using type = typename LocalBasis::Traits::RangeType;
};

template <class LocalBasis, class Geometry>
struct StandardDerivativeRange<LocalBasis,Geometry,Derivatives::Hessian>
{
  using type = HessianTensor<LocalBasis,Geometry::coorddimension>;
};

template <class LocalBasis, class Geometry>
struct StandardDerivativeRange<LocalBasis,Geometry,Derivatives::Laplacian>
{
  using type = typename LocalBasis::Traits::RangeType;
};

template <class LocalBasis, class Geometry, class Derivative>
struct ScalarDerivativeRange;

template <class LocalBasis, class Geometry>
struct ScalarDerivativeRange<LocalBasis,Geometry,Derivatives::Value>
{
  using type = typename LocalBasis::Traits::RangeType::value_type;
};

template <class LocalBasis, class Geometry>
struct ScalarDerivativeRange<LocalBasis,Geometry,Derivatives::Jacobian>
{
  using K = typename LocalBasis::Traits::RangeFieldType;
  using type = FieldVector<K,Geometry::coorddimension>;
};

template <class LocalBasis, class Geometry>
struct ScalarDerivativeRange<LocalBasis,Geometry,Derivatives::Partial>
{
  using type = typename LocalBasis::Traits::RangeType::value_type;
};

template <class LocalBasis, class Geometry>
struct ScalarDerivativeRange<LocalBasis,Geometry,Derivatives::Hessian>
{
  using type = HessianMatrix<LocalBasis,Geometry::coorddimension>;
};

template <class LocalBasis, class Geometry>
struct ScalarDerivativeRange<LocalBasis,Geometry,Derivatives::Laplacian>
{
  using type = typename LocalBasis::Traits::RangeType::value_type;
};

template <class LocalBasis, class Derivative>
struct PullbackPrecomputeBuffer;

template <class LocalBasis>
struct PullbackPrecomputeBuffer<LocalBasis,Derivatives::Value>
{
  using type = std::vector<typename LocalBasis::Traits::RangeType>;
};

template <class LocalBasis>
struct PullbackPrecomputeBuffer<LocalBasis,Derivatives::Jacobian>
{
  using type = std::vector<typename LocalBasis::Traits::JacobianType>;
};

template <class LocalBasis>
struct PullbackPrecomputeBuffer<LocalBasis,Derivatives::Partial>
{
  using type = std::vector<typename LocalBasis::Traits::JacobianType>;
};

template <class LocalBasis>
  requires requires { typename LocalBasis::Traits::HessianType; }
struct PullbackPrecomputeBuffer<LocalBasis,Derivatives::Hessian>
{
  using type = std::vector<typename LocalBasis::Traits::HessianType>;
};

template <class LocalBasis>
  requires (!requires { typename LocalBasis::Traits::HessianType; })
struct PullbackPrecomputeBuffer<LocalBasis,Derivatives::Hessian>
{
  using type = std::vector<HessianTensor<LocalBasis,LocalBasis::Traits::dimDomain>>;
};

template <class LocalBasis>
struct PullbackPrecomputeBuffer<LocalBasis,Derivatives::Laplacian>
{
  using type = typename PullbackPrecomputeBuffer<LocalBasis,Derivatives::Hessian>::type;
};

} // end namespace Impl

/**
 * \brief Derivative range traits for vector-valued pullback transformed bases.
 *
 * This traits class maps a derivative tag to the range type produced by
 * PullbackTransformedLocalBasis after the geometry pullback has been applied.
 * The nested alias template Range<Derivative>::type is the customization point
 * used by PullbackTransformedLocalBasis.
 *
 * The default implementation supports Derivatives::Value,
 * Derivatives::Jacobian, Derivatives::Partial, Derivatives::Hessian, and
 * Derivatives::Laplacian.  Values, partial derivatives, and Laplacians use the
 * range type of the reference local basis.  Jacobians are mapped to matrices
 * with one row per range component and one column per physical coordinate
 * direction.  Hessians are represented by an array of matrices, one for each
 * range component.
 *
 * \tparam LocalBasis Reference local basis type.
 * \tparam Geometry Bound element geometry type.
 */
template <class LocalBasis, class Geometry>
struct StandardDerivativeTraits
{
  //! Range traits for the quantity selected by Derivative.
  template <class Derivative>
  using Range = Impl::StandardDerivativeRange<LocalBasis,Geometry,Derivative>;
};

/**
 * \brief Derivative range traits for scalar components of pullback transformed bases.
 *
 * This traits class is useful when a vector-valued reference basis is used as a
 * source for scalar-valued component shape functions.  Values and partial
 * derivatives are represented by the scalar value_type of the reference range,
 * Jacobians are represented by physical-coordinate gradients, Hessians are
 * represented by matrices in physical coordinates, and Laplacians are
 * represented by scalar values.
 *
 * The nested alias template Range<Derivative>::type is the customization point
 * used by PullbackTransformedLocalBasis.
 *
 * \tparam LocalBasis Reference local basis type.
 * \tparam Geometry Bound element geometry type.
 */
template <class LocalBasis, class Geometry>
struct ScalarDerivativeTraits
{
  //! Range traits for the quantity selected by Derivative.
  template <class Derivative>
  using Range = Impl::ScalarDerivativeRange<LocalBasis,Geometry,Derivative>;
};

/**
 * \brief Local basis wrapper applying the geometry part of a pullback transform.
 *
 * PullbackTransformedLocalBasis adapts a reference-element local basis to an
 * element geometry.  Values and derivatives are requested with tags from the
 * Derivatives namespace.  For Derivatives::Jacobian and Derivatives::Partial,
 * the wrapper evaluates the reference-element Jacobian and transforms it with
 * the inverse element Jacobian, producing derivatives with respect to physical
 * coordinates.  For Derivatives::Hessian and Derivatives::Laplacian, second
 * reference partial derivatives are assembled into Hessian matrices and
 * transformed with the inverse element Jacobian from both sides.
 *
 * The class also exposes a two-stage evaluation protocol:
 *
 * - precompute() computes reference-element data that depends on the reference
 *   position and the local basis, but not on the physical geometry transform.
 * - finalize() converts such precomputed data into the public derivative range,
 *   using geometry information if needed.
 *
 * This split is intended as an experimental caching hook for local finite
 * element transformations.  evaluate() is the convenience interface combining
 * both stages.
 *
 * The object is non-owning: bind(LocalBasis const&) and bind(Geometry const&)
 * store pointers to externally managed objects.  The bound objects must outlive
 * this wrapper.
 *
 * \tparam LocalBasis Reference-element local basis type.
 * \tparam Geometry Element geometry type used for the pullback.
 * \tparam DerivativeTraits Traits mapping derivative tags to output range
 *                          types.  It must provide
 *                          DerivativeTraits::Range<Tag>::type.
 */
template <class LocalBasis, class Geometry,
          class DerivativeTraits = StandardDerivativeTraits<LocalBasis,Geometry>>
class PullbackTransformedLocalBasis
{
  public:

    //! Buffer type used by precompute() for the quantity selected by D.
    template <class D>
    using PrecomputeBuffer = typename Impl::PullbackPrecomputeBuffer<LocalBasis,D>::type;

    //! Public output range type for the quantity selected by D.
    template <class D>
    using DerivativeRange = typename DerivativeTraits::template Range<D>::type;

    //! Reference-element coordinate type.
    using Domain = typename LocalBasis::Traits::DomainType;

    //! Output range type for Derivatives::Value.
    using Range = DerivativeRange<Derivatives::Value>;

  public:
    /**
     * \brief Return the number of shape functions.
     *
     * The current pullback implementation forwards to the bound reference local
     * basis.
     */
    std::size_t size () const
    {
      return localBasis_->size();
    }

    /**
     * \brief Bind the reference local basis.
     *
     * The basis is stored by pointer and must remain valid while this object is
     * used.
     */
    void bind (LocalBasis const& localBasis)
    {
      localBasis_ = &localBasis;
    }

    /**
     * \brief Bind the element geometry used by finalize() and evaluate().
     *
     * The geometry is stored by pointer and must remain valid while this object
     * is used.
     */
    void bind (Geometry const& geometry)
    {
      geometry_ = &geometry;
    }

    /**
     * \brief Precompute reference-element shape-function values.
     *
     * This stage only requires the bound local basis.  The output can be reused
     * with finalize(Derivatives::Value, ...) for compatible requests.
     */
    void precompute (Derivatives::Value,
                     Domain const& x,
                     PrecomputeBuffer<Derivatives::Value>& out) const
    {
      assert(!!localBasis_);
      localBasis_->evaluateFunction(x,out);
    }

    /**
     * \brief Precompute reference-element shape-function Jacobians.
     *
     * This stage only requires the bound local basis.  The output stores
     * derivatives with respect to reference coordinates and can be reused with
     * finalize(Derivatives::Jacobian, ...) or
     * finalize(Derivatives::Partial, ...) for a bound geometry.
     */
    void precompute (Derivatives::Jacobian,
                     Domain const& x,
                     PrecomputeBuffer<Derivatives::Jacobian>& out) const
    {
      assert(!!localBasis_);
      localBasis_->evaluateJacobian(x,out);
    }

    /**
     * \brief Precompute reference-element shape-function Jacobians for a partial derivative.
     *
     * Global partial derivatives are selected after the geometry pullback.
     * Hence this precompute stage is identical to the Jacobian precomputation
     * and stores derivatives with respect to reference coordinates.
     */
    void precompute (Derivatives::Partial,
                     Domain const& x,
                     PrecomputeBuffer<Derivatives::Partial>& out) const
    {
      assert(!!localBasis_);
      localBasis_->evaluateJacobian(x,out);
    }

    /**
     * \brief Precompute reference-element shape-function Hessians.
     *
     * If the bound local basis provides evaluateHessian() and
     * Traits::HessianType, this method is used directly.  Otherwise this
     * implementation falls back to the local basis partial() interface,
     * evaluates all second-order partial derivatives, and assembles symmetric
     * Hessian matrices for all range components.
     */
    void precompute (Derivatives::Hessian,
                     Domain const& x,
                     PrecomputeBuffer<Derivatives::Hessian>& out) const
    {
      assert(!!localBasis_);

      if constexpr (requires { localBasis_->evaluateHessian(x,out); }) {
        localBasis_->evaluateHessian(x,out);
      }
      else {
        out.resize(localBasis_->size());

        std::vector<typename LocalBasis::Traits::RangeType> partialValues;
        std::array<unsigned int,LocalBasis::Traits::dimDomain> order;

        for (auto i : Dune::range(LocalBasis::Traits::dimDomain)) {
          for (auto j : Dune::range(i,LocalBasis::Traits::dimDomain)) {
            order.fill(0);
            ++order[i];
            ++order[j];
            localBasis_->partial(order,x,partialValues);

            for (auto k : Dune::range(out.size())) {
              for (auto r : Dune::range(LocalBasis::Traits::dimRange)) {
                out[k][r][i][j] = partialValues[k][r];
                out[k][r][j][i] = partialValues[k][r];
              }
            }
          }
        }
      }
    }

    /**
     * \brief Precompute reference-element Hessians for Laplacians.
     *
     * The Laplacian is the trace of the physical-coordinate Hessian, so the
     * precompute stage is identical to Hessian precomputation.
     */
    void precompute (Derivatives::Laplacian,
                     Domain const& x,
                     PrecomputeBuffer<Derivatives::Laplacian>& out) const
    {
      precompute(Derivatives::Hessian{}, x, out);
    }

    /**
     * \brief Convert precomputed reference values to the public value range.
     *
     * For StandardDerivativeTraits this is typically a copy.  For
     * ScalarDerivativeTraits this extracts/converts to the scalar component
     * value type.
     */
    void finalize (Derivatives::Value d,
                   Domain const& x,
                   PrecomputeBuffer<Derivatives::Value> const& in,
                   std::vector<DerivativeRange<Derivatives::Value>>& out) const
    {
      out.resize(in.size());
      using T = DerivativeRange<Derivatives::Value>;
      for (auto i : Dune::range(in.size()))
        out[i] = T(in[i]);
    }

    /**
     * \brief Transform precomputed reference Jacobians to physical derivatives.
     *
     * This stage requires a bound geometry.  For vector-valued ranges the
     * reference Jacobian is multiplied with geometry.jacobianInverse(x).  For
     * scalar-valued component ranges the inverse transposed Jacobian is applied
     * to the single reference-gradient row.
     */
    void finalize (Derivatives::Jacobian d,
                   Domain const& x,
                   PrecomputeBuffer<Derivatives::Jacobian> const& in,
                   std::vector<DerivativeRange<Derivatives::Jacobian>>& out) const
    {
      assert(!!geometry_);
      out.resize(in.size());

      if constexpr (requires{out[0][0][0];}) {
        auto&& Jinv = geometry_->jacobianInverse(x);
        for (auto i : Dune::range(in.size()))
          out[i] = in[i] * Jinv;
      }
      else {
        auto&& JinvT = geometry_->jacobianInverseTransposed(x);
        for (auto i : Dune::range(in.size()))
          JinvT.mv(in[i][0], out[i]);
      }
    }

    /**
     * \brief Transform reference Jacobians and select one physical partial derivative.
     *
     * This stage requires a bound geometry.  The reference Jacobians are first
     * transformed to physical-coordinate derivatives.  Only the column selected
     * by d.i is written to the output.
     */
    void finalize (Derivatives::Partial d,
                   Domain const& x,
                   PrecomputeBuffer<Derivatives::Partial> const& in,
                   std::vector<DerivativeRange<Derivatives::Partial>>& out) const
    {
      assert(!!geometry_);
      assert(0 <= d.i && d.i < Geometry::coorddimension);

      out.resize(in.size());

      if constexpr (requires{out[0][0];}) {
        auto&& Jinv = geometry_->jacobianInverse(x);
        for (auto i : Dune::range(in.size())) {
          DerivativeRange<Derivatives::Jacobian> jac = in[i] * Jinv;
          for (auto j : Dune::range(LocalBasis::Traits::dimRange))
            out[i][j] = jac[j][d.i];
        }
      }
      else {
        auto&& JinvT = geometry_->jacobianInverseTransposed(x);
        FieldVector<typename LocalBasis::Traits::RangeFieldType,Geometry::coorddimension> gradient;
        for (auto i : Dune::range(in.size())) {
          JinvT.mv(in[i][0], gradient);
          out[i] = gradient[d.i];
        }
      }
    }

    /**
     * \brief Transform reference Hessians to physical-coordinate Hessians.
     *
     * This stage requires a bound geometry.  Each reference Hessian is
     * transformed as \f$ J^{-T} H J^{-1} \f$, where J^{-1} is the inverse
     * element Jacobian returned by geometry.jacobianInverse(x).  For
     * vector-valued ranges this transformation is applied to each range
     * component separately.
     *
     * \note This is the affine-geometry Hessian transformation.  Curved
     * geometries require additional terms involving second derivatives of the
     * geometry map.
     */
    void finalize (Derivatives::Hessian d,
                   Domain const& x,
                   PrecomputeBuffer<Derivatives::Hessian> const& in,
                   std::vector<DerivativeRange<Derivatives::Hessian>>& out) const
    {
      assert(!!geometry_);
      out.resize(in.size());

      auto&& Jinv = geometry_->jacobianInverse(x);
      auto&& JinvT = Jinv.transposed();

      if constexpr (requires{out[0][0][0][0];}) {
        if constexpr (requires{in[0][0][0][0];}) {
          for (auto i : Dune::range(in.size()))
            for (auto r : Dune::range(LocalBasis::Traits::dimRange))
              out[i][r] = JinvT * in[i][r] * Jinv;
        }
        else {
          static_assert(LocalBasis::Traits::dimRange == 1,
            "A scalar HessianType can only be used with scalar local basis ranges.");
          for (auto i : Dune::range(in.size()))
            out[i][0] = JinvT * in[i] * Jinv;
        }
      }
      else {
        if constexpr (requires{in[0][0][0][0];}) {
          for (auto i : Dune::range(in.size()))
            out[i] = JinvT * in[i][0] * Jinv;
        }
        else {
          for (auto i : Dune::range(in.size()))
            out[i] = JinvT * in[i] * Jinv;
        }
      }
    }

    /**
     * \brief Transform reference Hessians and take their traces.
     *
     * This stage requires a bound geometry.  It first computes the
     * physical-coordinate Hessian and then writes the sum of its diagonal
     * entries to the output.
     *
     * \note This uses the same affine-geometry Hessian transformation as
     * finalize(Derivatives::Hessian, ...).
     */
    void finalize (Derivatives::Laplacian d,
                   Domain const& x,
                   PrecomputeBuffer<Derivatives::Laplacian> const& in,
                   std::vector<DerivativeRange<Derivatives::Laplacian>>& out) const
    {
      assert(!!geometry_);
      out.resize(in.size());

      auto&& Jinv = geometry_->jacobianInverse(x);
      auto&& JinvT = Jinv.transposed();

      if constexpr (requires{out[0][0];}) {
        for (auto i : Dune::range(in.size())) {
          out[i] = DerivativeRange<Derivatives::Laplacian>{};
          if constexpr (requires{in[0][0][0][0];}) {
            for (auto r : Dune::range(LocalBasis::Traits::dimRange)) {
              auto hessian = JinvT * in[i][r] * Jinv;
              for (auto j : Dune::range(Geometry::coorddimension))
                out[i][r] += hessian[j][j];
            }
          }
          else {
            static_assert(LocalBasis::Traits::dimRange == 1,
              "A scalar HessianType can only be used with scalar local basis ranges.");
            auto hessian = JinvT * in[i] * Jinv;
            for (auto j : Dune::range(Geometry::coorddimension))
              out[i][0] += hessian[j][j];
          }
        }
      }
      else {
        for (auto i : Dune::range(in.size())) {
          auto hessian = [&] {
            if constexpr (requires{in[0][0][0][0];})
              return JinvT * in[i][0] * Jinv;
            else
              return JinvT * in[i] * Jinv;
          }();
          out[i] = 0;
          for (auto j : Dune::range(Geometry::coorddimension))
            out[i] += hessian[j][j];
        }
      }
    }

    /**
     * \brief Evaluate shape-function values in one step.
     *
     * This is equivalent to precompute(Derivatives::Value, ...) followed by
     * finalize(Derivatives::Value, ...) unless the precomputed value type
     * already matches the public range type.
     */
    void evaluate (Derivatives::Value d,
                   Domain const& x,
                   std::vector<Range>& out) const
    {
      if constexpr(std::is_same_v<Range,typename LocalBasis::Traits::RangeType>)
        precompute(d,x,out);
      else {
        precompute(d,x,valueBuffer_);
        finalize(d,x,valueBuffer_,out);
      }
    }

    /**
     * \brief Evaluate physical-coordinate Jacobians in one step.
     *
     * This combines reference Jacobian precomputation with the geometry
     * pullback transformation.
     */
    void evaluate (Derivatives::Jacobian d,
                   Domain const& x,
                   std::vector<DerivativeRange<Derivatives::Jacobian>>& out) const
    {
      precompute(d,x,jacobianBuffer_);
      finalize(d,x,jacobianBuffer_,out);
    }

    /**
     * \brief Evaluate one physical-coordinate partial derivative in one step.
     *
     * This combines reference Jacobian precomputation with the geometry
     * pullback transformation and writes only the component selected by d.i.
     */
    void evaluate (Derivatives::Partial d,
                   Domain const& x,
                   std::vector<DerivativeRange<Derivatives::Partial>>& out) const
    {
      precompute(d,x,jacobianBuffer_);
      finalize(d,x,jacobianBuffer_,out);
    }

    /**
     * \brief Evaluate physical-coordinate Hessians in one step.
     *
     * This combines reference Hessian precomputation from second-order partial
     * derivatives with the geometry pullback transformation.
     */
    void evaluate (Derivatives::Hessian d,
                   Domain const& x,
                   std::vector<DerivativeRange<Derivatives::Hessian>>& out) const
    {
      precompute(d,x,hessianBuffer_);
      finalize(d,x,hessianBuffer_,out);
    }

    /**
     * \brief Evaluate physical-coordinate Laplacians in one step.
     *
     * This combines reference Hessian precomputation from second-order partial
     * derivatives with the geometry pullback transformation and trace.
     */
    void evaluate (Derivatives::Laplacian d,
                   Domain const& x,
                   std::vector<DerivativeRange<Derivatives::Laplacian>>& out) const
    {
      precompute(d,x,hessianBuffer_);
      finalize(d,x,hessianBuffer_,out);
    }

    //! Return the polynomial order of the bound reference local basis.
    int order () const
    {
      return localBasis_->order();
    }

  private:
    LocalBasis const* localBasis_ = nullptr;
    Geometry const* geometry_ = nullptr;

    mutable PrecomputeBuffer<Derivatives::Value> valueBuffer_ = {};
    mutable PrecomputeBuffer<Derivatives::Jacobian> jacobianBuffer_ = {};
    mutable PrecomputeBuffer<Derivatives::Hessian> hessianBuffer_ = {};
};


/**
 * \brief Local finite element wrapper exposing a pullback transformed basis.
 *
 * PullbackTransformedLocalFiniteElement is a small non-owning adapter around a
 * local finite element and an element geometry.  It binds the local finite
 * element's local basis into a PullbackTransformedLocalBasis and records the
 * bound geometry type.  The object intentionally only models the local finite
 * element aspect needed for element-local computations.
 *
 * The wrapper currently exposes the transformed basis, size, and geometry type.
 * Interpolation and local coefficients are outside the current experimental
 * interface.
 *
 * \tparam LocalFiniteElement Reference local finite element type.
 * \tparam Geometry Element geometry type used for the pullback.
 * \tparam DerivativeTraits Traits mapping derivative tags to transformed basis
 *                          output range types.
 */
template <class LocalFiniteElement, class Geometry,
          class DerivativeTraits = StandardDerivativeTraits<typename LocalFiniteElement::Traits::LocalBasisType,Geometry>>
class PullbackTransformedLocalFiniteElement
{
    using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;

  public:
    //! Type of the transformed local basis exposed by basis().
    using Basis = PullbackTransformedLocalBasis<LocalBasis, Geometry, DerivativeTraits>;

    /**
     * \brief Bind the reference local finite element.
     *
     * The local finite element is not stored directly.  Its local basis is
     * stored by pointer inside the transformed basis and must remain valid while
     * this object is used.
     */
    void bind (LocalFiniteElement const& lfe)
    {
      basis_.bind(lfe.localBasis());
    }

    /**
     * \brief Bind the element geometry used by the transformed basis.
     *
     * The geometry is stored by pointer inside the transformed basis and must
     * remain valid while this object is used.
     */
    void bind (Geometry const& geometry)
    {
      basis_.bind(geometry);
      type_ = geometry.type();
    }

    //! Return the transformed local basis.
    Basis const& basis () const
    {
      return basis_;
    }

    //! Return the number of shape functions of the transformed basis.
    std::size_t size () const
    {
      return basis_.size();
    }

    //! Return the geometry type of the bound element geometry.
    GeometryType type () const
    {
      return type_;
    }

  protected:
    Basis basis_;
    GeometryType type_;
};

} // end namespace Dune::Functions

#endif
