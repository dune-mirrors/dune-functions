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

} // end namespace Impl

/**
 * \brief Derivative range traits for vector-valued pullback transformed bases.
 *
 * This traits class maps a derivative tag to the range type produced by
 * PullbackTransformedLocalBasis after the geometry pullback has been applied.
 * The nested alias template Range<Derivative>::type is the customization point
 * used by PullbackTransformedLocalBasis.
 *
 * The default implementation supports Derivatives::Value and
 * Derivatives::Jacobian.  Values keep the range type of the reference local
 * basis.  Jacobians are mapped to matrices with one row per range component and
 * one column per physical coordinate direction.
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
 * source for scalar-valued component shape functions.  Values are represented
 * by the scalar value_type of the reference range, and Jacobians are represented
 * by physical-coordinate gradients.
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
 * Derivatives namespace.  For Derivatives::Jacobian, the wrapper evaluates the
 * reference-element Jacobian and transforms it with the inverse element
 * Jacobian, producing derivatives with respect to physical coordinates.
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
     * finalize(Derivatives::Jacobian, ...) for a bound geometry.
     */
    void precompute (Derivatives::Jacobian,
                     Domain const& x,
                     PrecomputeBuffer<Derivatives::Jacobian>& out) const
    {
      assert(!!localBasis_);
      localBasis_->evaluateJacobian(x,out);
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
