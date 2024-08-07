// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_FUNCTIONS_BASIX_FINITEELEMENT_HH
#define DUNE_FUNCTIONS_BASIX_FINITEELEMENT_HH

#if HAVE_BASIX

#include <optional>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/functions/functionspacebases/basix/localfiniteelement.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>

namespace Dune::Functions {

/**
 * \brief Definition of a global finite-element based on the basix library.
 *
 * \tparam G  The geometry of the global element the finite-element is defined on.
 * \tparam rangeClass  The class of the basis function range, i.e., RangeClass:scalar, vector or matrix.
 * \tparam F  A floating-point type used for the domain and range types of the basis functions.
 */
template <class G, RangeClass rangeClass, class F = double>
class BasixFiniteElement
{
public:
  using Geometry = G;
  using FieldType = F;
  using BasixType = ::basix::FiniteElement<FieldType>;
  static constexpr int dimDomain = Geometry::mydimension;
  using LocalFiniteElement = BasixLocalFiniteElement<dimDomain,rangeClass,FieldType>;

  struct Basis
  {
    static constexpr std::size_t dimRange ()
    {
      switch (rangeClass) {
        case RangeClass::scalar: return 1;
        case RangeClass::vector: return Geometry::coorddimension;
        case RangeClass::matrix: return Geometry::coorddimension * Geometry::coorddimension;
        default: return 0;
      }
    }

    using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;
    using Traits = LocalBasisTraits<
      F,dimDomain,FieldVector<F,dimDomain>,   // domain
      F,dimRange(),FieldVector<F,dimRange()>, // range
      FieldMatrix<F,dimRange(),dimDomain>     // jacobian
      >;

    using Domain = typename Traits::DomainType;
    using Range = typename Traits::RangeType;
    using Jacobian = typename Traits::JacobianType;

    /// \brief Return the number of basis functions
    std::size_t size () const
    {
      return lb_->size();
    }

    /// \brief Degree of the minimal polynomial space all basis functions are embedded
    std::size_t order () const
    {
      return lb_->order();
    }

    /// \brief Evaluate all shape functions in a point x
    void evaluateFunction (const Domain& x, std::vector<Range>& out) const
    {
      out.resize(size());

      using GlobalRange = Std::mdspan<F, Std::extents<std::size_t,Std::dynamic_extent,dimRange()>>;
#if USE_OUT_VECTOR_FOR_MDSPAN
      GlobalRange _global{&out[0][0], out.size()};
#else
      globalEvaluationBuffer_.resize(out.size() * dimRange());
      GlobalRange _global{globalEvaluationBuffer_.data(), out.size()};
#endif

      if constexpr (rangeClass == RangeClass::scalar)
      {
        if constexpr(std::is_same_v<Range, typename LocalBasis::Range>)
          lb_->evaluateFunction(x, _global);
        else
          DUNE_THROW(Dune::Exception, "Range type of global and local finite-element must match.");
      }
      else
      {
        localEvaluationBuffer_.resize(out.size() * LocalBasis::dimRange());
        using LocalRange = Std::mdspan<F, Std::extents<std::size_t,Std::dynamic_extent,LocalBasis::dimRange()>>;
        LocalRange _local{localEvaluationBuffer_.data(), out.size()};
        lb_->evaluateFunction(x, _local);

        auto J = geometry_->jacobian(x);
        auto K = geometry_->jacobianInverse(x);
        auto detJ = geometry_->integrationElement(x);

        using Jacobian = Std::mdspan<typename Geometry::ctype, Std::extents<std::size_t,Geometry::coorddimension,Geometry::mydimension>>;
        Jacobian _J{&J[0][0]};

        using JacobianInverse = Std::mdspan<typename Geometry::ctype, Std::extents<std::size_t,Geometry::mydimension,Geometry::coorddimension>>;
        JacobianInverse _K{&K[0][0]};

        auto map = lb_->basix().template map_fn<GlobalRange, LocalRange, Jacobian, JacobianInverse>();
        map(_global, _local, _J, detJ, _K);
      }

#if !USE_OUT_VECTOR_FOR_MDSPAN
      // copy the output values back into the output variable
      for (std::size_t i = 0; i < out.size(); ++i)
        for (std::size_t j = 0; j < dimRange(); ++j)
          out[i][j] = _global(i,j);
#endif
    }

    /// \brief Evaluate all shape function jacobians in a point x
    void evaluateJacobian (const Domain& x, std::vector<Jacobian>& out) const
    {
      out.resize(size());

      using GlobalJacobian = Std::mdspan<F, Std::extents<std::size_t,Std::dynamic_extent,dimRange(),dimDomain>>;
#if USE_OUT_VECTOR_FOR_MDSPAN
      GlobalJacobian _global{&out[0][0][0], out.size()};
#else
      globalEvaluationBuffer_.resize(out.size() * dimRange() * dimDomain);
      GlobalJacobian _global{globalEvaluationBuffer_.data(), out.size()};
#endif

      if constexpr (rangeClass == RangeClass::scalar)
      {
        if constexpr(std::is_same_v<Jacobian, typename LocalBasis::Jacobian>)
          lb_->evaluateJacobian(x, _global);
        else
          DUNE_THROW(Dune::Exception, "Jacobian type of global and local finite-element must match.");
      }
      else
      {
        localEvaluationBuffer_.resize(out.size() * LocalBasis::dimRange() * dimDomain);
        using LocalJacobian = Std::mdspan<F, Std::extents<std::size_t,Std::dynamic_extent,LocalBasis::dimRange(),dimDomain>>;
        LocalJacobian _local{localEvaluationBuffer_.data(), out.size()};
        lb_->evaluateJacobian(x, _local);

        // TODO: Transformation needs to be implemented.
        DUNE_THROW(Dune::NotImplemented, "Transform not yet implemented.");
      }

#if !USE_OUT_VECTOR_FOR_MDSPAN
      // copy the output values back into the output variable
      for (std::size_t i = 0; i < out.size(); ++i)
        for (std::size_t j = 0; j < dimRange(); ++j)
          for (std::size_t k = 0; k < dimDomain; ++k)
            out[i][j][k] = _global(i,j,k);
#endif
    }

    /// \brief Evaluate all shape function partial derivatives with given orders in a point x
    void partial (const std::array<unsigned int,dimDomain>& order,
                  const Domain& x, std::vector<Range>& out) const
    {
      out.resize(size());

      using GlobalPartials = Std::mdspan<F, Std::extents<std::size_t,Std::dynamic_extent,dimRange()>>;
#if USE_OUT_VECTOR_FOR_MDSPAN
      GlobalPartials _global{&out[0][0], out.size()};
#else
      globalEvaluationBuffer_.resize(out.size() * dimRange());
      GlobalPartials _global{globalEvaluationBuffer_.data(), out.size()};
#endif

      if constexpr (rangeClass == RangeClass::scalar)
      {
        if constexpr(std::is_same_v<Range, typename LocalBasis::Range>)
          lb_->partial(order, x, _global);
        else
          DUNE_THROW(Dune::Exception, "Range type of global and local finite-element must match.");
      }
      else
      {
        localEvaluationBuffer_.resize(out.size() * LocalBasis::dimRange());
        using LocalPartials = Std::mdspan<F, Std::extents<std::size_t,Std::dynamic_extent,LocalBasis::dimRange()>>;
        LocalPartials _local{localEvaluationBuffer_.data(), out.size()};
        lb_->partial(order, x, _local);

        // TODO: Transformation needs to be implemented.
        DUNE_THROW(Dune::NotImplemented, "Transform not yet implemented.");
      }

#if !USE_OUT_VECTOR_FOR_MDSPAN
      // copy the output values back into the output variable
      for (std::size_t i = 0; i < out.size(); ++i)
        for (std::size_t j = 0; j < dimRange(); ++j)
          out[i][j] = _global(i,j);
#endif
    }

    void bind (const Geometry& geometry)
    {
      geometry_ = &geometry;
    }

    const LocalBasis* lb_;
    const Geometry* geometry_ = nullptr;
    mutable std::vector<F> localEvaluationBuffer_ = {};
    mutable std::vector<F> globalEvaluationBuffer_ = {};
  };


  struct Coefficients
  {
    using LocalCoefficients = typename LocalFiniteElement::Traits::LocalCoefficientsType;

    /// \brief Return the number of local keys associated to local basis functions.
    std::size_t size () const
    {
      return lc_->size();
    }

    /// \brief Obtain the LocalKey associated to the `i`th basis function.
    const LocalKey& localKey (std::size_t i) const
    {
      return lc_->localKey(i);
    }

    void bind (const Geometry& geometry)
    {
      geometry_ = &geometry;
    }

    const LocalCoefficients* lc_;
    const Geometry* geometry_ = nullptr;
  };


  struct Interpolation
  {
    using LocalInterpolation = typename LocalFiniteElement::Traits::LocalInterpolationType;

    /// \brief Determine coefficients interpolating a given function `f`
    /// and store them in the output vector `out`.
    template<class Func, class C>
    void interpolate (const Func& f, std::vector<C>& out) const
    {
      li_->interpolate(f,out);
    }

    void bind (const Geometry& geometry)
    {
      geometry_ = &geometry;
    }

    const LocalInterpolation* li_;
    const Geometry* geometry_ = nullptr;
  };


  struct Traits {
    using Basis = typename BasixFiniteElement::Basis;
    using LocalBasisType = Basis;
    using Coefficients = typename BasixFiniteElement::Coefficients;
    using LocalCoefficientsType = Coefficients;
    using Interpolation = typename BasixFiniteElement::Interpolation;
    using LocalInterpolationType = Interpolation;
  };

public:
  /// \brief Construct a global finite-element from an associated local finite-element
  explicit BasixFiniteElement (LocalFiniteElement lfe)
    : lfe_(std::move(lfe))
    , basis_{&lfe_.localBasis()}
    , coefficients_{&lfe_.localCoefficients()}
    , interpolation_{&lfe_.localInterpolation()}
  {}

  /// \brief Construct the local finite-element from the basix library.
  explicit BasixFiniteElement (BasixType basix)
    : BasixFiniteElement(LocalFiniteElement(std::move(basix)))
  {}

  /// \brief Move constructor, needs to re-assign the internal pointers.
  BasixFiniteElement (const BasixFiniteElement& other)
    : BasixFiniteElement(other.lfe_)
  {
    if (other.geometry_.has_value())
      bind(*other.geometry_, other.cellInfo_);
  }

  /// \brief Move constructor, needs to re-assign the internal pointers.
  BasixFiniteElement (BasixFiniteElement&& other)
    : BasixFiniteElement(std::move(other.lfe_))
  {
    if (other.geometry_.has_value())
      bind(*other.geometry_, other.cellInfo_);
  }

  /// \brief Bind the global finite-element to a global element geometry and
  /// a cell permutation information.
  void bind (const Geometry& geometry, std::uint32_t cellInfo = 0)
  {
    geometry_.emplace(geometry);
    cellInfo_ = cellInfo;

    basis_.bind(*geometry_);
    coefficients_.bind(*geometry_);
    interpolation_.bind(*geometry_);
  }

  /// \brief Obtain a reference to the basis.
  const Basis& localBasis () const { return basis_; }
  const Basis& basis () const { return basis_; }

  /// \brief Obtain a reference to the coefficients.
  const Coefficients& localCoefficients () const { return coefficients_; }
  const Coefficients& coefficients () const { return coefficients_; }

  /// \brief Obtain a reference to the interpolation.
  const Interpolation& localInterpolation () const { return  interpolation_; }
  const Interpolation& interpolation () const { return  interpolation_; }

  /// \brief Return the dimension of the finite-element
  std::size_t size () const
  {
    return lfe_.size();
  }

  /// \brief Return the GeometryType the local finite-element is defined on
  GeometryType type () const
  {
    assert(geometry_.has_value());
    assert(lfe_.type() == geometry_->type());
    return lfe_.type();
  }

  std::uint32_t cellInfo () const
  {
    return cellInfo_;
  }

  /// \brief Obtain a reference to the basix implementation
  const BasixType& basix () const
  {
    return lfe_.basix();
  }

private:
  LocalFiniteElement lfe_;
  std::optional<Geometry> geometry_ = std::nullopt;
  std::uint32_t cellInfo_ = 0;

  Basis basis_;
  Coefficients coefficients_;
  Interpolation interpolation_;
};

} // end namespace Dune::Functions

#endif // HAVE_BASIX
#endif // DUNE_FUNCTIONS_BASIX_FINITEELEMENT_HH
