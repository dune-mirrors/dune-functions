// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PULLBACKTRANSFORMEDLOCALFINITEELEMENT_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PULLBACKTRANSFORMEDLOCALFINITEELEMENT_HH

#include <cstddef>

#include <dune/geometry/type.hh>

#include <dune/functions/functionspacebases/transformed/bindcontext.hh>
#include <dune/functions/functionspacebases/transformed/derivativetraits.hh>
#include <dune/functions/functionspacebases/transformed/localfiniteelement.hh>
#include <dune/functions/functionspacebases/transformed/pullback.hh>

namespace Dune::Functions {

/**
 * \brief Local basis wrapper applying the geometry part of a pullback transform.
 *
 * This is the compatibility facade for the pullback-transformed local-basis
 * interface.  The implementation is delegated to the generic
 * TransformedLocalBasis adapter and the DerivativePullback transformation
 * policy.
 *
 * The object is non-owning: bind(LocalBasis const&) and bind(Geometry const&)
 * store pointers to externally managed objects.  The bound objects must outlive
 * this wrapper.
 */
template <class LocalBasis, class Geometry,
          class DerivativeTraits = StandardDerivativeTraits<LocalBasis,Geometry>>
class PullbackTransformedLocalBasis
  : public TransformedLocalBasis<
      LocalBasis,
      GeometryBindContext<Geometry>,
      DerivativePullback<Geometry,FixedDerivativeTraitsFactory<DerivativeTraits>>>
{
    using Transformation = DerivativePullback<Geometry,FixedDerivativeTraitsFactory<DerivativeTraits>>;
    using Context = GeometryBindContext<Geometry>;
    using Base = TransformedLocalBasis<LocalBasis,Context,Transformation>;

  public:
    using Base::Base;
    using Base::bind;

    /**
     * \brief Bind the element geometry used by finalize() and evaluate().
     */
    void bind (Geometry const& geometry)
    {
      context_.bind(geometry);
      Base::bind(context_);
    }

  private:
    Context context_;
};


/**
 * \brief Local finite element wrapper exposing a pullback transformed basis.
 *
 * PullbackTransformedLocalFiniteElement is a small non-owning adapter around a
 * local finite element and an element geometry.  It binds the local finite
 * element's local basis into a PullbackTransformedLocalBasis and records the
 * bound geometry type.
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
     */
    void bind (LocalFiniteElement const& lfe)
    {
      basis_.bind(lfe.localBasis());
    }

    /**
     * \brief Bind the element geometry used by the transformed basis.
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
    GeometryType type_ = {};
};

} // end namespace Dune::Functions

#endif
