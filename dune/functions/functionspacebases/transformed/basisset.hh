// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_BASISSET_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TRANSFORMED_BASISSET_HH

#include <utility>
#include <vector>

namespace Dune::Functions {

/**
 * \brief Pipeline stage that linearly mixes complete basis-function sets.
 *
 * The wrapped policy has to provide bind(context) and
 *
 * \code{.cpp}
 * template<class In, class Out>
 * void transformBasisSet(In const& in, Out& out) const;
 * \endcode
 */
template<class Transformation>
class BasisSetTransformationStage
{
  public:
    template<class Derivative, class LocalBasis, class Context, class InputRange>
    using OutputRange = InputRange;

    BasisSetTransformationStage() = default;

    explicit BasisSetTransformationStage(Transformation transformation)
      : transformation_(std::move(transformation))
    {}

    template<class Context>
    void bind(Context const& context)
    {
      transformation_.bind(context);
    }

    template<class Derivative, class LocalBasis, class InputRange, class OutputRange>
    void transform(Derivative,
                   LocalBasis const&,
                   typename LocalBasis::Traits::DomainType const&,
                   std::vector<InputRange> const& in,
                   std::vector<OutputRange>& out) const
    {
      transformation_.transformBasisSet(in,out);
    }

    Transformation const& transformation() const
    {
      return transformation_;
    }

  private:
    Transformation transformation_;
};

} // end namespace Dune::Functions

#endif
