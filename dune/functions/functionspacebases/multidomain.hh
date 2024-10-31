// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_MULTIDOMAIN_HH
#define DUNE_FUNCTIONS_MULTIDOMAIN_HH

#include "compositebasis.hh"
#include "multidomain/restricted.hh"
#include "multidomain/domaininfo.hh"
#include "multidomain/restrictedgridview.hh"
#include "multidomain/multidomaingridview.hh"

namespace Dune::Functions::MultiDomain
{

// just make things more readable ... could we change this into a literal?!
int subdomain(int i) { return i; }

} // end namespace Dune::Functions::MultiDomain

namespace Dune::Functions::BasisFactory {

/**
 * \brief Create a pre-basis factory that can build a PowerPreBasis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam ChildPreBasisFactory Types of child pre-basis factory
 * \tparam IndexMergingStrategy An IndexMergingStrategy type
 * \param childPreBasisFactory Child pre-basis factory
 * \param ims IndexMergingStrategy to be used
 *
 * This overload can be used to explicitly supply an IndexMergingStrategy.
 */
// template<std::size_t k, class ChildPreBasisFactory, class IndexMergingStrategy>
// auto subdomain(ChildPreBasisFactory&& childPreBasisFactory, const IndexMergingStrategy&)
// {
//   return [childPreBasisFactory](const auto& gridView) {
//     auto childPreBasis = childPreBasisFactory(gridView);
//     return PowerPreBasis<IndexMergingStrategy, decltype(childPreBasis), k>(std::move(childPreBasis));
//   };
// }

/**
 * \brief Create a factory builder that can build a PowerPreBasis
 *
 * \ingroup FunctionSpaceBasesImplementations
 *
 * \tparam ChildPreBasisFactory Types of child pre-basis factory
 * \param childPreBasisFactory Child pre-basis factory
 *
 * This overload will select the BasisFactory::BlockedInterleaved strategy.
 */
template<class ChildPreBasisFactory, typename SubDomain>
auto restrict(ChildPreBasisFactory&& childPreBasisFactory, SubDomain&& subdomain)
{
  return [childPreBasisFactory, subdomain](const auto& gridView) {
#warning (LATER) static_assert(... MultiDomainGridView ...)

    // the gridView must be a MultiDomainGridView that simply wraps
    // a gridView and allows to get access to the individual
    // subdomainGridView's.
    const auto & subdomainGridView = gridView.subdomainGridView(subdomain);
    // instantiate child with restricted gridview
    auto childPreBasis = childPreBasisFactory(subdomainGridView);
    return RestrictedPreBasis<decltype(childPreBasis)>(std::move(childPreBasis),subdomain);
  };
}

#warning (LATER) support multiDomainPower

template<class... Args>
#warning (LATER) support index-merging-strategy
auto multiDomainComposite(
  const std::shared_ptr<MultiDomain::PartitionedDomainInfo> & domainInfo,
  Args&&... subDomainPreBasisFactories)
{
#warning (LATER) static_assert(... restricted nodes ...)
  auto compositeFactory = composite(std::forward<Args>(subDomainPreBasisFactories)...);
  return [compositeFactory,domainInfo](const auto& gridView) {
    using GridView = std::decay_t<decltype(gridView)>;
    using MultiDomainGridView = Dune::Functions::MultiDomain::MultiDomainGridView<GridView>;
    auto mdgv = std::make_shared<MultiDomainGridView>(gridView,domainInfo);
    auto compositePreBasis = compositeFactory(*mdgv);
    using CompositePreBasis = std::decay_t<decltype(compositePreBasis)>;
    return MultiDomainPreBasis<CompositePreBasis, MultiDomainGridView>(std::move(compositePreBasis), mdgv);
  };
}

} // end namespace BasisFactory

#endif // DUNE_FUNCTIONS_MULTIDOMAIN_HH
