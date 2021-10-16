// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LEAFBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LEAFBASIS_HH

#include <cassert>
#include <cstddef>

namespace Dune {
namespace Functions {

/**
 * \brief CRTP base class for leaf pre-bases
 *
 * \tparam Derived  The actual implementation of a leaf pre-basis
 */
template<class Derived>
class LeafPreBasis
{
public:
  //! Type used for indices and size information
  using size_type = std::size_t;

  static constexpr size_type maxMultiIndexSize = 1;
  static constexpr size_type minMultiIndexSize = 1;
  static constexpr size_type multiIndexBufferSize = 1;

  //! Return number of possible values for next position in multi index
  template<class SizePrefix>
  size_type size(const SizePrefix& prefix) const
  {
    assert(prefix.size() == 0 || prefix.size() == 1);
    return (prefix.size() == 0) ?  derived().dimension() : 0;
  }

  //! Get the total dimension of the space spanned by this basis
  size_type size() const
  {
    return derived().dimension();
  }

private:
  const Derived& derived() const
  {
    return static_cast<const Derived&>(*this);
  }

  Derived& derived()
  {
    return static_cast<Derived&>(*this);
  }
};


} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_LEAFBASIS_HH
