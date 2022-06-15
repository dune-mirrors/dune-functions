#ifndef DUNE_FUNCTIONS_COMMON_VECTORSPAN
#define DUNE_FUNCTIONS_COMMON_VECTORSPAN

namespace Dune::Functions::Imp {

/**
 * \brief provide access to a continous span of an existing vector
 *
 * This is intended to be used to compute `jit.mv(x, y)` where `jit`
 * is the transposed inverse Jacobian of a grid element's geometry.
 *
 * \tparam V underlying vector type
 * \tparam n size of span
 */
template<class V, std::size_t n>
class VectorSpan
{
  V& v_;
  const std::size_t offset_;

public:
  constexpr VectorSpan(V& v, std::size_t offset)
    : v_(v)
    , offset_(offset)
  {
    /* Nothing. */
  }

  decltype(auto) operator[](std::size_t i)
  {
    return v_[offset_ + i];
  }

  decltype(auto) operator[](std::size_t i) const
  {
    return v_[offset_ + i];
  }

  constexpr static std::size_t size()
  {
    return n;
  }

  constexpr static std::size_t N()
  {
    return n;
  }
};

} /* namespace */

#endif
