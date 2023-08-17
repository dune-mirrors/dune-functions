#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_VECTORFLOATCMP_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_VECTORFLOATCMP_HH

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

namespace Dune::Functions::Impl
{

  /**
   * @brief lexicographic float-comparison equal for same size vectors
   *
   * @param left
   * @param right
   * @return
   */
  template <class GlobalCoordinate>
  inline bool vectorEqual(GlobalCoordinate const &left, GlobalCoordinate const &right)
  {
    assert(left.size() == right.size()); // Coordinate must be of same length
    return Dune::FloatCmp::eq(left, right);
  }

  /**
   * @brief lexicographic float-comparison Less for same size vectors
   *
   * @param left
   * @param right
   * @return
   */
  template <class GlobalCoordinate>
  bool vectorLess(GlobalCoordinate const &left, GlobalCoordinate const &right)
    {
    assert(left.size() == right.size()); // Coordinate must be of same length

    for (std::size_t i = 0; i < left.size(); ++i)
    {
      if (Dune::FloatCmp::ge(left[i], right[i]))
      {
        if (Dune::FloatCmp::eq(left[i], right[i]))
          continue;
        else
        return false;
    }
      else
    return true;
  }
    return false;
  }

  /**
   * @brief lexicographic float-comparison Greater for same size vectors
   *
   * @param left
   * @param right
   * @return
   */
  template <class GlobalCoordinate>
  bool vectorGreater(GlobalCoordinate const &left, GlobalCoordinate const &right)
  {
    assert(left.size() == right.size()); // Coordinate must be of same length

    for (std::size_t i = 0; i < left.size(); ++i)
    {
      if (Dune::FloatCmp::le(left[i], right[i]))
      {
        if (Dune::FloatCmp::eq(left[i], right[i]))
          continue;
        else
          return false;
      }
      else
        return true;
    }
    return false; // full equality
  }

} // namespace Dune::Functions::Impl

#endif