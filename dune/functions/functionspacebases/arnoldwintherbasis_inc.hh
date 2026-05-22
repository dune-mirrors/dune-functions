#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ARNOLDWINTHERREFERENCE_INC_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_ARNOLDWINTHERREFERENCE_INC_HH
namespace Dune::Functions
{
namespace Impl
{
template<class D, class R, int dim, unsigned int k>
void ArnoldWintherReferenceLocalBasis<D, R, dim, k>::evaluateFunction(
    typename Traits::DomainType const &in, std::vector<typename Traits::RangeType> &out) const
{
  out.resize(size());
  auto iter = out.begin();

  // generated with sympy from symfem library
  auto const &x = in[0], y = in[1];
  if constexpr (k == 2) {

    // 0th basis function
    *(iter++) = sym<Range>(x * x * (2 * x + 9 * y - 3) + y * (y * (18 - 10 * y) - 9) + 1,
                           x * (-6 * x * y + y * (6 - 9 * y)), 6 * x * y * y + y * y * (3 * y - 3));

    // 1th basis function
    *(iter++) = sym<Range>(x * (x * (15 * x + 45 * y - 24) + y * (30 * y - 36) + 9),
                           x * (x * (-10 * x - 45 * y + 18) + y * (48 - 45 * y) - 9) +
                               y * (y * (18 - 10 * y) - 9) + 1,
                           x * (30 * x * y + y * (45 * y - 36)) + y * (y * (15 * y - 24) + 9));

    // 2th basis function
    *(iter++) = sym<Range>(x * x * (3 * x + 6 * y - 3), x * (-9 * x * y + y * (6 - 6 * y)),
                           x * (x * (18 - 10 * x) + 9 * y * y - 9) + y * y * (2 * y - 3) + 1);

    // 3th basis function
    *(iter++) = sym<Range>(x * x * (-2 * x - 9 * y + 3), x * (6 * x * y + y * (9 * y - 6)),
                           -6 * x * y * y + y * y * (3 - 3 * y));

    // 4th basis function
    *(iter++) = sym<Range>(x * x * (-9 * x - 18 * y + 9),
                           x * (x * (10 * x + 27 * y - 12) + y * (18 * y - 18) + 3),
                           x * (-30 * x * y + y * (24 - 27 * y)) + y * (y * (9 - 6 * y) - 3));

    // 5th basis function
    *(iter++) = sym<Range>(0, 0, x * (x * (10 * x - 12) + 3));

    // 6th basis function
    *(iter++) = sym<Range>(y * (y * (10 * y - 12) + 3), 0, 0);

    // 7th basis function
    *(iter++) = sym<Range>(x * (x * (-6 * x - 27 * y + 9) + y * (24 - 30 * y) - 3),
                           x * (18 * x * y + y * (27 * y - 18)) + y * (y * (10 * y - 12) + 3),
                           -18 * x * y * y + y * y * (9 - 9 * y));

    // 8th basis function
    *(iter++) = sym<Range>(x * x * (-3 * x - 6 * y + 3), x * (9 * x * y + y * (6 * y - 6)),
                           -9 * x * y * y + y * y * (3 - 2 * y));

    // 9th basis function
    *(iter++) =
        sym<Range>(x * x * (-18 * x - 36 * y + 18), x * (54 * x * y + y * (36 * y - 36)),
                   x * (x * (60 * x - 96) + y * (24 - 54 * y) + 36) + y * (y * (30 - 12 * y) - 18));

    // 10th basis function
    *(iter++) = sym<Range>(x * (x * (-60 * x - 90 * y + 78) + 24 * y - 18),
                           x * (x * (60 * x + 180 * y - 96) + y * (90 * y - 132) + 36),
                           x * (-180 * x * y + y * (192 - 180 * y)) + y * (y * (66 - 30 * y) - 36));

    // 11th basis function
    *(iter++) =
        sym<Range>(x * x * (18 * x + 36 * y - 18), x * (-54 * x * y + y * (36 - 36 * y)),
                   x * (x * (84 - 60 * x) + y * (54 * y - 24) - 24) + y * (y * (12 * y - 18) + 6));

    // 12th basis function
    *(iter++) = sym<Range>(x * (x * (60 * x + 90 * y - 66) + 6),
                           x * (x * (-60 * x - 180 * y + 84) + y * (108 - 90 * y) - 24),
                           x * (180 * x * y + y * (180 * y - 168)) + y * (y * (30 * y - 54) + 24));

    // 13th basis function
    *(iter++) =
        sym<Range>(x * (x * (-12 * x - 54 * y + 30) + 24 * y - 18) + y * (y * (60 * y - 96) + 36),
                   x * (36 * x * y + y * (54 * y - 36)), -36 * x * y * y + y * y * (18 - 18 * y));

    // 14th basis function
    *(iter++) = sym<Range>(x * (x * (30 * x + 180 * y - 66) + y * (180 * y - 192) + 36),
                           x * (-90 * x * y + y * (132 - 180 * y)) + y * (y * (96 - 60 * y) - 36),
                           x * y * (90 * y - 24) + y * (y * (60 * y - 78) + 18));

    // 15th basis function
    *(iter++) =
        sym<Range>(x * (x * (12 * x + 54 * y - 18) - 24 * y + 6) + y * (y * (84 - 60 * y) - 24),
                   x * (-36 * x * y + y * (36 - 54 * y)), 36 * x * y * y + y * y * (18 * y - 18));

    // 16th basis function
    *(iter++) = sym<Range>(x * (x * (-30 * x - 180 * y + 54) + y * (168 - 180 * y) - 24),
                           x * (90 * x * y + y * (180 * y - 108)) + y * (y * (60 * y - 84) + 24),
                           -90 * x * y * y + y * (y * (66 - 60 * y) - 6));

    // 17th basis function
    *(iter++) = sym<Range>(x * (x * (-3 * x + 9 * y + 6) - 3), x * (9 * x * y - 9 * y * y),
                           x * y * (12 - 9 * y) + y * (3 * y * y - 3));

    // 18th basis function
    *(iter++) =
        sym<Range>(x * (x * (15 * x + 45 * y - 12) - 3), x * (-45 * x * y + y * (36 - 45 * y)),
                   x * y * (45 * y - 12) + y * (y * (15 * y - 18) + 3));

    // 19th basis function
    *(iter++) = sym<Range>(x * (x * (3 * x - 9 * y) + 12 * y - 3), x * (-9 * x * y + 9 * y * y),
                           9 * x * y * y + y * (y * (6 - 3 * y) - 3));

    // 20th basis function
    *(iter++) = sym<Range>(x * (x * (-15 * x - 45 * y + 18) + 12 * y - 3),
                           x * (45 * x * y + y * (45 * y - 36)),
                           -45 * x * y * y + y * (y * (12 - 15 * y) + 3));

    // 21th basis function
    *(iter++) = sym<Range>(x * (-24 * x - 24 * y + 24), 0, 0);

    // 22th basis function
    *(iter++) =
        sym<Range>(x * (-24 * x - 48 * y + 24), 24 * x * y, -48 * x * y + y * (24 - 24 * y));

    // 23th basis function
    *(iter++) = sym<Range>(0, 0, -24 * x * y + y * (24 - 24 * y));
  }
  else
    DUNE_THROW(NotImplemented, "Higher order not implemented");
}
template<class D, class R, int dim, unsigned int k>
void ArnoldWintherReferenceLocalBasis<D, R, dim, k>::evaluateDivergence(
    typename Traits::DomainType const &in, std::vector<typename Traits::DivergenceType> &out) const
{
  out.resize(size());
  auto iter = out.begin();

  // generated with sympy from symfem library
  auto const &x = in[0], y = in[1];
  if constexpr (k == 2) {
    // 0th basis function
    *(iter++) = {0, 0};

    // 1th basis function
    *(iter++) = {0, 0};

    // 2th basis function
    *(iter++) = {0, 0};

    // 3th basis function
    *(iter++) = {0, 0};

    // 4th basis function
    *(iter++) = {0, 0};

    // 5th basis function
    *(iter++) = {0, 0};

    // 6th basis function
    *(iter++) = {0, 0};

    // 7th basis function
    *(iter++) = {0, 0};

    // 8th basis function
    *(iter++) = {0, 0};

    // 9th basis function
    *(iter++) = {0, 24 * x + 24 * y - 18};

    // 10th basis function
    *(iter++) = {24 * x + 24 * y - 18, 0};

    // 11th basis function
    *(iter++) = {0, 6 - 24 * x};

    // 12th basis function
    *(iter++) = {6 - 24 * x, 0};

    // 13th basis function
    *(iter++) = {24 * x + 24 * y - 18, 0};

    // 14th basis function
    *(iter++) = {0, -24 * x - 24 * y + 18};

    // 15th basis function
    *(iter++) = {6 - 24 * y, 0};

    // 16th basis function
    *(iter++) = {0, 24 * y - 6};

    // 17th basis function
    *(iter++) = {12 * x - 3, 12 * x - 3};

    // 18th basis function
    *(iter++) = {12 * x - 3, 3 - 12 * x};

    // 19th basis function
    *(iter++) = {12 * y - 3, 12 * y - 3};

    // 20th basis function
    *(iter++) = {12 * y - 3, 3 - 12 * y};

    // 21th basis function
    *(iter++) = {-48 * x - 24 * y + 24, 0};

    // 22th basis function
    *(iter++) = {-24 * x - 48 * y + 24, -48 * x - 24 * y + 24};

    // 23th basis function
    *(iter++) = {0, -24 * x - 48 * y + 24};
  }
}
} // namespace Impl
} // namespace Dune::Functions
#endif
