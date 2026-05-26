#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HELLANHERRMANNJOHNSONREFERENCE0_INC_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HELLANHERRMANNJOHNSONREFERENCE0_INC_HH

namespace Dune::Functions {
namespace Impl {

template<class D, class R, int dim, unsigned int k>
void HellanHerrmannJohnsonReferenceLocalBasis<D, R, dim, k>::evaluateFunction(
    typename Traits::DomainType const& x, std::vector<typename Traits::RangeType>& out) const
{
   out.resize(size());
   // ...
   using Range = typename Traits::RangeType;
   if constexpr(k == 0) {
     R half = R(1)/2;
     out[0] = Range({{0,-half},{-half,1}});
     out[1] = Range({{1,-half},{-half,0}});
     out[2] = Range({{0, half},{ half,0}});
   }
   else if constexpr(k == 1) {
     out[0] = sym<Range>(0, 3*x[0]+3*x[1]-2, -6*x[0]-6*x[1]+4);
     out[1] = sym<Range>(0, 1-3*x[0], 6*x[0]-2);
     out[2] = sym<Range>(-6*x[0]-6*x[1]+4, 3*x[0]+3*x[1]-2, 0);
     out[3] = sym<Range>(6*x[1]-2, 1-3*x[1], 0);
     out[4] = sym<Range>(0, 3*x[0]-1, 0);
     out[5] = sym<Range>(0, 3*x[1]-1, 0);
     out[6] = sym<Range>(3*x[0], -15*x[0]/2-15*x[1]/2+6, 3*x[1]);
     out[7] = sym<Range>(-3*x[0], 3*x[0]+3*x[1]/2-R(3)/2, 0);
     out[8] = sym<Range>(0, -3*x[0]/2-3*x[1]+R(3)/2, 3*x[1]);
  }
  else if constexpr(k == 2) {
    R xx = x[0]*x[0];
    R xy = x[0]*x[1];
    R yy = x[1]*x[1];

    // edge 0
    out[0] = sym<Range>(0, -15*xx - 30*xy + 18*x[0] - 15*yy + 18*x[1] - R(9)/2,
                         30*xx + 60*xy - 36*x[0] + 30*yy - 36*x[1] + 9);
    out[1] = sym<Range>(0, -15*xx + 12*x[0] - R(3)/2, 30*xx - 24*x[0] + 3);
    out[2] = sym<Range>(0, 15*xx/2 + 15*xy/2 - 15*x[0]/2 - 15*yy/4 + 3*x[1]/2 + R(3)/4,
                         -15*xx - 15*xy + 15*x[0] + 15*yy/2 - 3*x[1] - R(3)/2);

    // edge 1
    out[3] = sym<Range>(30*xx + 60*xy - 36*x[0] + 30*yy - 36*x[1] + 9, -15*xx - 30*xy + 18*x[0] - 15*yy + 18*x[1] - R(9)/2, 0);
    out[4] = sym<Range>(30*yy - 24*x[1] + 3, -15*yy + 12*x[1] - R(3)/2, 0);
    out[5] = sym<Range>(15*xx/2 - 15*xy - 3*x[0] - 15*yy + 15*x[1] - R(3)/2, -15*xx/4 + 15*xy/2 + 3*x[0]/2 + 15*yy/2 - 15*x[1]/2 + R(3)/4, 0);

    // edge 2
    out[6] = sym<Range>(0, 15*xx - 12*x[0] + R(3)/2, 0);
    out[7] = sym<Range>(0, 15*yy - 12*x[1] + R(3)/2, 0);
    out[8] = sym<Range>(0, 15*xx/4 + 15*xy - 6*x[0] + 15*yy/4 - 6*x[1] + R(3)/2, 0);

    // cell
    out[9] = sym<Range>(-60*xx - 60*xy + 48*x[0], 90*xx + 180*xy - 120*x[0] + 90*yy - 120*x[1] + 36, -60*xy - 60*yy + 48*x[1]);
    out[10] = sym<Range>(60*xx + 60*xy - 48*x[0], -45*xx - 60*xy + 48*x[0] - 15*yy + 24*x[1] - 9, 0);
    out[11] = sym<Range>(0, 15*xx + 60*xy - 24*x[0] + 45*yy - 48*x[1] + 9, -60*xy - 60*yy + 48*x[1]);
    out[12] = sym<Range>(30*xx - 12*x[0], -135*xx - 150*xy + 150*x[0] + 30*x[1] - 24, 60*xy - 12*x[1]);
    out[13] = sym<Range>(-30*xx + 12*x[0], 45*xx + 30*xy - 42*x[0] - 6*x[1] + 6, 0);
    out[14] = sym<Range>(0, -30*xx - 60*xy + 36*x[0] + 12*x[1] - 6, 60*xy - 12*x[1]);
    out[15] = sym<Range>(60*xy - 12*x[0], -150*xy + 30*x[0] - 135*yy + 150*x[1] - 24, 30*yy - 12*x[1]);
    out[16] = sym<Range>(-60*xy + 12*x[0], 60*xy - 12*x[0] + 30*yy - 36*x[1] + 6, 0);
    out[17] = sym<Range>(0, -30*xy + 6*x[0] - 45*yy + 42*x[1] - 6, 30*yy - 12*x[1]);
  }
}

template<class D, class R, int dim, unsigned int k>
void HellanHerrmannJohnsonReferenceLocalBasis<D, R, dim, k>::evaluateDivDiv(
    typename Traits::DomainType const& x, std::vector<typename Traits::DivDivType>& out) const
{
   out.resize(size());
   if constexpr(k < 2)
     std::fill(out.begin(), out.end(), R(0));
   else if constexpr(k == 2) {
     out[0] = 0; out[1] = 0; out[2] = 30;
     out[3] = 0; out[4] = 0; out[5] = 30;
     out[6] = 0; out[7] = 0; out[8] = 30;
     out[9]  =  120; out[10] = 0;   out[11] = 0;
     out[12] = -240; out[13] = 0;   out[14] = -120;
     out[15] = -240; out[16] = 120; out[17] = 0;
   }
}

} // namespace Impl
} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_HELLANHERRMANNJOHNSONREFERENCE0_INC_HH
