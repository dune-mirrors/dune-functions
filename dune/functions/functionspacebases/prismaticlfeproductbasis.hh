// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PRISMATICLFEPRODUCTBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PRISMATICLFEPRODUCTBASIS_HH

#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dune/localfunctions/meta/product.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lfeprebasismixin.hh>

namespace Dune {
  namespace Functions {

    /**
     * \brief A pre-basis for a prismatic product of two local finite-element bases
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     * \tparam GV  The grid view that the FE basis is defined on
     * \tparam LFE1  Type of the first local finite-element
     * \tparam LFE2  Type of the second local finite-element
     */
    template<class GV, class LFE1, class LFE2>
    class PrismaticLFEProductPreBasis
        : public LFEPreBasisMixin<GV, PrismaticProduct<LFE1, LFE2>>
    {
      using Base = LFEPreBasisMixin<GV, PrismaticProduct<LFE1, LFE2>>;

    public:
      PrismaticLFEProductPreBasis (const GV& gridView, const LFE1& lfe1, const LFE2& lfe2) :
        PrismaticLFEProductPreBasis(gridView, PrismaticProduct<LFE1, LFE2>{lfe1, lfe2})
      {}

    protected:
      PrismaticLFEProductPreBasis (const GV& gridView, const PrismaticProduct<LFE1, LFE2>& lfe) :
        Base(gridView, lfe, lfe.localCoefficients().layout())
      {}
    };


    namespace BasisFactory {

      /**
       * \brief A factory that can create a HierarchicalLagrange pre-basis
       *
       * \ingroup FunctionSpaceBasesImplementations
       *
       * \tparam LFE1  Type of the first local finite-element
       * \tparam LFE2  Type of the second local finite-element
       */
      template<class LFE1, class LFE2>
      auto prismaticLFEProduct(const LFE1& lfe1, const LFE2& lfe2)
      {
        return [&](const auto& gridView) {
          return PrismaticLFEProductPreBasis{gridView, lfe1, lfe2};
        };
      }

    } // end namespace BasisFactory

    /** \brief Basis of a prismatic product of two bases
     *
     * \ingroup FunctionSpaceBasesImplementations
     *
     * \tparam GV The GridView that the space is defined on
     * \tparam LFE1  Type of the first local finite-element
     * \tparam LFE2  Type of the second local finite-element
     */
    template<class GV, class LFE1, class LFE2>
    using PrismaticLFEProductBasis = DefaultGlobalBasis<PrismaticLFEProductPreBasis<GV, LFE1, LFE2> >;

  } // end namespace Functions
} // end namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_PRISMATICLFEPRODUCTBASIS_HH
