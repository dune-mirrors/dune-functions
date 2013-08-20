// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_GRIDVIEWFUNCTION_HH
#define DUNE_FUNCTIONS_COMMON_GRIDVIEWFUNCTION_HH

#include <memory>
#include <dune/functions/common/differentiablefunction.hh>

namespace Dune {

namespace Functions {



template<typename GV, typename RT>
class GridViewFunction
  : public DifferentiableFunction<typename GV::template Codim<0>::Geometry::GlobalCoordinate,RT>
{

  typedef DifferentiableFunction<
    typename GV::template Codim<0>::Geometry::GlobalCoordinate,
    RT
    > Base;

  typedef GV GridView;
  typedef typename Base::Domain Domain;
  typedef typename Base::Range Range;
  typedef typename Base::DerivativeRange DerivativeRange;

  typedef GridViewFunction<GV,DerivativeRange> Derivative;


  typedef typename GV::template Codim<0>::Geometry::LocalCoordinate LocalDomain;
  typedef typename GV::template Codim<0>::Entity Element;

  virtual void evaluate(const Element& e, const LocalDomain& coord, Range& r) const = 0;

  virtual Derivative* derivative() const = 0;

};



} // end of namespace Dune::Functions
} // end of namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_GRIDVIEWFUNCTION_HH
