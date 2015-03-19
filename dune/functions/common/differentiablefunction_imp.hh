// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH
#define DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH

#include <tuple>

#include <dune/common/exceptions.hh>

#include <dune/functions/common/signature.hh>
#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/interfaces.hh>
#include <dune/functions/common/derivativedirection.hh>

#include "concept.hh"

namespace Dune {
namespace Functions {
namespace Imp {

/**
 * A concept describing types that have a derivative() method found by ADL
 */
struct HasFreeDerivative
{
  template<class F, int D>
  auto require(F&& f, DerivativeDirection<D>&& d) -> decltype(
    derivative(f,d)
  );
};



template<class Dummy, class F, int D,
  typename std::enable_if<
    Dune::Functions::Concept::models< HasFreeDerivative, F>() , int>::type = 0>
auto derivativeIfImplemented(const F& f) -> decltype(derivative(f))
{
  return derivative(f);
}



template<class Dummy, class F, int D,
  typename std::enable_if<
    not(Dune::Functions::Concept::models< HasFreeDerivative, F>()) , int>::type = 0>
Dummy derivativeIfImplemented(const F& f, DerivativeDirection<D> d)
{
  DUNE_THROW(Dune::NotImplemented, "Derivative not implemented");
}

template<typename DerivativeInterfaces, int P>
class PartialDerivativeWrapperBase;
template<typename... DerivativeInterfaces>
class PartialDerivativeWrapperBase<std::tuple<DerivativeInterfaces...>, 0> {};
template<typename... DerivativeInterfaces, int P>
class PartialDerivativeWrapperBase<std::tuple<DerivativeInterfaces...>, P> :
  public PartialDerivativeWrapperBase<std::tuple<DerivativeInterfaces...>, P-1>
{

  /**
   * \brief Wrapper type of returned derivatives
   */
  using DerivativeInterface = typename std::tuple_element< P-1, std::tuple<DerivativeInterfaces...> >::type;

public:

  /**
   * \brief compute partial derivative wrt P'th parameter
   */
  virtual DerivativeInterface derivative(DerivativeDirection<P> d) const = 0;

};

template<class Signature, class DerivativePack>
class DifferentiableFunctionWrapperBase;

template<typename Range, typename... Domain, typename... DerivativeInterfaces>
class DifferentiableFunctionWrapperBase<Range(Domain...), std::tuple<DerivativeInterfaces...> > :
  public PolymorphicType<DifferentiableFunctionWrapperBase<Range(Domain...), std::tuple<DerivativeInterfaces...> > >,
  public PartialDerivativeWrapperBase<std::tuple<DerivativeInterfaces...>, sizeof...(DerivativeInterfaces)>
{
  static_assert( sizeof...(Domain) == sizeof...(DerivativeInterfaces), "Type Mismatch");
public:
  virtual Range operator() (const Domain&... x) const = 0;
};



template<typename Signature, typename DerivativeInterfaces, int P, class WrapperImp>
class PartialDerivativeWrapper;

template<typename Signature, typename... DerivativeInterfaces, class WrapperImp>
class PartialDerivativeWrapper<Signature, std::tuple<DerivativeInterfaces...>, 0, WrapperImp> :
    public DifferentiableFunctionWrapperBase<Signature, std::tuple<DerivativeInterfaces...> >
{};

template<typename Signanture, typename... DerivativeInterfaces, int P, class WrapperImp>
class PartialDerivativeWrapper<Signanture, std::tuple<DerivativeInterfaces...>, P, WrapperImp> :
  public PartialDerivativeWrapper<Signanture, std::tuple<DerivativeInterfaces...>, P-1, WrapperImp>
{

  /**
   * \brief Wrapper type of returned derivatives
   */
  using DerivativeInterface = typename std::tuple_element< P-1, std::tuple<DerivativeInterfaces...> >::type;

public:

  /**
   * \brief compute partial derivative wrt P'th parameter
   */
  virtual DerivativeInterface derivative(DerivativeDirection<P> d) const
  {
    auto f_ = static_cast<const WrapperImp*>(this)->f_;
    using FImp = decltype(f_);
    return derivativeIfImplemented<DerivativeInterface, FImp>(f_,d);
  };

};

template<class Signature, class DerivativeInterfaces, class FImp>
class DifferentiableFunctionWrapper;

template<typename Range, typename... Domain, typename... DerivativeInterfaces, class FImp>
class DifferentiableFunctionWrapper< Range(Domain...), std::tuple<DerivativeInterfaces...>, FImp> :
    public PartialDerivativeWrapper< Range(Domain...), std::tuple<DerivativeInterfaces...>, sizeof...(DerivativeInterfaces),
                                    DifferentiableFunctionWrapper< Range(Domain...), std::tuple<DerivativeInterfaces...>, FImp> >
{
  static_assert( sizeof...(Domain) == sizeof...(DerivativeInterfaces), "Type Mismatch");

public:

  template<class F, disableCopyMove<DifferentiableFunctionWrapper, F> = 0>
  DifferentiableFunctionWrapper(F&& f) :
    f_(std::forward<F>(f))
  {}

  virtual Range operator() (const Domain&... x) const
  {
    return f_(x...);
  }

  virtual DifferentiableFunctionWrapper* clone() const
  {
    return new DifferentiableFunctionWrapper(*this);
  }

  virtual DifferentiableFunctionWrapper* clone(void* buffer) const
  {
    return new (buffer) DifferentiableFunctionWrapper(*this);
  }

  virtual DifferentiableFunctionWrapper* move(void* buffer)
  {
    return new (buffer) DifferentiableFunctionWrapper(std::move(*this));
  }

// private:
  FImp f_;
};



}}} // namespace Dune::Functions::Imp



#endif // DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_IMP_HH
