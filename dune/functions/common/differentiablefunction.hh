// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_HH
#define DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_HH

#include <type_traits>

#include <dune/functions/common/type_traits.hh>
#include <dune/functions/common/defaultderivativetraits.hh>
#include <dune/functions/common/derivativedirection.hh>
#include <dune/functions/common/differentiablefunction_imp.hh>
#include <dune/functions/common/polymorphicsmallobject.hh>
#include <dune/functions/common/concept.hh>
#include <dune/functions/common/signature.hh>

namespace Dune {
namespace Functions {

/**
 * Default implementation is empty
 * The actual implementation is only given if Signature is an type
 * describing a function signature as Range(Domain).
 */
template<class Signature, template<class> class DerivativeTraits=DefaultDerivativeTraits, size_t bufferSize=56>
class DifferentiableFunction
{};



/**
 * \brief Class storing differentiable functions using type erasure
 *
 */
template<typename Range, typename... Domain, template<class> class DerivativeTraits, size_t bufferSize>
class DifferentiableFunction< Range(Domain...), DerivativeTraits, bufferSize>
{
public:

  /**
   * \brief Signature of wrapped functions
   */
  using Signature = Range(Domain...);

  template<int P>
  struct PartialDomain
  {
    using type = typename std::tuple_element< P-1, std::tuple<Domain...> >::type;
  };
  template<typename P>
  struct DerivativeInterface
  {
    using PartialDomain = P;
    using PartialSignature = Range(PartialDomain);
    using RawSignature = typename SignatureTraits<PartialSignature>::RawSignature;
    using DerivativeRange = typename DerivativeTraits<RawSignature>::Range;
    using DerivativeSignature = DerivativeRange(Domain...);
    using type = DifferentiableFunction<DerivativeSignature, DerivativeTraits, bufferSize>;
  };

  using DerivativeInterfaces = std::tuple<typename DerivativeInterface<Domain>::type...>;

  /**
   * \brief Construct from function
   *
   * \tparam F Function type
   *
   * \param f Function of type F
   *
   * Calling derivative(DifferentiableFunction) will result in an exception
   * if the passed function does provide a free derivative() function
   * found via ADL.
   */
  template<class F, disableCopyMove<DifferentiableFunction, F> = 0 >
  DifferentiableFunction(F&& f) :
    f_(Imp::DifferentiableFunctionWrapper<Signature, DerivativeInterfaces, typename std::decay<F>::type>(std::forward<F>(f)))
  {}

  DifferentiableFunction() = default;

  /**
   * \brief Evaluation of wrapped function
   */
  Range operator() (const Domain&... x) const
  {
    return f_.get().operator()(x...);
  }

  /**
   * \brief Get derivative of wrapped function
   *
   * This is a free function that will be found by ADL.
   */
  template<int D = 1>
  friend
  typename DerivativeInterface< typename PartialDomain<D>::type >::type
  derivative(const DifferentiableFunction& t, DerivativeDirection<D> dir = DerivativeDirection<D>())
  {
    return t.f_.get().derivative(dir);
  }

private:
  PolymorphicSmallObject<Imp::DifferentiableFunctionWrapperBase<Signature, DerivativeInterfaces>, bufferSize > f_;
};



}} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_COMMON_DIFFERENTIABLE_FUNCTION_HH
