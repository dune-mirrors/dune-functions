#include "config.h"
#include <iostream>

#include <dune/common/hybridutilities.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/typetree/treepath.hh>

#include "basisfactories.hh"

using namespace Dune;

template <class T>
void printType()
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
}

template <class T = double, class Traits = ISTLTraits, class PreBasis,
          class Prefix = Dune::TypeTree::HybridTreePath<>>
auto vectorGenerator(PreBasis const& preBasis, Prefix prefix = {})
{
  using namespace Dune::Indices;

  const auto size = preBasis.size(prefix);

  // the coefficient type
  if constexpr(Dune::Functions::isStaticConstant(size) && size == 0) { return T{}; }
  else {
    const auto isUniform = preBasis.isUniform(prefix);
    if constexpr(Dune::Functions::isStaticConstant(isUniform)) {
      if constexpr(isUniform) {
        if constexpr(Dune::Functions::isStaticConstant(size)) {
          // fixed-size vector type
          using V = decltype(vectorGenerator<T,Traits>(preBasis, push_back(prefix, _0)));
          return typename Traits::template PowerVector<V, std::size_t(size)>{};
        }
        else {
          // dynamic-size vector type
          auto subPrefix = push_back(prefix, _0);
          auto subSize = preBasis.size(subPrefix);
          if constexpr (Dune::Functions::isStaticConstant(subSize)) {
            if constexpr(subSize == 0)
              // scalar coefficients
              return typename Traits::template DynamicVector<T>{};
            else
              // leaf-blocked
              return typename Traits::template LeafBlockVector<T, std::size_t(subSize)>{};
          }
          else {
            // hierarchic blocked vector
            using V = decltype(vectorGenerator<T,Traits>(preBasis, subPrefix));
            return typename Traits::template DynamicVector<V>{};
          }
        }
      }
      else {
        if constexpr(Dune::Functions::isStaticConstant(size)) {
          // multi-type vector type
          return Dune::unpackIntegerSequence([&](auto... i) {
            return typename Traits::template CompositeVector<
              decltype(vectorGenerator<T,Traits>(preBasis, push_back(prefix,i)))...>{};
          }, std::make_index_sequence<std::size_t(size)>{});
        }
        else {
          DUNE_THROW(Dune::NotImplemented, "!isUniform && !isStaticConstant(size) not supported");
          return 0;
        }
      }
    }
    else {
      DUNE_THROW(Dune::NotImplemented, "!isStaticConstant(isUniform) not supported");
      return 0;
    }
  }
}


template <class Factory, std::size_t I>
void test(Dune::TestSuite& testSuite, Factory const& factory, index_constant<I> ii)
{
  using namespace Dune::Indices;
  auto basis = factory.basis(ii);
  auto const& preBasis = basis.preBasis();

  using Vector1 = decltype(vectorGenerator<double>(preBasis));    // generate vector
  using Vector2 = decltype(factory.template vector<double>(ii));   // expected type

  std::cout << ii << ":" << std::endl;
  printType<Vector1>();
  testSuite.check(std::is_same_v<Vector1,Vector2>);
  if (!std::is_same_v<Vector1,Vector2>) {
    printType<Vector2>();
  }
  std::cout << std::endl;
}

int main(int argc, char** argv)
{
  MPIHelper::instance(argc, argv);

  using Factory = Dune::Functions::BasisFactories<2>;
  Factory factory{};
  Dune::TestSuite testSuite;

  Hybrid::forEach(range(index_constant<Factory::num_bases>{}),
    [&](auto ii) {
      Dune::TestSuite subTestSuite("basis " + std::to_string(std::size_t(ii)));
      test(subTestSuite, factory, ii);
      testSuite.subTest(subTestSuite);
    });

  return testSuite.exit();
}
