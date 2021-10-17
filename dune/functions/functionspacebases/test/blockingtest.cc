#include "config.h"
#include <iostream>

#include <dune/common/hybridutilities.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/typetree/treepath.hh>

#include "basisfactories.hh"

using namespace Dune;

template <class T>
void printType(const T&)
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
}

template <class PreBasis, class Prefix>
void printInfo(const PreBasis preBasis, const Prefix& prefix)
{
  auto size = preBasis.size(prefix);
  auto isUniform = preBasis.isUniform(prefix);

  std::cout << "  size(" << prefix << ")=" << size
            << " (static=" << Dune::Functions::isIntegralConstant(size) << ")" << std::endl;
  std::cout << "  isUniform(" << prefix << ")=" << isUniform
            << " (static=" << Dune::Functions::isIntegralConstant(isUniform) << ")" << std::endl;
}

template <class Factory, std::size_t I>
void test(Dune::TestSuite& testSuite, Factory const& factory, index_constant<I> ii)
{
  using namespace Dune::Indices;
  using namespace Dune::TypeTree;

  auto basis = factory.basis(ii);
  using Basis = decltype(basis);
  using SizePrefix = typename Basis::SizePrefix;
  using S = typename SizePrefix::value_type;

  auto const& pb = basis.preBasis();

  std::cout << ii << ":" << std::endl;
  printType(pb);

  // check whether size functions are consistent
  testSuite.check(S(pb.size(hybridTreePath())) == pb.size(SizePrefix{}), "{}");
  if constexpr(SizePrefix::capacity() > 0)
    testSuite.check(S(pb.size(hybridTreePath(_0))) == pb.size(SizePrefix{0}), "{0}");
  if constexpr(SizePrefix::capacity() > 1)
    testSuite.check(S(pb.size(hybridTreePath(_0,_0))) == pb.size(SizePrefix{0,0}), "{0,0}");
  if constexpr(SizePrefix::capacity() > 2)
    testSuite.check(S(pb.size(hybridTreePath(_0,_0,_0))) == pb.size(SizePrefix{0,0,0}), "{0,0,0}");
  if constexpr(SizePrefix::capacity() > 1)
    testSuite.check(S(pb.size(hybridTreePath(_0,_1))) == pb.size(SizePrefix{0,1}), "{0,1}");
  if constexpr(SizePrefix::capacity() > 1)
    testSuite.check(S(pb.size(hybridTreePath(_1,_0))) == pb.size(SizePrefix{1,0}), "{1,0}");
  if constexpr(SizePrefix::capacity() > 2)
    testSuite.check(S(pb.size(hybridTreePath(_1,_0,_0))) == pb.size(SizePrefix{1,0,0}), "{1,0,0}");
  if constexpr(SizePrefix::capacity() > 1)
    testSuite.check(S(pb.size(hybridTreePath(_1,_1))) == pb.size(SizePrefix{1,1}), "{1,1}");

  std::cout << std::endl;
}

int main(int argc, char** argv)
{
  MPIHelper::instance(argc, argv);

  Dune::TestSuite testSuite;

  using Factory = Dune::Functions::BasisFactories<2>;
  Factory factory{};
  Hybrid::forEach(range(index_constant<Factory::num_bases>{}),
    [&](auto ii) {
      Dune::TestSuite subTestSuite("basis " + std::to_string(std::size_t(ii)));
      test(subTestSuite, factory, ii);
      testSuite.subTest(subTestSuite);
    });

  return testSuite.exit();
}
