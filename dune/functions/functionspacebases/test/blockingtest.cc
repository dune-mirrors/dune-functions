#include "config.h"
#include <iostream>

#include <dune/common/hybridutilities.hh>
#include <dune/typetree/treepath.hh>

#include "basisfactories.hh"

using namespace Dune;

template <class T>
void printType()
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

template <class Tester, std::size_t I>
void test(Tester const& tester, index_constant<I> ii)
{
  using namespace Dune::Indices;
  auto basis = tester.basis(ii);
  auto const& preBasis = basis.preBasis();

  std::cout << ii << ":" << std::endl;
  printType<decltype(basis)>();
  printInfo(preBasis, TypeTree::hybridTreePath());
  printInfo(preBasis, TypeTree::hybridTreePath(_0));
  printInfo(preBasis, TypeTree::hybridTreePath(_0,_0));
  printInfo(preBasis, TypeTree::hybridTreePath(_0,_0,_0));
  printInfo(preBasis, TypeTree::hybridTreePath(_0,_1));
  printInfo(preBasis, TypeTree::hybridTreePath(_1,_0));
  printInfo(preBasis, TypeTree::hybridTreePath(_1,_0,_0));
  printInfo(preBasis, TypeTree::hybridTreePath(_1,_1));
  std::cout << std::endl;
}

int main(int argc, char** argv)
{
  MPIHelper::instance(argc, argv);

  using Tester = Dune::Functions::BasisFactories<2>;
  Tester tester{};
  Hybrid::forEach(range(index_constant<Tester::num_bases>{}),
    [&tester](auto ii) {
      test(tester, ii);
    });
}
