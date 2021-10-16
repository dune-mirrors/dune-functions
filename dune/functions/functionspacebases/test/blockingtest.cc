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
            << " (static=" << Dune::Functions::isStaticConstant(size) << ")" << std::endl;
  std::cout << "  isUniform(" << prefix << ")=" << isUniform
            << " (static=" << Dune::Functions::isStaticConstant(isUniform) << ")" << std::endl;
}

template <class Tester, std::size_t I>
void test(Tester const& tester, index_constant<I> ii)
{
  auto basis = tester.basis(ii);
  auto const& preBasis = basis.preBasis();

  printType<decltype(basis)>();
  printInfo(preBasis, TypeTree::hybridTreePath());
}



template <class Tester>
void test2(Tester const& tester)
{
  using namespace Dune::Indices;
  auto basis = tester.basis(index_constant<3>{});
  auto const& preBasis = basis.preBasis();

  printType<decltype(preBasis)>();
  printInfo(preBasis, TypeTree::hybridTreePath());
  printInfo(preBasis, TypeTree::hybridTreePath(_0));
  printInfo(preBasis, TypeTree::hybridTreePath(_0,_0));
  printInfo(preBasis, TypeTree::hybridTreePath(_0,_0,_0));
  printInfo(preBasis, TypeTree::hybridTreePath(_0,_1));
  printInfo(preBasis, TypeTree::hybridTreePath(_1,_0));
  printInfo(preBasis, TypeTree::hybridTreePath(_1,_0,_0));
  printInfo(preBasis, TypeTree::hybridTreePath(_1,_1));
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

  test2(tester);
}
