#include "config.h"
#include <iostream>

#include <dune/common/hybridutilities.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dune/typetree/treepath.hh>

#include "basisfactories.hh"

using namespace Dune;

template <class T>
void printType()
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
}

template <class T>
void printType(const T&)
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
}

template <class T = double, class Traits = Dune::StdTraits, class SizeTree>
auto vectorGenerator(SizeTree const& sizeTree)
{
  using Properties = Dune::Functions::SizeTreeProperties<SizeTree>;

  if constexpr(Properties::isLeaf) {
    if constexpr(Properties::isStatic)
      return FieldVector<T,SizeTree::size()>{};
    else
      return typename Traits::template DynamicVector<T>(sizeTree.size());
  }
  else {
    auto block = [&](auto i) { return vectorGenerator(sizeTree[i]); };
    if constexpr(Properties::isTypeUniform) {
      using Block = decltype(block(Dune::index_constant<0>{}));
      if constexpr(Properties::isStatic) {
        return Dune::unpackIntegerSequence(
          [block,&](auto... ii) {
            return typename Traits::template PowerVector<Block, std::size_t(SizeTree::size())>{block(ii)...};
          },
          std::make_index_sequence<std::size_t(SizeTree::size())>{});
      }
      else {
        typename Traits::template DynamicVector<Block> container;
        container.reserve(sizeTree.size());
        for (std::size_t i = 0; i < sizeTree.size(); ++i)
          container.emplace_back(block(i));
        return container;
      }
    }
    else {
      return Dune::unpackIntegerSequence(
        [block,&](auto... ii) {
          return typename Traits::template CompositeVector<decltype(block(ii))...>{block(ii)...};
        }, std::make_index_sequence<std::size_t(SizeTree::size())>());
    }
  }
}


struct ISTLTraits
{
  template <class T>
  using DynamicVector = BlockVector<T>;

  template <class B, std::size_t N>
  using PowerVector = BlockVector<B>;

  template <class... Rows>
  using CompositeVector = MultiTypeBlockVector<Rows...>;
};


template <class Factory, std::size_t I>
void test(Dune::TestSuite& testSuite, Factory const& factory, index_constant<I> ii)
{
  using namespace Dune::Indices;
  auto basis = factory.basis(ii);
  auto const& preBasis = basis.preBasis();

  using Vector1 = decltype(vectorGenerator<double>(sizeTree(preBasis)));    // generate vector
  using Vector2 = decltype(factory.template vector<double>(ii));   // expected type

  std::cout << ii << ":" << std::endl;
  printType(preBasis);
  printType(sizeTree(preBasis));
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
