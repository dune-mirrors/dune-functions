#include "config.h"
#include <iostream>

#include <dune/common/classname.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/test/testsuite.hh>

#include "basisfactories.hh"

using namespace Dune;

template <class ST>
inline constexpr bool isFlat = std::is_same_v<
  std::decay_t<decltype(std::declval<ST>()[Dune::Indices::_0])>, Dune::Functions::EmptyIndexTree>;

template <class T = double, class IndexTree>
auto vectorGenerator(IndexTree const& indexTree)
{
  using Traits = Dune::StdTraits;
  using namespace Dune::Functions;
  if constexpr(isFlat<IndexTree>) {
    if constexpr(HasStaticSize_v<IndexTree>)
      return FieldVector<T,IndexTree::size()>{};
    else
      return Traits::DynamicVector<T>(indexTree.size());
  }
  else { // !flat
    auto block = [&](auto i) { return vectorGenerator<T>(indexTree[i]); };
    if constexpr(IndexTree::isTypeUniform) {
      using Block = decltype(block(Dune::index_constant<0>{}));
      if constexpr(HasStaticSize_v<IndexTree>) {
        return Dune::unpackIntegerSequence([block](auto... ii) {
            return Traits::PowerVector<Block, IndexTree::size()>{block(ii)...};
          },
          std::make_index_sequence<IndexTree::size()>{});
      }
      else { // !IndexTree::isStatic
        Traits::DynamicVector<Block> container;
        container.reserve(indexTree.size());
        for (std::size_t i = 0; i < indexTree.size(); ++i)
          container.emplace_back(block(i));
        return container;
      }
    }
    else {
      if constexpr(HasStaticSize_v<IndexTree>) {
        return Dune::unpackIntegerSequence([block](auto... ii) {
            return Traits::CompositeVector<decltype(block(ii))...>{block(ii)...};
          }, std::make_index_sequence<IndexTree::size()>());
      } else {
        DUNE_THROW(Dune::NotImplemented, "Dynamic NonUniformIndexTrees not implemented.");
      }
    }
  }
}

template <class Factory, std::size_t I>
void test(Dune::TestSuite& testSuite, Factory const& factory, index_constant<I> ii)
{
  using namespace Dune::Indices;
  auto basis = factory.basis(ii);
  auto const& preBasis = basis.preBasis();

  using Vector1 = decltype(vectorGenerator<double>(preBasis.indexTree()));    // generate vector
  using Vector2 = decltype(factory.template vector<double>(ii));   // expected type

  std::cout << ii << ":" << std::endl;
  std::cout << Dune::className(preBasis) << std::endl;
  std::cout << Dune::className(preBasis.indexTree()) << std::endl;;
  std::cout << Dune::className<Vector1>() << std::endl;;
  testSuite.check(std::is_same_v<Vector1,Vector2>);
  if (!std::is_same_v<Vector1,Vector2>) {
    std::cout << Dune::className<Vector2>() << std::endl;;
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
