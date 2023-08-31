#include "config.h"
#include <iostream>

#include <dune/common/classname.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/functions/functionspacebases/containerdescriptors.hh>

#include "basisfactories.hh"

using namespace Dune;
using namespace Dune::Functions;

template <class CD>
struct VectorGenerator;

template <class T = double, class Traits = Dune::StdTraits, class CD>
auto vectorGenerator(CD const& containerDescriptor)
{
  return VectorGenerator<CD>::template apply<T,Traits>(containerDescriptor);
}


template <std::size_t n>
struct VectorGenerator<ContainerDescriptors::FlatArray<n>>
{
  template <class T, class Traits, class CD>
  static auto apply(CD const& containerDescriptor)
  {
    return FieldVector<T,CD::size()>{};
  }
};

template <>
struct VectorGenerator<ContainerDescriptors::FlatVector>
{
  template <class T, class Traits, class CD>
  static auto apply(CD const& containerDescriptor)
  {
    return Traits::template makeVector<T>(containerDescriptor.size());
  }
};

template <class C, std::size_t n>
struct VectorGenerator<ContainerDescriptors::Array<C,n>>
{
  template <class T, class Traits, class CD>
  static auto apply(CD const& containerDescriptor)
  {
    auto block = [&](auto i) { return vectorGenerator<T,Traits>(containerDescriptor[i]); };
    using Block0 = decltype(block(Dune::index_constant<0>{}));
    return Dune::unpackIntegerSequence([block](auto... ii) {
        return Traits::template makeArray<Block0,CD::size()>(block(ii)...);
      },
      std::make_index_sequence<CD::size()>{});
  }
};

template <class C, std::size_t n>
struct VectorGenerator<ContainerDescriptors::UniformArray<C,n>>
{
  template <class T, class Traits, class CD>
  static auto apply(CD const& containerDescriptor)
  {
    auto block = [&](auto i) { return vectorGenerator<T,Traits>(containerDescriptor[i]); };
    using Block0 = decltype(block(Dune::index_constant<0>{}));
    return Traits::template makeUniformArray<n>(block(0));
  }
};

template <class C>
struct VectorGenerator<ContainerDescriptors::Vector<C>>
{
  template <class T, class Traits, class CD>
  static auto apply(CD const& containerDescriptor)
  {
    auto block = [&](auto i) { return vectorGenerator<T,Traits>(containerDescriptor[i]); };
    using Block0 = decltype(block(Dune::index_constant<0>{}));
    auto container = Traits::template makeVector<Block0>(containerDescriptor.size());
    for (std::size_t i = 0; i < containerDescriptor.size(); ++i)
      container[i] = block(i);
    return container;
  }
};

template <class C>
struct VectorGenerator<ContainerDescriptors::UniformVector<C>>
{
  template <class T, class Traits, class CD>
  static auto apply(CD const& containerDescriptor)
  {
    auto block = [&](auto i) { return vectorGenerator<T,Traits>(containerDescriptor[i]); };
    using Block0 = decltype(block(Dune::index_constant<0>{}));
    return Traits::makeUniformVector(containerDescriptor.size(), block(0));
  }
};

template <class... Cs>
struct VectorGenerator<ContainerDescriptors::Tuple<Cs...>>
{
  template <class T, class Traits, class CD>
  static auto apply(CD const& containerDescriptor)
  {
    auto block = [&](auto i) { return vectorGenerator<T,Traits>(containerDescriptor[i]); };
    return Dune::unpackIntegerSequence([block](auto... ii) {
        return Traits::makeTuple(block(ii)...);
      }, std::make_index_sequence<CD::size()>());
  }
};


template <class GridView, std::size_t I>
void test(Dune::TestSuite& testSuite, GridView const& gridView, index_constant<I> ii)
{
  using namespace Dune::Indices;
  using namespace Dune::Functions::BasisFactory;

  auto basis = makeBasis(gridView, BasisFactories::basis(ii));
  auto const& preBasis = basis.preBasis();

  using Vector1 = decltype(vectorGenerator<double>(preBasis.containerDescriptor()));    // generate vector
  using Vector2 = decltype(BasisFactories::template vector<double>(ii));   // expected type

  std::cout << ii << ":" << std::endl;
  std::cout << Dune::className(preBasis) << std::endl;
  std::cout << Dune::className(preBasis.containerDescriptor()) << std::endl;;
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
  Dune::TestSuite testSuite;

  using Grid = Dune::YaspGrid<2>;
  Grid grid({1.0, 1.0}, {2, 2});

  using namespace Dune::Functions::BasisFactory;
  Hybrid::forEach(range(index_constant<BasisFactories::size>{}),
    [&](auto ii) {
      Dune::TestSuite subTestSuite("basis " + std::to_string(std::size_t(ii)));
      test(subTestSuite, grid.leafGridView(), ii);
      testSuite.subTest(subTestSuite);
    });

  return testSuite.exit();
}
