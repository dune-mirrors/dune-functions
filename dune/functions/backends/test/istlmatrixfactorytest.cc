// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/backends/istlmatrixfactory.hh>
#include <dune/functions/functionspacebases/containerdescriptors.hh>


namespace CD = Dune::Functions::ContainerDescriptors;

template <class T>
void print_type (const T&)
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
}

void checkMatrixFactory (Dune::TestSuite& test)
{
  using namespace Dune::Indices;

  auto mat1 = Dune::Functions::istlMatrixFactory(CD::FlatVector{10},CD::FlatVector{10});
  static_assert(std::is_same_v<decltype(mat1),Dune::BCRSMatrix<double>>);

  auto mat2 = Dune::Functions::istlMatrixFactory(CD::UniformVector<CD::FlatArray<2>>{10}, CD::UniformVector<CD::FlatArray<2>>{10});
  static_assert(std::is_same_v<decltype(mat2),Dune::BCRSMatrix<Dune::FieldMatrix<double,2,2>>>);

  // more complicated test
  CD::Tuple<CD::Array<CD::FlatVector,3>,CD::FlatVector> stokes{
    CD::Array<CD::FlatVector,3>{CD::FlatVector{10},CD::FlatVector{10},CD::FlatVector{10}}, CD::FlatVector{5} };
  auto mat3 = Dune::Functions::istlMatrixFactory(stokes, stokes);
  print_type(mat3);

  auto pattern3 = Dune::Functions::istlMatrixFactory<Dune::Functions::ContainerDescriptors::ISTLSparsityPatternFactory>(stokes, stokes);
  print_type(pattern3);

}

int main (int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;
  checkMatrixFactory(test);

  return test.exit();
}
