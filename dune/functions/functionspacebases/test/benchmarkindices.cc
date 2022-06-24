// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <array>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/common/timer.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/backends/istlvectorbackend.hh>

#include <dune/functions/functionspacebases/cachedglobalbasis.hh>

template<class Basis>
void benchmarkIndices(int n, const Basis& basis)
{
  std::cout << "MultiIndex type : " << Dune::className<typename Basis::MultiIndex>() << std::endl;

  {
    auto localView = basis.localView();
    auto timer = Dune::Timer();
    long unsigned int sum = 0;
    for (std::size_t i = 0; i<n; ++i)
    {
      for (const auto& e : Dune::elements(basis.gridView()))
      {
        localView.bind(e);
        for (std::size_t j = 0; j<localView.size(); ++j)
          sum += localView.index(j)[0];
      }
    }
    std::cout << "Using local index cache               :" << n << " times took: " << timer.elapsed()
      << "       " << sum << std::endl;
  }

  if constexpr (std::is_same_v<typename Basis::MultiIndex, Dune::Functions::StaticMultiIndex<std::size_t, 1>>)
  {
    auto localView = basis.localView();
    auto timer = Dune::Timer();
    long unsigned int sum = 0;
    for (std::size_t i = 0; i<n; ++i)
    {
      for (const auto& e : Dune::elements(basis.gridView()))
      {
        localView.bind(e);
        auto index = localView.index(0)[0];
        for (std::size_t j = 0; j<localView.size(); ++j)
        {
          index++;
          sum += index;
        }
      }
    }
    std::cout << "Incrementing from local index cache   :" << n << " times took: " << timer.elapsed()
      << "       " << sum << std::endl;
  }

  auto cachedBasis = Dune::Functions::CachedGlobalBasis(basis);
  {
    auto localView = cachedBasis.localView();
    //  double sum = 0;
    auto timer = Dune::Timer();
    long unsigned int sum = 0;
    for (std::size_t i = 0; i<n; ++i)
    {
      for (const auto& e : Dune::elements(cachedBasis.gridView()))
      {
        localView.bind(e);
        for (std::size_t j = 0; j<localView.size(); ++j)
          sum += localView.index(j)[0];
      }
    }
    std::cout << "Using global index cache              :" << n << " times took: " << timer.elapsed()
      << "       " << sum << std::endl;
  }

  if constexpr (std::is_same_v<typename Basis::MultiIndex, Dune::Functions::StaticMultiIndex<std::size_t, 1>>)
  {
    auto localView = cachedBasis.localView();
    auto timer = Dune::Timer();
    long unsigned int sum = 0;
    for (std::size_t i = 0; i<n; ++i)
    {
      for (const auto& e : Dune::elements(cachedBasis.gridView()))
      {
        localView.bind(e);
        auto index = localView.index(0)[0];
        for (std::size_t j = 0; j<localView.size(); ++j)
        {
          index++;
          sum += index;
        }
      }
    }
    std::cout << "Incrementing from global index cache  :" << n << " times took: " << timer.elapsed()
      << "       " << sum << std::endl;
  }
}


int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test;

  // Generate grid for testing
  const int dim = 2;
  using Grid = Dune::YaspGrid<dim>;
  Dune::FieldVector<double,dim> l(1);
  std::array<int,dim> elements = {{10, 10}};
  Grid grid(l,elements);

  using namespace Dune::Functions::BasisFactory;

  grid.globalRefine(5);
  auto gridView = grid.leafGridView();

  benchmarkIndices(10,
      makeBasis(gridView,
        power<3>(
          power<3>(
            composite(
              lagrange<3>(),
              lagrange<1>(),
              flatLexicographic()),
            flatLexicographic()),
          blockedInterleaved())
      )
    );

  benchmarkIndices(10,
      makeBasis(gridView,
        composite(
          power<3>(
            power<3>(
              composite(
                lagrange<3>(),
                lagrange<1>(),
                flatLexicographic()),
              flatLexicographic()),
            blockedInterleaved()),
          lagrange<1>()
        )
      )
    );

  benchmarkIndices(10,
      makeBasis(gridView,
        power<dim>(lagrange<2>(), flatLexicographic())
      )
    );

  benchmarkIndices(10,
      makeBasis(gridView,
        power<dim>(lagrange<2>(), flatInterleaved())
      )
    );

  benchmarkIndices(10,
      makeBasis(gridView,
        power<dim>(lagrange<2>(), blockedLexicographic())
      )
    );

  benchmarkIndices(10,
      makeBasis(gridView,
        power<dim>(lagrange<2>(), blockedInterleaved())
      )
    );

  benchmarkIndices(10,
      makeBasis(gridView,
        composite(
          power<dim>(lagrange<2>(), flatLexicographic()),
          lagrange<1>(),
          flatLexicographic()
        )
      )
    );

  benchmarkIndices(10,
      makeBasis(gridView,
        composite(
          power<dim>(lagrange<2>(), flatInterleaved()),
          lagrange<1>(),
          flatLexicographic()
        )
      )
    );

  benchmarkIndices(10,
      makeBasis(gridView,
        composite(
          power<dim>(lagrange<2>(), blockedLexicographic()),
          lagrange<1>(),
          flatLexicographic()
        )
      )
    );

  benchmarkIndices(10,
      makeBasis(gridView,
        composite(
          power<dim>(lagrange<2>(), blockedInterleaved()),
          lagrange<1>(),
          flatLexicographic()
        )
      )
    );

  benchmarkIndices(10,
      makeBasis(gridView,
        composite(
          power<dim>(lagrange<2>(), flatLexicographic()),
          lagrange<1>(),
          blockedLexicographic()
        )
      )
    );

  benchmarkIndices(10,
      makeBasis(gridView,
        composite(
          power<dim>(lagrange<2>(), flatInterleaved()),
          lagrange<1>(),
          blockedLexicographic()
        )
      )
    );

  benchmarkIndices(10,
      makeBasis(gridView,
        composite(
          power<dim>(lagrange<2>(), blockedLexicographic()),
          lagrange<1>(),
          blockedLexicographic()
        )
      )
    );

  benchmarkIndices(10,
      makeBasis(gridView,
        composite(
          power<dim>(lagrange<2>(), blockedInterleaved()),
          lagrange<1>(),
          blockedLexicographic()
        )
      )
    );

  benchmarkIndices(10,
      makeBasis(gridView, lagrangeDG<5>())
    );

  return test.exit();
}
