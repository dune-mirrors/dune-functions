
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>

#include <dune/functions/functionspacebases/hermitebasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/albertagrid.hh>

// #include <dune/functions/functionspacebases/test/enabledifferentiabilitycheck.hh>

using namespace Dune;
using namespace Dune::Functions;

int main(int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite test_1d("1d"), test_2d("2d"), test_3d("3d");

  using namespace Dune::Functions::BasisFactory;

  { // 1d
    std::cout << "Hermite test in 1d" << std::endl;
    std::unique_ptr<OneDGrid> grid = StructuredGridFactory<OneDGrid>::createSimplexGrid({0.}, {1.}, {10});

    auto gridView = grid->levelGridView(0);

    {
      auto basis = makeBasis(gridView, hermite());
      test_1d.subTest(checkBasis(basis, EnableContinuityCheck(), EnableDifferentiabilityCheck(),
                                 EnableVertexDifferentiabilityCheck(),
                                 CheckLocalFiniteElementFlag()));
    }
  }

  { //2d
    std::cout << "Hermite test in 2d" << std::endl;

    auto grid = StructuredGridFactory<UGGrid<2>>::createSimplexGrid({0., 0.}, {1., 1.}, {{10, 10}});

    auto gridView = grid->leafGridView();
    std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1) << " facettes and " << gridView.size(2) << " vertices" << std::endl;
    // using GridView = decltype(gridView);
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, hermite());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test_2d.subTest(checkBasis(basis, EnableContinuityCheck(),
                                 EnableVertexDifferentiabilityCheck(),
                                 CheckLocalFiniteElementFlag()));
    }
  }

  { // 2d specialization
    std::cout << "Hermite test in 2d" << std::endl;

    auto grid = StructuredGridFactory<UGGrid<2>>::createSimplexGrid({0., 0.}, {1., 1.}, {{10, 10}});

    auto gridView = grid->leafGridView();
    std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1)
              << " facettes and " << gridView.size(2) << " vertices" << std::endl;
    // using GridView = decltype(gridView);
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, hermite<double,true>());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test_2d.subTest(checkBasis(basis, EnableContinuityCheck(),
                                 EnableVertexDifferentiabilityCheck(),
                                 CheckLocalFiniteElementFlag()));
    }
  }

  { // 2d  with tangential map
    std::cout << "Hermite test in 2d with transcribed tangential" << std::endl;

    auto grid = StructuredGridFactory<UGGrid<2>>::createSimplexGrid({0., 0.}, {1., 1.}, {{10, 10}});

    auto gridView = grid->leafGridView();
    std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1)
              << " facettes and " << gridView.size(2) << " vertices" << std::endl;
    // using GridView = decltype(gridView);
    {
      using namespace Dune::Functions::BasisFactory;
      auto f = [](auto const &x) -> FieldVector<double, 2> {
        return {0., 1.};
      };
      auto basis = makeBasis(gridView, hermite(f));
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test_2d.subTest(checkBasis(basis, EnableContinuityCheck(),
                                 EnableVertexDifferentiabilityCheck(),
                                 CheckLocalFiniteElementFlag()));
    }
  }

  { // 3d
    std::cout << "Hermite test in 3d" << std::endl;

    auto grid = StructuredGridFactory<UGGrid<3>>::createSimplexGrid({0., 0., 0.}, {1., 1., 1.},
                                                                    {{3, 3, 3}});

    auto gridView = grid->leafGridView();
    std::cout << "Grid has " << gridView.size(0) << " elementes and " << gridView.size(1)
              << " facettes and " << gridView.size(2) << " edges and " << gridView.size(3)
              << " vertices " << std::endl;
    {
      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, hermite());
      std::cout << "Basis has " << basis.size() << " Dofs" << std::endl;

      test_3d.subTest(checkBasis(basis, EnableContinuityCheck(),
                                 EnableVertexDifferentiabilityCheck(),
                                 CheckLocalFiniteElementFlag()));
    }
  }

  return test_1d.exit() + test_2d.exit() + test_3d.exit();
}
