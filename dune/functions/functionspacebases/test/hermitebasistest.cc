// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>
//#define HERMITE_INTERPOLATION_VARIANT_A 1
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/stringutility.hh>

#include <dune/grid/onedgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/common/differentiablefunctionfromcallables.hh>

#include <dune/functions/functionspacebases/boundarydofs.hh>
//#include <dune/functions/functionspacebases/reducedcubichermitetrianglebasis.hh>
#include <dune/functions/functionspacebases/test/cubichermitebasis.hh>
#include <dune/functions/functionspacebases/test/interpolatetest.hh>
#include <dune/common/timer.hh>

#include <dune/istl/matrixindexset.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/io.hh>



using namespace Dune::Functions;



template<class Geometry>
auto quadratureRule(const Geometry& geometry, unsigned int order)
{
  constexpr auto dim = Geometry::mydimension;
  const auto& rule = Dune::QuadratureRules<double, dim>::rule(geometry.type(), order);
  return Dune::transformedRangeView(rule, [&](auto&& qp) {
    return std::pair(qp.position(), qp.weight());
  });
}



template<class LocalView>
auto indices(const LocalView& localView)
{
  auto size = localView.size();
  return Dune::transformedRangeView(Dune::range(size), [&](auto localIndex) {
    return std::pair(localIndex, localView.index(localIndex));
  });
}



template<class Basis, class K>
void setupMatrixPattern(const Basis& basis, Dune::BCRSMatrix<K>& matrix) {
  auto matrixIndexSet = Dune::MatrixIndexSet(basis.size(), basis.size());
  auto localView = basis.localView();
  for(const auto& e : Dune::elements(basis.gridView()))
  {
    localView.bind(e);
    for (auto [local_i, i] : indices(localView))
      for (auto [local_j, j] : indices(localView))
        matrixIndexSet.add(i,j);
  }
  matrixIndexSet.exportIdx(matrix);
  matrix = 0;
}



template<class Basis, class K>
void assembleMassMatrix(const Basis& basis, Dune::BCRSMatrix<K>& matrix) {
  using Range = typename Basis::LocalView::Tree::FiniteElement::Traits::LocalBasisType::Traits::RangeType;
  std::vector<Range> values;
  auto localView = basis.localView();
  matrix = 0;
  for(const auto& e : Dune::elements(basis.gridView()))
  {
    localView.bind(e);
    auto geometry = e.geometry();
    const auto& localBasis = localView.tree().finiteElement().localBasis();
    auto order = 2*localBasis.order();
    for(auto [position, weight] : quadratureRule(geometry, order))
    {
      localBasis.evaluateFunction(position, values);
      auto integrationElement = geometry.integrationElement(position);
      for (auto [local_i, i] : indices(localView))
        for (auto [local_j, j] : indices(localView))
          matrix[i][j] += values[local_i] * values[local_j] * integrationElement * weight;
    }
  }
}



template<class Basis, class K>
void assembleDirichletStiffnessMatrix(const Basis& basis, Dune::BCRSMatrix<K>& matrix) {
  using Jacobian = typename Basis::LocalView::Tree::FiniteElement::Traits::LocalBasisType::Traits::JacobianType;
  std::vector<Jacobian> jacobians;
  auto localView = basis.localView();
  matrix = 0;
  for(const auto& e : Dune::elements(basis.gridView()))
  {
    localView.bind(e);
    auto geometry = e.geometry();
    const auto& localBasis = localView.tree().finiteElement().localBasis();
    auto order = 2*localBasis.order();
    for(auto [position, weight] : quadratureRule(geometry, order))
    {
      localBasis.evaluateJacobian(position, jacobians);
      auto integrationElement = geometry.integrationElement(position);
      auto jacobianInverse = geometry.jacobianInverse(position);
      for (auto [local_i, i] : indices(localView))
        for (auto [local_j, j] : indices(localView))
          matrix[i][j] += (jacobians[local_i] * jacobianInverse)[0] * (jacobians[local_j] * jacobianInverse)[0] * integrationElement * weight;
    }
  }
  auto isBoundary = std::vector<bool>(basis.size(), false);
  Dune::Functions::forEachBoundaryDOF(basis, [&] (auto&& index) {
    isBoundary[index] = true;
  });
  for(auto&& [m_i, i] : Dune::sparseRange(matrix))
  {
    if (isBoundary[i])
      for(auto&& [m_ij, j] : Dune::sparseRange(m_i))
        if (isBoundary[j])
          m_ij = (i==j);
  }
}



template<class Basis, class F>
void benchmarkBasis(const Basis& basis, const F& f, int repeat, std::string name)
{
  std::cout << "------------------------------------------------------------" << std::endl;
  std::cout << "Benchmarks for " << name << ":" << std::endl;
  auto timer = Dune::Timer();

  timer.reset();
  auto x = std::vector<double>();
  {
    for (auto i : Dune::range(repeat))
      Dune::Functions::interpolate(basis, x, f);
  }
  std::cout << "Time spend on doing " << Dune::formatString("%4d", repeat) << " interpolations:      " << timer.elapsed() << std::endl;

  auto matrix = Dune::BCRSMatrix<double>();

  timer.reset();
  setupMatrixPattern(basis, matrix);
  std::cout << "Time spend on setting up sparsity pattern:    " << timer.elapsed() << std::endl;

  timer.reset();
  assembleMassMatrix(basis, matrix);
  std::cout << "Time spend on assembling mass matrix:         " << timer.elapsed() << std::endl;

  writeMatrixToMatlab(matrix, std::string("massmatrix_")+name);

  timer.reset();
  assembleDirichletStiffnessMatrix(basis, matrix);
  std::cout << "Time spend on assembling H1 stiffness matrix: " << timer.elapsed() << std::endl;

  writeMatrixToMatlab(matrix, std::string("stiffmatrix_")+name);
}




int main (int argc, char* argv[])
{
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite suite;

  using namespace Dune::Functions::BasisFactory;

  {
    auto gridPtr = Dune::StructuredGridFactory<Dune::OneDGrid>::createCubeGrid({0}, {1}, {10});

    auto gridView = gridPtr->leafGridView();
    auto basis = makeBasis(gridView, cubicHermite());

//    suite.subTest(checkBasis(basis, EnableContinuityCheck(), EnableVertexJacobianContinuityCheck()));
    suite.subTest(checkBasis(basis, EnableContinuityCheck()));

    {
      using Domain = Dune::FieldVector<double,1>;
      using Range = Dune::FieldVector<double,1>;

      using SignatureTag = Dune::Functions::SignatureTag<Range(Domain)>;
      auto f = [](const auto& x) {
        return x*x;
      };
      auto df = [](const auto& x) {
        return 2*x;
      };
      auto ff = Dune::Functions::makeDifferentiableFunctionFromCallables(SignatureTag(), f, df);
      auto coefficients = std::vector<double>();
      Dune::Functions::interpolate(basis, coefficients, ff);

      suite.subTest(checkInterpolateConsistency<Range>(basis, coefficients));

#if NDEBUG
      auto fGridFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<Range>(basis, coefficients);

      benchmarkBasis(basis, fGridFunction, 100, "1d");
#endif
    }
  }
  {
    auto factory = Dune::GridFactory<Dune::UGGrid<2>>{};
    factory.insertVertex({0,0});
    factory.insertVertex({0,1});
    factory.insertVertex({2,0});
    factory.insertVertex({1,2});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {0,1,2});
    factory.insertElement(Dune::GeometryTypes::simplex(2), {1,2,3});
    auto gridPtr = factory.createGrid();
    auto& grid = *gridPtr;
#if NDEBUG
    grid.globalRefine(7);
#else
    grid.globalRefine(3);
#endif

    using Domain = Dune::FieldVector<double,2>;
    using Range = Dune::FieldVector<double,1>;
    using Jacobian = Dune::FieldMatrix<double,1,2>;
    using SignatureTag = Dune::Functions::SignatureTag<Range(Domain)>;

    auto f = [](const auto& x) {
      return x[0]*x[0] + x[1]*x[1];
    };
    auto df = [](const auto& x) {
      return Jacobian({{2*x[0],2*x[1]}});
    };
    auto ff = Dune::Functions::makeDifferentiableFunctionFromCallables(SignatureTag(), f, df);


    {
      auto basis = makeBasis(grid.leafGridView(), cubicHermite());
//      suite.subTest(checkBasis(basis, EnableContinuityCheck(), EnableVertexJacobianContinuityCheck()));
      suite.subTest(checkBasis(basis, EnableContinuityCheck()));

      auto coefficients = std::vector<double>();
      Dune::Functions::interpolate(basis, coefficients, ff);
      suite.subTest(checkInterpolateConsistency<Range>(basis, coefficients));

#if NDEBUG
      auto fGridFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<Range>(basis, coefficients);

      benchmarkBasis(basis, fGridFunction, 10, "2d");
#endif
    }

    {
      auto basis = makeBasis(grid.leafGridView(), reducedCubicHermite());
//      suite.subTest(checkBasis(basis, EnableContinuityCheck(), EnableVertexJacobianContinuityCheck()));
      suite.subTest(checkBasis(basis, EnableContinuityCheck()));

      auto coefficients = std::vector<double>();
      Dune::Functions::interpolate(basis, coefficients, ff);
      suite.subTest(checkInterpolateConsistency<Range>(basis, coefficients));

#if NDEBUG
      auto fGridFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<Range>(basis, coefficients);

      benchmarkBasis(basis, fGridFunction, 10, "2d_reduced");
#endif
    }

//    {
//      auto basis = makeBasis(grid.leafGridView(), reducedCubicHermiteTriangle());
//      suite.subTest(checkBasis(basis, EnableContinuityCheck(), EnableVertexJacobianContinuityCheck()));
//
//      auto coefficients = std::vector<double>();
//      Dune::Functions::interpolate(basis, coefficients, ff);
//      suite.subTest(checkInterpolateConsistency<Range>(basis, coefficients));
//
//#if NDEBUG
//      auto fGridFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<Range>(basis, coefficients);
//
//      benchmarkBasis(basis, fGridFunction, 10, "2d_reduced_397");
//#endif
//    }

  }

  return suite.exit();
}
