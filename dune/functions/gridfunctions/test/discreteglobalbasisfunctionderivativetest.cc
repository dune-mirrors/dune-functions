// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/functions/common/functionconcepts.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/argyrisbasis.hh>
#include <dune/functions/functionspacebases/cubichermitebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/morleybasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/functions/gridfunctions/composedgridfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

struct Difference2
{
  template<class K, int n>
  double operator()(Dune::FieldVector<K, n> x, Dune::FieldVector<K, n> y) const
  {
    return (x -= y).two_norm2();
  }

  template<class K, int n, int m>
  double operator()(Dune::FieldMatrix<K, n, m> x, Dune::FieldMatrix<K, n, m> y) const
  {
    return (x -= y).frobenius_norm2();
  }
};

template<class Function, class Reference>
Dune::TestSuite compare(const Function& function, const Reference& reference, const int quadOrder, const double tol)
{
  Dune::TestSuite test;

  double err2 = 0.0;
  double err3 = 0.0;
  const auto& gridView = function.basis().gridView();
  const int dim = std::decay_t<decltype(gridView)>::dimension;

  const auto f = Dune::Functions::makeComposedGridFunction(
    Difference2(),
    function, reference
    );
  auto flocal = localFunction(f);

  for (const auto& e : elements(gridView)) {
    const auto geometry = e.geometry();
    const auto& quad = Dune::QuadratureRules<double, dim>::rule(e.type(), quadOrder);
    flocal.bind(e);
    for (const auto& qp : quad) {
      const auto x = qp.position();
      const auto integrationElement = geometry.integrationElement(x);
      err2 += flocal(x) * qp.weight() * integrationElement;

      // evaluate the functions in global coordinates
      const auto X = geometry.global(x);
      err3 += f(X) * qp.weight() * integrationElement;
    }
  }

  const double err = 0.5*(std::sqrt(err2) + std::sqrt(err3));

  std::cout << "err = " << err << "\n";
  test.check(err <= tol);

  return test;
}

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite test;

  const int dim = 2;
  const int order = 3;
  auto grid = Dune::YaspGrid<dim>({1., 1.}, {10, 42});
  const auto gridView = leafGridView(grid);

  using namespace Dune::Functions::BasisFactory;

  // scalar Lagrange basis with scalar coefficients
  {
    const auto f = [](auto&& x) -> double {
      return 42. * x[0] * x[0] + 13. * x[1] * x[1] + 7 * x[0] * x[1];
    };
    const auto fprime = Dune::Functions::makeAnalyticGridViewFunction(
      [](auto&& x) -> Dune::FieldVector<double, dim> {
        return {
          84. * x[0] + 7 * x[1],
          26. * x[1] + 7 * x[0]
        };
      },
      gridView);
    const auto basis = makeBasis(gridView, lagrange<order>());

    auto coefficients = std::vector<double>();
    Dune::Functions::interpolate(basis, coefficients, f);

    auto f2 = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(basis, coefficients);
    auto f2prime = derivative(f2);
    static_assert(Dune::Functions::Concept::isDifferentiableGridViewFunction<
                  decltype(f2),
                  double(Dune::FieldVector<double, 2>),
                  std::decay_t<decltype(gridView)>>());
    static_assert(Dune::Functions::Concept::isDifferentiableGridViewFunction<
                  decltype(f2prime),
                  Dune::FieldVector<double, 2>(Dune::FieldVector<double, 2>),
                  std::decay_t<decltype(gridView)>>());

    // `order` should be enough; `order+1` is more than enough.
    // The tolerance is ~100 times the error observed when writing this test.
    test.subTest(compare(f2prime, fprime, order+1, 6.6e-11));
  }

  // scalar Lagrange basis with vector coefficients
  {
    const auto f = [](auto&& x) -> Dune::FieldVector<double, 3> {
      return {
        x[0] * x[0] * x[0] + 23. * x[0] * x[1],
        6. * x[0] + 9000. * x[1],
        9001. * x[0] * x[0] + 17. * x[1]
      };
    };
    const auto fprime = Dune::Functions::makeAnalyticGridViewFunction(
      [](auto&& x) -> Dune::FieldMatrix<double, 3, 2> {
        return {
          { 3. * x[0] * x[0] + 23. * x[1], 23. * x[0] },
          { 6., 9000. },
          { 18002. * x[0], 17. }
        };
      },
      gridView);
    const auto basis = makeBasis(gridView, power<3>(lagrange<order>(), blockedInterleaved()));

    auto coefficients = std::vector< Dune::FieldVector<double, 3> >();
    Dune::Functions::interpolate(basis, coefficients, f);

    auto f2 = Dune::Functions::makeDiscreteGlobalBasisFunction< Dune::FieldVector<double, 3> >(basis, coefficients);
    auto f2prime = derivative(f2);
    static_assert(Dune::Functions::Concept::isDifferentiableGridViewFunction<
                  decltype(f2),
                  Dune::FieldVector<double, 3>(Dune::FieldVector<double, 2>),
                  std::decay_t<decltype(gridView)>>());
    static_assert(Dune::Functions::Concept::isDifferentiableGridViewFunction<
                  decltype(f2prime),
                  Dune::FieldMatrix<double, 3, 2>(Dune::FieldVector<double, 2>),
                  std::decay_t<decltype(gridView)>>());

    // `order` should be enough; `order+1` is more than enough.
    // The tolerance is ~100 times the error observed when writing this test.
    test.subTest(compare(f2prime, fprime, order+1, 1.7e-8));
  }

  // Scalar transformed bases already provide physical Jacobians. Their
  // discrete derivatives must not apply the geometry chain rule a second time.
  {
    using Grid = Dune::UGGrid<2>;
    auto transformedGrid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(
      {0.,0.},{1.,1.},{2,2});
    auto transformedGridView = transformedGrid->leafGridView();
    auto checkTransformedBasis = [&](auto basis) {
      std::vector<double> coefficients(basis.size());
      for (std::size_t i = 0; i < coefficients.size(); ++i)
        coefficients[i] = 0.25 + 0.125*i;
      auto discrete = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(
        basis,coefficients);
      auto discreteDerivative = derivative(discrete);
      auto localDerivative = localFunction(discreteDerivative);
      constexpr double epsilon = 1e-6;
      for (auto const& element : elements(transformedGridView)) {
        localDerivative.bind(element);
        auto const& rule = Dune::QuadratureRules<double,2>::rule(
          element.type(),2);
        for (auto const& point : rule) {
          auto x = element.geometry().global(point.position());
          auto finiteDifference = Dune::FieldVector<double,2>{};
          for (int direction = 0; direction < 2; ++direction) {
            auto up = x;
            auto down = x;
            up[direction] += epsilon;
            down[direction] -= epsilon;
            finiteDifference[direction] =
              (discrete(up)-discrete(down))/(2*epsilon);
          }
          auto transformedDerivative = localDerivative(point.position());
          test.check(
            (transformedDerivative-finiteDifference).two_norm() < 1e-7,
            "transformed discrete derivative agrees with finite differences");
        }
      }
    };

    checkTransformedBasis(makeBasis(transformedGridView,argyris()));
    checkTransformedBasis(makeBasis(transformedGridView,morley()));
    checkTransformedBasis(makeBasis(transformedGridView,cubicHermite()));
  }

  return test.exit();
}
