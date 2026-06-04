// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_TESTTRANSFORMEDLOCALBASIS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_TESTTRANSFORMEDLOCALBASIS_HH

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/functions/common/densevectorview.hh>
#include <dune/functions/functionspacebases/transformed/concepts.hh>
#include <dune/functions/functionspacebases/transformed/derivative.hh>

namespace Dune::Functions::Test {

namespace Impl {

template<class Basis>
using TestField = std::common_type_t<typename Basis::Domain::value_type,
                                     typename Basis::Traits::RangeFieldType>;

/**
 * \brief Return a scalar norm for arithmetic, vector, or matrix-like range values.
 */
template<class T>
double norm(T const& value)
{
  if constexpr (std::is_arithmetic_v<T>)
    return std::abs(value);
  else if constexpr (requires { value.frobenius_norm(); })
    return value.frobenius_norm();
  else
    return value.two_norm();
}

/**
 * \brief Return the norm of the difference of two derivative values.
 */
template<class T>
double differenceNorm(T const& a, T const& b)
{
  if constexpr (std::is_arithmetic_v<T>)
    return std::abs(a-b);
  else {
    auto difference = a;
    difference -= b;
    return norm(difference);
  }
}

/**
 * \brief Human-readable name for derivative tags used in TestSuite messages.
 */
template<class Derivative>
std::string derivativeName(Derivative)
{
  if constexpr (std::is_same_v<Derivative,Derivatives::Value>)
    return "Value";
  else if constexpr (std::is_same_v<Derivative,Derivatives::Jacobian>)
    return "Jacobian";
  else if constexpr (std::is_same_v<Derivative,Derivatives::Divergence>)
    return "Divergence";
  else if constexpr (std::is_same_v<Derivative,Derivatives::Curl>)
    return "Curl";
  else
    return "derivative";
}

/**
 * \brief Default physical step length for central finite differences.
 *
 * For first derivatives, the central-difference truncation error is O(h^2)
 * while roundoff behaves like eps/h.  Thus eps^(1/3) is a more balanced
 * default than sqrt(eps).
 */
template<class Basis>
double finiteDifferenceEpsilon()
{
  using Field = TestField<Basis>;
  return std::cbrt(static_cast<double>(std::numeric_limits<Field>::epsilon()));
}

/**
 * \brief Default relative tolerance for finite-difference derivative checks.
 */
template<class Basis>
double finiteDifferenceTolerance()
{
  using Field = TestField<Basis>;
  return 100.0 * std::cbrt(static_cast<double>(std::numeric_limits<Field>::epsilon()));
}

/**
 * \brief Default tolerance for comparing split and combined evaluation.
 */
template<class Basis>
double stagedEvaluationTolerance()
{
  using Field = TestField<Basis>;
  return 1000.0 * static_cast<double>(std::numeric_limits<Field>::epsilon());
}

/**
 * \brief Map a physical coordinate perturbation to a local-coordinate step.
 */
template<class Basis, class Geometry>
typename Basis::Domain localStepForPhysicalDirection(Geometry const& geometry,
                                                     typename Basis::Domain const& x,
                                                     int direction,
                                                     double epsilon)
{
  using ctype = typename Basis::Domain::value_type;

  FieldVector<ctype,Geometry::coorddimension> physicalStep(0);
  physicalStep[direction] = epsilon;

  typename Basis::Domain localStep;
  geometry.jacobianInverse(x).mv(physicalStep, localStep);
  return localStep;
}

/**
 * \brief Build central-difference points and check that they stay in the reference element.
 */
template<class Basis, class Geometry>
bool shiftedPointsInside(Geometry const& geometry,
                         GeometryType type,
                         typename Basis::Domain const& x,
                         int direction,
                         double epsilon,
                         typename Basis::Domain& up,
                         typename Basis::Domain& down)
{
  auto localStep = localStepForPhysicalDirection<Basis>(geometry, x, direction, epsilon);
  up = x;
  down = x;
  up += localStep;
  down -= localStep;

  auto const& referenceElement = ReferenceElements<typename Basis::Domain::value_type,Basis::Traits::dimDomain>::general(type);
  return referenceElement.checkInside(up) and referenceElement.checkInside(down);
}

/**
 * \brief Approximate one component of the physical derivative of one shape function.
 */
template<class Basis, class Geometry>
double componentDerivative(Basis const& basis,
                           Geometry const& geometry,
                           GeometryType type,
                           typename Basis::Domain const& x,
                           std::size_t shapeFunction,
                           int component,
                           int direction,
                           double epsilon)
{
  typename Basis::Domain up;
  typename Basis::Domain down;
  if (not shiftedPointsInside<Basis>(geometry, type, x, direction, epsilon, up, down))
    return 0.0;

  std::vector<typename Basis::template DerivativeRange<Derivatives::Value>> upValues;
  std::vector<typename Basis::template DerivativeRange<Derivatives::Value>> downValues;
  basis.evaluate(Derivatives::Value{}, up, upValues);
  basis.evaluate(Derivatives::Value{}, down, downValues);

  auto upView = Dune::Functions::Impl::DenseVectorView(upValues[shapeFunction]);
  auto downView = Dune::Functions::Impl::DenseVectorView(downValues[shapeFunction]);
  return (upView[component] - downView[component])/(2*epsilon);
}

/**
 * \brief Approximate the full Jacobian by central differences of transformed values.
 */
template<class Basis, class Geometry>
typename Basis::template DerivativeRange<Derivatives::Jacobian>
finiteDifference(Basis const& basis,
                 Geometry const& geometry,
                 GeometryType type,
                 Derivatives::Jacobian,
                 typename Basis::Domain const& x,
                 std::size_t shapeFunction,
                 double epsilon)
{
  typename Basis::template DerivativeRange<Derivatives::Jacobian> derivative{};

  for (auto direction : Dune::range(Geometry::coorddimension)) {
    typename Basis::Domain up;
    typename Basis::Domain down;
    if (not shiftedPointsInside<Basis>(geometry, type, x, direction, epsilon, up, down))
      continue;

    std::vector<typename Basis::template DerivativeRange<Derivatives::Value>> upValues;
    std::vector<typename Basis::template DerivativeRange<Derivatives::Value>> downValues;
    basis.evaluate(Derivatives::Value{}, up, upValues);
    basis.evaluate(Derivatives::Value{}, down, downValues);

    auto upView = Dune::Functions::Impl::DenseVectorView(upValues[shapeFunction]);
    auto downView = Dune::Functions::Impl::DenseVectorView(downValues[shapeFunction]);

    if constexpr (requires { derivative[0][0]; }) {
      for (auto component : Dune::range(Basis::Traits::dimRange))
        derivative[component][direction] = (upView[component] - downView[component])/(2*epsilon);
    }
    else
      derivative[direction] = (upView[0] - downView[0])/(2*epsilon);
  }

  return derivative;
}

/**
 * \brief Approximate divergence by summing directional component derivatives.
 */
template<class Basis, class Geometry>
typename Basis::template DerivativeRange<Derivatives::Divergence>
finiteDifference(Basis const& basis,
                 Geometry const& geometry,
                 GeometryType type,
                 Derivatives::Divergence,
                 typename Basis::Domain const& x,
                 std::size_t shapeFunction,
                 double epsilon)
{
  typename Basis::template DerivativeRange<Derivatives::Divergence> divergence{};
  for (auto direction : Dune::range(Geometry::coorddimension))
    divergence += componentDerivative(basis, geometry, type, x, shapeFunction, direction, direction, epsilon);
  return divergence;
}

/**
 * \brief Approximate curl by central differences of transformed values.
 */
template<class Basis, class Geometry>
typename Basis::template DerivativeRange<Derivatives::Curl>
finiteDifference(Basis const& basis,
                 Geometry const& geometry,
                 GeometryType type,
                 Derivatives::Curl,
                 typename Basis::Domain const& x,
                 std::size_t shapeFunction,
                 double epsilon)
{
  typename Basis::template DerivativeRange<Derivatives::Curl> curl{};

  if constexpr (Geometry::coorddimension == 2) {
    curl = componentDerivative(basis, geometry, type, x, shapeFunction, 1, 0, epsilon)
         - componentDerivative(basis, geometry, type, x, shapeFunction, 0, 1, epsilon);
  }
  else {
    static_assert(Geometry::coorddimension == 3,
      "Curl finite-difference test supports only dimension 2 or 3.");

    curl[0] = componentDerivative(basis, geometry, type, x, shapeFunction, 2, 1, epsilon)
            - componentDerivative(basis, geometry, type, x, shapeFunction, 1, 2, epsilon);
    curl[1] = componentDerivative(basis, geometry, type, x, shapeFunction, 0, 2, epsilon)
            - componentDerivative(basis, geometry, type, x, shapeFunction, 2, 0, epsilon);
    curl[2] = componentDerivative(basis, geometry, type, x, shapeFunction, 1, 0, epsilon)
            - componentDerivative(basis, geometry, type, x, shapeFunction, 0, 1, epsilon);
  }

  return curl;
}

/**
 * \brief Check that precompute()+finalize() agrees with evaluate().
 */
template<class Basis, class Derivative>
void checkStagedEvaluation(TestSuite& test,
                           Basis const& basis,
                           Derivative derivative,
                           typename Basis::Domain const& x)
{
  typename Basis::template PrecomputeBuffer<Derivative> precomputed;
  std::vector<typename Basis::template DerivativeRange<Derivative>> finalized;
  std::vector<typename Basis::template DerivativeRange<Derivative>> evaluated;

  basis.precompute(derivative, x, precomputed);
  basis.finalize(derivative, x, precomputed, finalized);
  basis.evaluate(derivative, x, evaluated);

  test.check(finalized.size() == basis.size(), derivativeName(derivative) + " finalize size");
  test.check(evaluated.size() == basis.size(), derivativeName(derivative) + " evaluate size");

  for (std::size_t i = 0; i < std::min(finalized.size(), evaluated.size()); ++i)
    test.check(differenceNorm(finalized[i], evaluated[i]) < stagedEvaluationTolerance<Basis>(),
      derivativeName(derivative) + " staged evaluation consistency");
}

/**
 * \brief Test one derivative tag on all quadrature points of one bound element.
 *
 * Value checks currently verify size and staged evaluation consistency.  Other
 * supported derivative tags are also compared to central finite differences of
 * the transformed values.
 */
template<class Basis, class Element, class Derivative>
TestSuite testDerivative(Basis const& basis,
                         Element const& element,
                         Derivative derivative,
                         unsigned int quadratureOrder,
                         double epsilon,
                         double tolerance)
{
  TestSuite test(derivativeName(derivative));

  auto geometry = element.geometry();
  auto const& quad = QuadratureRules<typename Basis::Domain::value_type,Basis::Traits::dimDomain>::rule(element.type(), quadratureOrder);

  for (auto const& quadraturePoint : quad) {
    auto const& x = quadraturePoint.position();
    checkStagedEvaluation(test, basis, derivative, x);

    if constexpr (not std::is_same_v<Derivative,Derivatives::Value>) {
      std::vector<typename Basis::template DerivativeRange<Derivative>> derivatives;
      basis.evaluate(derivative, x, derivatives);

      for (auto i : Dune::range(derivatives.size())) {
        auto finiteDifferenceValue = finiteDifference(basis, geometry, element.type(), derivative, x, i, epsilon);
        auto scale = std::max({1.0, norm(derivatives[i]), norm(finiteDifferenceValue)});
        auto error = differenceNorm(derivatives[i], finiteDifferenceValue);

        test.check(error <= tolerance*scale, derivativeName(derivative) + " finite difference")
          << "shape function " << i
          << " at " << x
          << ": error " << error
          << ", tolerance " << tolerance*scale;
      }
    }
  }

  return test;
}

} // namespace Impl

/**
 * \brief Test a transformed local finite element for a selected list of derivative tags.
 *
 * The test uses the finite element's physicalBasis(), verifies the staged
 * transformation protocol for every tag, and compares derivative tags to finite
 * differences of transformed values where a finite-difference implementation is
 * provided.
 */
template<class LocalFiniteElement, class Element, class... Derivatives>
TestSuite testTransformedLocalFiniteElement(LocalFiniteElement const& finiteElement,
                                            Element const& element,
                                            Derivatives... derivatives)
{
  TestSuite test("transformed local finite element");

  auto const& basis = finiteElement.physicalBasis();
  using Basis = std::remove_cvref_t<decltype(basis)>;
  static_assert((Concept::TransformedLocalBasis<Basis,Derivatives> && ...));

  test.check(finiteElement.size() == basis.size(), "finite element and physical basis sizes agree");
  (test.subTest(Impl::testDerivative(basis, element, derivatives, 4,
    Impl::finiteDifferenceEpsilon<Basis>(),
    Impl::finiteDifferenceTolerance<Basis>())), ...);

  return test;
}

} // end namespace Dune::Functions::Test

#endif
