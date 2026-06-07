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
  else if constexpr (std::is_same_v<Derivative,Derivatives::Gradient>)
    return "Gradient";
  else if constexpr (std::is_same_v<Derivative,Derivatives::Hessian>)
    return "Hessian";
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
 * \brief Approximate the ambient tangential gradient of a scalar shape function.
 */
template<class Basis, class Geometry>
typename Basis::template DerivativeRange<Derivatives::Gradient>
finiteDifference(Basis const& basis,
                 Geometry const& geometry,
                 GeometryType type,
                 Derivatives::Gradient,
                 typename Basis::Domain const& x,
                 std::size_t shapeFunction,
                 double epsilon)
{
  static_assert(Basis::Traits::dimRange == 1);
  auto jacobian = finiteDifference(
    basis,geometry,type,Derivatives::Jacobian{},x,shapeFunction,epsilon);
  typename Basis::template DerivativeRange<Derivatives::Gradient> gradient{};
  for (auto direction : Dune::range(Geometry::coorddimension))
    gradient[direction] = jacobian[0][direction];
  return gradient;
}

/**
 * \brief Approximate the physical Hessian of a scalar shape function.
 */
template<class Basis, class Geometry>
typename Basis::template DerivativeRange<Derivatives::Hessian>
finiteDifference(Basis const& basis,
                 Geometry const& geometry,
                 GeometryType type,
                 Derivatives::Hessian,
                 typename Basis::Domain const& x,
                 std::size_t shapeFunction,
                 double epsilon)
{
  static_assert(Basis::Traits::dimRange == 1,
    "Hessian finite-difference test currently supports scalar basis functions.");

  typename Basis::template DerivativeRange<Derivatives::Hessian> hessian{};
  std::vector<typename Basis::template DerivativeRange<Derivatives::Value>> centerValues;
  basis.evaluate(Derivatives::Value{}, x, centerValues);
  auto center = Dune::Functions::Impl::DenseVectorView(centerValues[shapeFunction])[0];

  for (auto i : Dune::range(Geometry::coorddimension)) {
    typename Basis::Domain up;
    typename Basis::Domain down;
    if (shiftedPointsInside<Basis>(geometry, type, x, i, epsilon, up, down)) {
      std::vector<typename Basis::template DerivativeRange<Derivatives::Value>> upValues;
      std::vector<typename Basis::template DerivativeRange<Derivatives::Value>> downValues;
      basis.evaluate(Derivatives::Value{}, up, upValues);
      basis.evaluate(Derivatives::Value{}, down, downValues);
      auto upValue = Dune::Functions::Impl::DenseVectorView(upValues[shapeFunction])[0];
      auto downValue = Dune::Functions::Impl::DenseVectorView(downValues[shapeFunction])[0];
      hessian[i][i] = (upValue - 2*center + downValue)/(epsilon*epsilon);
    }

    for (auto j : Dune::range(i)) {
      auto stepI = localStepForPhysicalDirection<Basis>(geometry, x, i, epsilon);
      auto stepJ = localStepForPhysicalDirection<Basis>(geometry, x, j, epsilon);
      std::array<typename Basis::Domain,4> points{x, x, x, x};
      points[0] += stepI;
      points[0] += stepJ;
      points[1] += stepI;
      points[1] -= stepJ;
      points[2] -= stepI;
      points[2] += stepJ;
      points[3] -= stepI;
      points[3] -= stepJ;

      auto const& referenceElement = ReferenceElements<
        typename Basis::Domain::value_type,
        Basis::Traits::dimDomain>::general(type);
      if (std::all_of(points.begin(), points.end(), [&](auto const& point) { return referenceElement.checkInside(point); })) {
        std::array<double,4> values;
        for (auto k : Dune::range(4)) {
          std::vector<typename Basis::template DerivativeRange<Derivatives::Value>> pointValues;
          basis.evaluate(Derivatives::Value{}, points[k], pointValues);
          values[k] = Dune::Functions::Impl::DenseVectorView(pointValues[shapeFunction])[0];
        }
        hessian[i][j] = (values[0] - values[1] - values[2] + values[3])/(4*epsilon*epsilon);
        hessian[j][i] = hessian[i][j];
      }
    }
  }

  return hessian;
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

  static_assert(Geometry::mydimension == 2 || Geometry::mydimension == 3,
    "Curl finite-difference test supports only intrinsic dimension 2 or 3.");

  using Field = typename Basis::Traits::RangeFieldType;
  using ReferenceVector = FieldVector<Field,Geometry::mydimension>;
  std::array<ReferenceVector,Geometry::mydimension> derivatives{};
  auto const& referenceElement = ReferenceElements<
    typename Basis::Domain::value_type,
    Basis::Traits::dimDomain>::general(type);

  for (auto direction : Dune::range(Geometry::mydimension)) {
    auto up = x;
    auto down = x;
    up[direction] += epsilon;
    down[direction] -= epsilon;
    if (!referenceElement.checkInside(up) || !referenceElement.checkInside(down))
      continue;

    auto covariantValue = [&](auto const& point) {
      std::vector<typename Basis::template DerivativeRange<Derivatives::Value>> values;
      basis.evaluate(Derivatives::Value{},point,values);
      auto value = Dune::Functions::Impl::DenseVectorView(values[shapeFunction]);
      ReferenceVector covariant;
      geometry.jacobianTransposed(point).mv(value,covariant);
      return covariant;
    };

    auto upValue = covariantValue(up);
    auto downValue = covariantValue(down);
    derivatives[direction] = upValue;
    derivatives[direction] -= downValue;
    derivatives[direction] /= 2*epsilon;
  }

  auto integrationElement = geometry.integrationElement(x);
  if constexpr (Geometry::mydimension == 2)
    curl = (derivatives[0][1] - derivatives[1][0])/integrationElement;
  else {
    ReferenceVector referenceCurl;
    referenceCurl[0] = derivatives[1][2] - derivatives[2][1];
    referenceCurl[1] = derivatives[2][0] - derivatives[0][2];
    referenceCurl[2] = derivatives[0][1] - derivatives[1][0];
    geometry.jacobian(x).mv(referenceCurl,curl);
    curl /= integrationElement;
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
