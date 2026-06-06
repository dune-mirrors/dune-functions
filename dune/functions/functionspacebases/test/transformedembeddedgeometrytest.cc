// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <array>
#include <cmath>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/common/localbasis.hh>

#include <dune/functions/functionspacebases/transformed/basisset.hh>
#include <dune/functions/functionspacebases/transformed/bindcontext.hh>
#include <dune/functions/functionspacebases/transformed/localfiniteelement.hh>
#include <dune/functions/functionspacebases/transformed/piola.hh>

namespace {

struct CustomDerivative {};

struct CustomReferenceEvaluator
{
  template<class Derivative, class LocalBasis>
  using Range = typename LocalBasis::Traits::RangeType;

  template<class LocalBasis, class Out>
  void evaluate(CustomDerivative,
                LocalBasis const& localBasis,
                typename LocalBasis::Traits::DomainType const& x,
                Out& out) const
  {
    localBasis.evaluateFunction(x,out);
    for (auto& value : out)
      value *= 2;
  }
};

struct PassThroughStage
{
  template<class Derivative, class LocalBasis, class Context, class InputRange>
  using OutputRange = InputRange;

  template<class Context>
  void bind(Context const&)
  {}

  template<class Derivative, class LocalBasis, class InputRange, class OutputRange>
  void transform(Derivative,
                 LocalBasis const&,
                 typename LocalBasis::Traits::DomainType const&,
                 std::vector<InputRange> const& in,
                 std::vector<OutputRange>& out) const
  {
    out.assign(in.begin(),in.end());
  }
};

template<int dim>
class AffineVectorLocalBasis
{
public:
  using Traits = Dune::LocalBasisTraits<
    double,dim,Dune::FieldVector<double,dim>,
    double,dim,Dune::FieldVector<double,dim>,Dune::FieldMatrix<double,dim,dim>>;

  std::size_t size() const { return 1; }
  int order() const { return 1; }

  void evaluateFunction(typename Traits::DomainType const& x,
                        std::vector<typename Traits::RangeType>& out) const
  {
    static_assert(dim == 2);
    out.resize(1);
    out[0][0] = x[0] + x[1];
    out[0][1] = 2*x[0] + 3*x[1];
  }

  void evaluateJacobian(typename Traits::DomainType const&,
                        std::vector<typename Traits::JacobianType>& out) const
  {
    static_assert(dim == 2);
    out.resize(1);
    out[0][0][0] = 1;
    out[0][0][1] = 1;
    out[0][1][0] = 2;
    out[0][1][1] = 3;
  }
};

class AffineScalarLocalBasis
{
public:
  using Traits = Dune::LocalBasisTraits<
    double,1,Dune::FieldVector<double,1>,
    double,1,Dune::FieldVector<double,1>,Dune::FieldMatrix<double,1,1>>;

  std::size_t size() const { return 1; }
  int order() const { return 1; }

  void evaluateFunction(typename Traits::DomainType const& x,
                        std::vector<typename Traits::RangeType>& out) const
  {
    out.resize(1);
    out[0][0] = x[0];
  }

  void evaluateJacobian(typename Traits::DomainType const&,
                        std::vector<typename Traits::JacobianType>& out) const
  {
    out.resize(1);
    out[0][0][0] = 1;
  }
};

class ScalarLocalBasis2
{
public:
  using Traits = Dune::LocalBasisTraits<
    double,2,Dune::FieldVector<double,2>,
    double,1,Dune::FieldVector<double,1>,Dune::FieldMatrix<double,1,2>>;
};

bool close(double a, double b)
{
  return std::abs(a-b) < 1e-12;
}

} // namespace

int main()
{
  Dune::TestSuite test;

  using SurfaceGeometry = Dune::MultiLinearGeometry<double,2,3>;
  std::vector<Dune::FieldVector<double,3>> surfaceCorners{
    Dune::FieldVector<double,3>{0,0,0},
    Dune::FieldVector<double,3>{2,0,0},
    Dune::FieldVector<double,3>{0,3,0}};
  SurfaceGeometry surface(Dune::GeometryTypes::simplex(2),surfaceCorners);
  Dune::Functions::GeometryBindContext<SurfaceGeometry> surfaceContext(surface);
  AffineVectorLocalBasis<2> vectorBasis;
  Dune::FieldVector<double,2> x{0.25,0.25};

  {
    using Transformation = Dune::Functions::ContravariantPiolaTransformation<SurfaceGeometry>;
    Dune::Functions::TransformedLocalBasis<
      AffineVectorLocalBasis<2>,decltype(surfaceContext),Transformation> basis;
    basis.bind(vectorBasis);
    basis.bind(surfaceContext);

    std::vector<decltype(basis)::DerivativeRange<Dune::Functions::Derivatives::Value>> values;
    std::vector<decltype(basis)::DerivativeRange<Dune::Functions::Derivatives::Divergence>> divergences;
    basis.evaluate(Dune::Functions::Derivatives::Value{},x,values);
    basis.evaluate(Dune::Functions::Derivatives::Divergence{},x,divergences);

    test.check(close(values[0][0],1.0/6.0));
    test.check(close(values[0][1],5.0/8.0));
    test.check(close(values[0][2],0.0));
    test.check(close(divergences[0],4.0/6.0));

    Transformation transformation;
    transformation.bind(surfaceContext);
    auto pullback = transformation.localFunctionPullback(
      [](auto const&) { return Dune::FieldVector<double,3>{1,2,7}; });
    auto local = pullback(x);
    test.check(close(local[0],3.0));
    test.check(close(local[1],4.0));
  }

  {
    using Transformation = Dune::Functions::CovariantPiolaTransformation<SurfaceGeometry>;
    Dune::Functions::TransformedLocalBasis<
      AffineVectorLocalBasis<2>,decltype(surfaceContext),Transformation> basis;
    basis.bind(vectorBasis);
    basis.bind(surfaceContext);

    std::vector<decltype(basis)::DerivativeRange<Dune::Functions::Derivatives::Value>> values;
    std::vector<decltype(basis)::DerivativeRange<Dune::Functions::Derivatives::Curl>> curls;
    basis.evaluate(Dune::Functions::Derivatives::Value{},x,values);
    basis.evaluate(Dune::Functions::Derivatives::Curl{},x,curls);

    test.check(close(values[0][0],1.0/4.0));
    test.check(close(values[0][1],5.0/12.0));
    test.check(close(values[0][2],0.0));
    test.check(close(curls[0],1.0/6.0));

    Transformation transformation;
    transformation.bind(surfaceContext);
    auto pullback = transformation.localFunctionPullback(
      [](auto const&) { return Dune::FieldVector<double,3>{1,2,7}; });
    auto local = pullback(x);
    test.check(close(local[0],2.0));
    test.check(close(local[1],6.0));
  }

  using CurveGeometry = Dune::MultiLinearGeometry<double,1,2>;
  std::vector<Dune::FieldVector<double,2>> curveCorners{
    Dune::FieldVector<double,2>{0,0},
    Dune::FieldVector<double,2>{2,2}};
  CurveGeometry curve(Dune::GeometryTypes::simplex(1),curveCorners);
  using CurveContext = Dune::Functions::GeometryBindContext<CurveGeometry>;
  CurveContext curveContext(curve);
  using CurveStage = Dune::Functions::GeometryDerivativePullbackStage<CurveGeometry>;
  static_assert(!Dune::Functions::Concept::TransformationStage<
    CurveStage,
    Dune::Functions::Derivatives::Hessian,
    AffineScalarLocalBasis,
    CurveContext,
    AffineScalarLocalBasis::Traits::JacobianType>);
  static_assert(!Dune::Functions::Concept::LocalBasisTransformation<
    Dune::Functions::CovariantPiolaTransformation<CurveGeometry>,
    AffineScalarLocalBasis,
    CurveContext,
    Dune::Functions::Derivatives::Curl>);
  using ScalarTransformation = Dune::Functions::TransformationPipeline<
    CurveContext,CurveStage>;
  AffineScalarLocalBasis scalarBasis;
  Dune::Functions::TransformedLocalBasis<
    AffineScalarLocalBasis,CurveContext,ScalarTransformation> transformedScalarBasis;
  transformedScalarBasis.bind(scalarBasis);
  transformedScalarBasis.bind(curveContext);

  std::vector<decltype(transformedScalarBasis)::DerivativeRange<
    Dune::Functions::Derivatives::Gradient>> gradients;
  transformedScalarBasis.evaluate(
    Dune::Functions::Derivatives::Gradient{},Dune::FieldVector<double,1>{0.5},gradients);
  test.check(close(gradients[0][0],0.25));
  test.check(close(gradients[0][1],0.25));

  using CustomPipeline = Dune::Functions::BasicTransformationPipeline<
    CustomReferenceEvaluator,CurveContext,PassThroughStage,PassThroughStage>;
  Dune::Functions::TransformedLocalBasis<
    AffineScalarLocalBasis,CurveContext,CustomPipeline> customBasis;
  customBasis.bind(scalarBasis);
  customBasis.bind(curveContext);
  typename decltype(customBasis)::template PrecomputeBuffer<CustomDerivative> customBuffer;
  std::vector<typename decltype(customBasis)::template DerivativeRange<CustomDerivative>> customValues;
  customBasis.precompute(CustomDerivative{},Dune::FieldVector<double,1>{0.5},customBuffer);
  auto referenceValues = customBuffer.reference;
  customBasis.finalize(
    CustomDerivative{},Dune::FieldVector<double,1>{0.5},customBuffer,customValues);
  auto* intermediateData = std::get<0>(customBuffer.intermediate).data();
  customBasis.finalize(
    CustomDerivative{},Dune::FieldVector<double,1>{0.5},customBuffer,customValues);
  test.check(close(customValues[0][0],1.0),
    "custom reference evaluator is used by the pipeline");
  test.check(std::get<0>(customBuffer.intermediate).data() == intermediateData,
    "pipeline intermediate storage is reused");
  test.check(customBuffer.reference == referenceValues,
    "pipeline finalization leaves cached reference values unchanged");

  using NonAffineGeometry = Dune::MultiLinearGeometry<double,2,2>;
  std::vector<Dune::FieldVector<double,2>> nonAffineCorners{
    {0,0},{1,0},{0,1},{1.2,1}};
  NonAffineGeometry nonAffineGeometry(
    Dune::GeometryTypes::cube(2),nonAffineCorners);
  Dune::Functions::GeometryBindContext<NonAffineGeometry> nonAffineContext(
    nonAffineGeometry);
  Dune::Functions::GeometryDerivativePullbackStage<NonAffineGeometry> hessianStage;
  hessianStage.bind(nonAffineContext);
  ScalarLocalBasis2 scalarBasis2;
  std::vector<Dune::FieldMatrix<double,2,2>> referenceHessians(1);
  std::vector<Dune::FieldMatrix<double,2,2>> physicalHessians;
  bool rejectedNonAffineHessian = false;
  try {
    hessianStage.transform(
      Dune::Functions::Derivatives::Hessian{},
      scalarBasis2,
      Dune::FieldVector<double,2>{0.5,0.5},
      referenceHessians,
      physicalHessians);
  }
  catch (Dune::NotImplemented const&) {
    rejectedNonAffineHessian = true;
  }
  test.check(rejectedNonAffineHessian,
    "non-affine physical Hessians are rejected explicitly");

  return test.exit();
}
