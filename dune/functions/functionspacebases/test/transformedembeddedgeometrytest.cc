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

#include <dune/functions/functionspacebases/transformed/bindcontext.hh>
#include <dune/functions/functionspacebases/transformed/geometryderivative.hh>
#include <dune/functions/functionspacebases/transformed/localfiniteelement.hh>
#include <dune/functions/functionspacebases/transformed/pipeline.hh>
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

  void evaluateDivergence(typename Traits::DomainType const&,
                          std::vector<double>& out) const
  {
    divergenceEvaluated = true;
    out = {4};
  }

  void evaluateCurl(typename Traits::DomainType const&,
                    std::vector<double>& out) const
  {
    curlEvaluated = true;
    out = {1};
  }

  mutable bool divergenceEvaluated = false;
  mutable bool curlEvaluated = false;
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
  struct Traits : Dune::LocalBasisTraits<
      double,2,Dune::FieldVector<double,2>,
      double,1,Dune::FieldVector<double,1>,Dune::FieldMatrix<double,1,2>>
  {
    using HessianType = Dune::FieldMatrix<double,2,2>;
  };

  std::size_t size() const { return 1; }
  int order() const { return 2; }

  void evaluateFunction(typename Traits::DomainType const& x,
                        std::vector<typename Traits::RangeType>& out) const
  {
    out.resize(1);
    out[0][0] = x[0]*x[0] + 3*x[1]*x[1];
  }

  void evaluateHessian(typename Traits::DomainType const&,
                       std::vector<typename Traits::HessianType>& out) const
  {
    out.resize(1);
    out[0] = 0;
    out[0][0][0] = 2;
    out[0][1][1] = 6;
  }

  void evaluateLaplacian(typename Traits::DomainType const&,
                         std::vector<typename Traits::RangeType>& out) const
  {
    out.resize(1);
    out[0][0] = 8;
  }
};

struct DummyLocalCoefficients {};
struct DummyLocalInterpolation {};

class ReferenceLocalFiniteElement
{
public:
  using Traits = Dune::LocalFiniteElementTraits<
    AffineScalarLocalBasis,DummyLocalCoefficients,DummyLocalInterpolation>;

  AffineScalarLocalBasis const& localBasis() const { return basis_; }
  DummyLocalCoefficients const& localCoefficients() const { return coefficients_; }
  DummyLocalInterpolation const& localInterpolation() const { return interpolation_; }

private:
  AffineScalarLocalBasis basis_;
  DummyLocalCoefficients coefficients_;
  DummyLocalInterpolation interpolation_;
};

struct InvalidContext {};
struct InvalidInterpolationTransformation {};

template<class LocalFiniteElement, class Context, class Transformation,
         class InterpolationTransformation, Dune::Functions::LocalBasisMode mode>
concept ValidTransformedLocalFiniteElement =
  requires {
    typename Dune::Functions::TransformedLocalFiniteElement<
      LocalFiniteElement,Context,Transformation,InterpolationTransformation,mode>;
  };

bool close(double a, double b)
{
  return std::abs(a-b) < 1e-12;
}

/**
 * \brief Check contravariant surface Piola values, divergence, and interpolation pullback.
 */
void testContravariantSurfacePiola(Dune::TestSuite& test)
{
  using SurfaceGeometry = Dune::MultiLinearGeometry<double,2,3>;
  std::vector<Dune::FieldVector<double,3>> surfaceCorners{
    Dune::FieldVector<double,3>{0,0,0},
    Dune::FieldVector<double,3>{2,0,0},
    Dune::FieldVector<double,3>{0,3,0}};
  SurfaceGeometry surface(Dune::GeometryTypes::simplex(2),surfaceCorners);
  Dune::Functions::GeometryBindContext<SurfaceGeometry> surfaceContext(surface);
  AffineVectorLocalBasis<2> vectorBasis;
  Dune::FieldVector<double,2> x{0.25,0.25};

  using Transformation = Dune::Functions::ContravariantPiolaTransformation<SurfaceGeometry>;
  Dune::Functions::TransformedLocalBasis<
    AffineVectorLocalBasis<2>,decltype(surfaceContext),Transformation> basis;
  basis.bind(vectorBasis);
  basis.bind(surfaceContext);

  std::vector<decltype(basis)::DerivativeRange<Dune::Functions::Derivatives::Value>> values;
  std::vector<decltype(basis)::DerivativeRange<Dune::Functions::Derivatives::Divergence>> divergences;
  using DivergenceBuffer = decltype(basis)::PrecomputeBuffer<
    Dune::Functions::Derivatives::Divergence>;
  static_assert(std::same_as<typename DivergenceBuffer::value_type,double>);
  basis.evaluate(Dune::Functions::Derivatives::Value{},x,values);
  basis.evaluate(Dune::Functions::Derivatives::Divergence{},x,divergences);

  test.check(close(values[0][0],1.0/6.0));
  test.check(close(values[0][1],5.0/8.0));
  test.check(close(values[0][2],0.0));
  test.check(close(divergences[0],4.0/6.0));
  test.check(vectorBasis.divergenceEvaluated,
    "contravariant Piola precomputation uses reference divergence directly");

  Transformation transformation;
  transformation.bind(surfaceContext);
  auto pullback = transformation.localFunctionPullback(
    [](auto const&) { return Dune::FieldVector<double,3>{1,2,7}; });
  auto local = pullback(x);
  test.check(close(local[0],3.0));
  test.check(close(local[1],4.0));
}

/**
 * \brief Check covariant surface Piola values, curl, and interpolation pullback.
 */
void testCovariantSurfacePiola(Dune::TestSuite& test)
{
  using SurfaceGeometry = Dune::MultiLinearGeometry<double,2,3>;
  std::vector<Dune::FieldVector<double,3>> surfaceCorners{
    Dune::FieldVector<double,3>{0,0,0},
    Dune::FieldVector<double,3>{2,0,0},
    Dune::FieldVector<double,3>{0,3,0}};
  SurfaceGeometry surface(Dune::GeometryTypes::simplex(2),surfaceCorners);
  Dune::Functions::GeometryBindContext<SurfaceGeometry> surfaceContext(surface);
  AffineVectorLocalBasis<2> vectorBasis;
  Dune::FieldVector<double,2> x{0.25,0.25};

  using Transformation = Dune::Functions::CovariantPiolaTransformation<SurfaceGeometry>;
  Dune::Functions::TransformedLocalBasis<
    AffineVectorLocalBasis<2>,decltype(surfaceContext),Transformation> basis;
  basis.bind(vectorBasis);
  basis.bind(surfaceContext);

  std::vector<decltype(basis)::DerivativeRange<Dune::Functions::Derivatives::Value>> values;
  std::vector<decltype(basis)::DerivativeRange<Dune::Functions::Derivatives::Curl>> curls;
  using CurlBuffer = decltype(basis)::PrecomputeBuffer<
    Dune::Functions::Derivatives::Curl>;
  static_assert(std::same_as<typename CurlBuffer::value_type,double>);
  basis.evaluate(Dune::Functions::Derivatives::Value{},x,values);
  basis.evaluate(Dune::Functions::Derivatives::Curl{},x,curls);

  test.check(close(values[0][0],1.0/4.0));
  test.check(close(values[0][1],5.0/12.0));
  test.check(close(values[0][2],0.0));
  test.check(close(curls[0],1.0/6.0));
  test.check(vectorBasis.curlEvaluated,
    "covariant Piola precomputation uses reference curl directly");

  Transformation transformation;
  transformation.bind(surfaceContext);
  auto pullback = transformation.localFunctionPullback(
    [](auto const&) { return Dune::FieldVector<double,3>{1,2,7}; });
  auto local = pullback(x);
  test.check(close(local[0],2.0));
  test.check(close(local[1],6.0));
}

/**
 * \brief Check tangential gradients and reusable pipeline storage.
 */
void testScalarGradientAndPipeline(Dune::TestSuite& test)
{
  using CurveGeometry = Dune::MultiLinearGeometry<double,1,2>;
  std::vector<Dune::FieldVector<double,2>> curveCorners{
    Dune::FieldVector<double,2>{0,0},
    Dune::FieldVector<double,2>{2,2}};
  CurveGeometry curve(Dune::GeometryTypes::simplex(1),curveCorners);
  using CurveContext = Dune::Functions::GeometryBindContext<CurveGeometry>;
  CurveContext curveContext(curve);
  using CurveStage = Dune::Functions::GeometryDerivativeStage<CurveGeometry>;
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
  using ScalarTransformation = Dune::Functions::GeometryDerivativePipeline<
    CurveContext,CurveStage>;
  AffineScalarLocalBasis scalarBasis;
  Dune::Functions::TransformedLocalBasis<
    AffineScalarLocalBasis,CurveContext,ScalarTransformation> transformedScalarBasis;
  transformedScalarBasis.bind(scalarBasis);
  transformedScalarBasis.bind(curveContext);

  std::vector<decltype(transformedScalarBasis)::DerivativeRange<
    Dune::Functions::Derivatives::Gradient>> gradients;
  typename decltype(transformedScalarBasis)::PrecomputeBuffer<
    Dune::Functions::Derivatives::Gradient> gradientBuffer;
  transformedScalarBasis.precompute(
    Dune::Functions::Derivatives::Gradient{},
    Dune::FieldVector<double,1>{0.5},
    gradientBuffer);
  static_assert(std::same_as<
    typename decltype(gradientBuffer.reference)::value_type,
    Dune::FieldVector<double,1>>);
  transformedScalarBasis.evaluate(
    Dune::Functions::Derivatives::Gradient{},Dune::FieldVector<double,1>{0.5},gradients);
  test.check(close(gradients[0][0],0.25));
  test.check(close(gradients[0][1],0.25));

  using CustomPipeline = Dune::Functions::BasicBasisEvaluationPipeline<
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
}

/**
 * \brief Check semantic reference Laplacians and Hessian-based physical Laplacians.
 */
void testReferenceOperatorSelection(Dune::TestSuite& test)
{
  ScalarLocalBasis2 localBasis;
  Dune::FieldVector<double,2> x{0.25,0.25};

  std::vector<Dune::FieldVector<double,1>> referenceLaplacians;
  Dune::Functions::ReferenceLocalBasisEvaluator{}.evaluate(
    Dune::Functions::Derivatives::Laplacian{},localBasis,x,referenceLaplacians);
  test.check(close(referenceLaplacians[0][0],8),
    "reference evaluator computes the named reference Laplacian");

  using Geometry = Dune::MultiLinearGeometry<double,2,2>;
  std::vector<Dune::FieldVector<double,2>> corners{
    {0,0},{2,0},{0,3}};
  Geometry geometry(Dune::GeometryTypes::simplex(2),corners);
  using Context = Dune::Functions::GeometryBindContext<Geometry>;
  Context context(geometry);
  using Transformation = Dune::Functions::GeometryDerivativePipeline<
    Context,Dune::Functions::GeometryDerivativeStage<Geometry>>;
  Dune::Functions::TransformedLocalBasis<
    ScalarLocalBasis2,Context,Transformation> basis;
  basis.bind(localBasis);
  basis.bind(context);

  using LaplacianBuffer = decltype(basis)::PrecomputeBuffer<
    Dune::Functions::Derivatives::Laplacian>;
  static_assert(std::same_as<
    typename decltype(std::declval<LaplacianBuffer>().reference)::value_type,
    Dune::FieldMatrix<double,2,2>>);

  std::vector<decltype(basis)::DerivativeRange<
    Dune::Functions::Derivatives::Laplacian>> physicalLaplacians;
  basis.evaluate(Dune::Functions::Derivatives::Laplacian{},x,physicalLaplacians);
  test.check(close(physicalLaplacians[0][0],7.0/6.0),
    "geometry pipeline obtains the physical Laplacian from the reference Hessian");
}

/**
 * \brief Check compile-time policy diagnostics and non-affine Hessian rejection.
 */
void testPolicyConstraintsAndHessianRejection(Dune::TestSuite& test)
{
  using NonAffineGeometry = Dune::MultiLinearGeometry<double,2,2>;
  std::vector<Dune::FieldVector<double,2>> nonAffineCorners{
    {0,0},{1,0},{0,1},{1.2,1}};
  NonAffineGeometry nonAffineGeometry(
    Dune::GeometryTypes::cube(2),nonAffineCorners);
  Dune::Functions::GeometryBindContext<NonAffineGeometry> nonAffineContext(
    nonAffineGeometry);
  Dune::Functions::GeometryDerivativeTransformation<NonAffineGeometry> hessianTransformation;
  hessianTransformation.bind(nonAffineContext);
  ScalarLocalBasis2 scalarBasis2;
  Dune::FieldMatrix<double,2,2> referenceHessian;
  Dune::FieldMatrix<double,2,2> physicalHessian;
  bool rejectedNonAffineHessian = false;
  try {
    hessianTransformation.transform(
      Dune::Functions::Derivatives::Hessian{},
      scalarBasis2,
      Dune::FieldVector<double,2>{0.5,0.5},
      referenceHessian,
      physicalHessian);
  }
  catch (Dune::NotImplemented const&) {
    rejectedNonAffineHessian = true;
  }
  test.check(rejectedNonAffineHessian,
    "non-affine physical Hessians are rejected explicitly");

  using Context = Dune::Functions::GeometryBindContext<NonAffineGeometry>;
  using Transformation = Dune::Functions::GeometryDerivativePipeline<
    Context,Dune::Functions::GeometryDerivativeStage<NonAffineGeometry>>;
  static_assert(Dune::Functions::Concept::ReferenceLocalFiniteElement<
    ReferenceLocalFiniteElement>);
  static_assert(ValidTransformedLocalFiniteElement<
    ReferenceLocalFiniteElement,
    Context,
    Transformation,
    Dune::Functions::NoInterpolationTransformation,
    Dune::Functions::LocalBasisMode::physical>);
  static_assert(!ValidTransformedLocalFiniteElement<
    ReferenceLocalFiniteElement,
    InvalidContext,
    Transformation,
    Dune::Functions::NoInterpolationTransformation,
    Dune::Functions::LocalBasisMode::physical>);
  static_assert(!ValidTransformedLocalFiniteElement<
    ReferenceLocalFiniteElement,
    Context,
    Transformation,
    InvalidInterpolationTransformation,
    Dune::Functions::LocalBasisMode::physical>);
}

} // namespace

int main()
{
  Dune::TestSuite test;

  testContravariantSurfacePiola(test);
  testCovariantSurfacePiola(test);
  testScalarGradientAndPipeline(test);
  testReferenceOperatorSelection(test);
  testPolicyConstraintsAndHessianRejection(test);

  return test.exit();
}
