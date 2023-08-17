// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BOUND_TEST_LOCALFE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_TEST_BOUND_TEST_LOCALFE_HH

/** \file \brief Unit tests for GlobalValuedLocalFiniteElement objects that utilize
 * global(!)derivative degrees of freedom
 *
 *  \note This file is an adapted copy from
 * dune-localfunctions/dune/localfunctions/test/test-localfe.hh The tests in this file do not
 * replace the tests in the file mentioned above, but apply them on the local finite element bound
 * to a grid element. For the same reason this file is part of dune-functions rather than dune-
 * localfunctions.
 *
 *  \note This header is not part of the official Dune API and might be subject to
 * change.  You can use this header to test external finite element implementations, but be warned
 * that your tests might break with future Dune versions.
 */

#include <iomanip>
#include <iostream>
#include <typeinfo>

#include <dune/common/classname.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/common/transpose.hh>
#include <dune/localfunctions/common/derivative.hh>
#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

double TOL = 1e-9;
// The FD approximation used for checking the Jacobian uses half of the
// precision -- so we have to be a little bit more tolerant here.
double jacobianTOL = 1e-5; // sqrt(TOL)

// Shapefunctions need to have a derivative
template <class LFE, class Element>
class ShapeFunctionDerivativeAsCallable;

// This class wraps one shape function of a local finite element as a callable
// that can be fed to the LocalInterpolation::interpolate method.
template <class FE, class Element>
class ShapeFunctionAsCallable
{
public:
  typedef typename FE::Traits::LocalBasisType::Traits::DomainType DomainType;
  typedef typename FE::Traits::LocalBasisType::Traits::RangeType RangeType;
  typedef typename FE::Traits::LocalBasisType::Traits::RangeFieldType CT;

  ShapeFunctionAsCallable(const FE &fe, int shapeFunction, Element const &element)
      : fe_(fe), shapeFunction_(shapeFunction), e_(element)
  {
  }

  RangeType operator()(DomainType x) const
  {
    std::vector<RangeType> yy;
    fe_.localBasis().evaluateFunction(x, yy);
    return yy[shapeFunction_];
  }

  friend ShapeFunctionDerivativeAsCallable<FE, Element> derivative(ShapeFunctionAsCallable const &t)
  {
    return ShapeFunctionDerivativeAsCallable(t.fe_, t.shapeFunction_, t.e_);
  }

private:
  const FE &fe_;
  int shapeFunction_;
  Element const &e_;
};

template <class LFE, class Element>
class ShapeFunctionDerivativeAsCallable
{
  typedef typename LFE::Traits::LocalBasisType::Traits::DomainType DomainType;
  typedef typename LFE::Traits::LocalBasisType::Traits::RangeType ValueType;
  typedef typename LFE::Traits::LocalBasisType::Traits::RangeFieldType R;
  static int constexpr dimRange = LFE::Traits::LocalBasisType::Traits::dimRange;

  typedef typename LFE::Traits::LocalBasisType::Traits::JacobianType JacobianType;
  static constexpr unsigned int dim = LFE::Traits::LocalBasisType::Traits::dimDomain;

public:
  ShapeFunctionDerivativeAsCallable(const LFE &fe, int shapeFunction, Element const &element)
      : fe_(fe), shapeFunction_(shapeFunction), e_(element)
  {
  }

  JacobianType operator()(DomainType const &x) const
  {
    std::vector<JacobianType> yy;
    fe_.localBasis().evaluateJacobian(x, yy);
    // JacobianType returnValue = 0.;
    // e_.geometry().jacobianInverseTransposed(x).mv(yy[shapeFunction_][0], returnValue[0]);
    // return returnValue;
    return yy[shapeFunction_] * transpose(e_.geometry().jacobianInverseTransposed(x));
  }

  auto friend derivative(ShapeFunctionDerivativeAsCallable const &t)
  {
    return [&t](DomainType const &x)
    {
      Dune::FieldVector<Dune::FieldMatrix<R, dim, dim>, dimRange> referenceHessian, hessian;
      std::vector<ValueType> partials;
      std::array<unsigned int, dim> dir;
      for (std::size_t i = 0; i < dim; ++i)
        for (std::size_t j = 0; j < dim; ++j)
        {
          dir = {};
          dir[i]++;
          dir[j]++;
          t.fe_.localBasis().partial(dir, x, partials);
          for (std::size_t k = 0; k < dimRange; ++k)
            referenceHessian[k][i][j] = partials[t.shapeFunction_][k];
        }
      hessian = 0;
      auto JIT = t.e_.geometry().jacobianInverseTransposed(x);
      for (std::size_t k = 0; k < dimRange; ++k)
        hessian[k] = JIT * referenceHessian[k] * transpose(JIT);
      return hessian;
    };
  }

public:
  const LFE &fe_;
  int shapeFunction_;
  Element const &e_;
};
// Check whether the degrees of freedom computed by LocalInterpolation
// are dual to the shape functions.  See Ciarlet, "The Finite Element Method
// for Elliptic Problems", 1978, for details.
template <class FE, class Element>
bool testLocalInterpolation(const FE &fe, Element const &element)
{
  bool ret = true;
  std::vector<typename ShapeFunctionAsCallable<FE, Element>::CT> coeff;
  for (size_t i = 0; i < fe.size(); ++i)
  {
    //////////////////////////////////////////////////////////////////////////////
    //  Part A: Feed the shape functions to the 'interpolate' method in form of
    //    a class providing an evaluate() method.
    //    This way is deprecated since dune-localfunctions 2.7 and hence removed.
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    //  Part B: Redo the same test, but feed the shape functions to the
    //    'interpolate' method in form of a callable.
    //////////////////////////////////////////////////////////////////////////////

    // The i-th shape function as a function that 'interpolate' can deal with
    ShapeFunctionAsCallable<FE, Element> sfAsCallable(fe, i, element);

    // Compute degrees of freedom for that shape function
    // We expect the result to be the i-th unit vector
    fe.localInterpolation().interpolate(sfAsCallable, coeff);

    // Check size of weight vector
    if (coeff.size() != fe.localBasis().size())
    {
      std::cout << "Bug in LocalInterpolation for finite element type " << Dune::className(fe)
                << std::endl;
      std::cout << "    Interpolation produces " << coeff.size() << " degrees of freedom"
                << std::endl;
      std::cout << "    Basis has size " << fe.localBasis().size() << std::endl;
      std::cout << std::endl;
      return false;
    }

    // Check if interpolation weights are equal to coefficients
    for (std::size_t j = 0; j < coeff.size(); ++j)
    {
      if (std::abs(coeff[j] - (i == j)) > TOL)
      {
        std::cout << std::setprecision(16);
        std::cout << "Bug in LocalInterpolation for finite element type " << Dune::className(fe)
                  << std::endl;
        std::cout << "    Degree of freedom " << j << " applied to shape function " << i
                  << " yields value " << coeff[j] << ", not the expected value " << (i == j)
                  << std::endl;
        std::cout << std::endl;
        ret = false;
      }
    }
  }
  return ret;
}

template <class Domain, class Range>
class Constant
{
public:
  Constant(double val) : value(val) {}

  Range operator()(Domain const &x) const { return value; }

  friend Constant<Domain, typename Dune::Impl::DerivativeTraits<Range(Domain)>::type>
  derivative(Constant const &t)
  {
    return {0.};
  }

private:
  Range value;
};

template<class K, int d>
auto norm(Dune::FieldVector<K, d> v) -> typename Dune::FieldTraits<Dune::FieldVector<K, d>>::real_type
{
  return v.two_norm();
}

template<class K, int n, int m>
auto norm(Dune::FieldMatrix<K, n, m> mat) -> typename Dune::FieldTraits<Dune::FieldMatrix<K, n, m>>::real_type
{
  return mat.frobenius_norm();
}

template<class K>
auto norm(K scalar) {return std::abs(scalar);}

template<class K, int n>
auto getDerivativeComponent(Dune::FieldVector<K, n> v, std::size_t i)->K{return v[i];}

template<class K, int n, int m>
auto getDerivativeComponent(Dune::FieldMatrix<K, n, m> mat, std::size_t i) {
  Dune::FieldVector<K, n> ret;
  for (std::size_t j = 0; j < n; ++j)
    ret[j] = mat[j][i];
  return ret;
 }

template<class K, int n, int m, int l>
auto getDerivativeComponent(Dune::FieldMatrix<Dune::FieldVector<K,l>, n, m> mat, std::size_t i)
{
  Dune::FieldMatrix<K, n, m> ret;
  for (std::size_t j = 0; j < n; ++j)
    for (std::size_t k = 0; k < m; ++k)
      ret[j][k] = mat[j][k][i];
  return ret;
}


// Check whether the space spanned by the shape functions
// contains the constant functions
template <class FE>
bool testCanRepresentConstants(const FE &fe, unsigned order = 5)
{
  typedef typename FE::Traits::LocalBasisType LB;
  using RangeType = typename LB::Traits::RangeType;
  using JacobianType = typename LB::Traits::JacobianType;
  bool success = true;

  // Construct the constant '1' function
  // auto constantOne = [](const typename LB::Traits::DomainType& xi) { return RangeType(1.0); };
  Constant<typename LB::Traits::DomainType, RangeType> constantOne(1.);
  // Project the constant function onto the FE space
  std::vector<double> coefficients;
  fe.localInterpolation().interpolate(constantOne, coefficients);

  // A set of test points
  const auto &quad = Dune::QuadratureRules<double, LB::Traits::dimDomain>::rule(fe.type(), order);

  // Loop over all quadrature points
  for (size_t i = 0; i < quad.size(); i++)
  {

    // Get a test point
    const auto &testPoint = quad[i].position();

    // Compute value of the representation of constantOne at the test point
    std::vector<RangeType> values;
    fe.localBasis().evaluateFunction(testPoint, values);

    RangeType sum(0);
    for (size_t j = 0; j < values.size(); j++)
      sum += coefficients[j] * values[j];

    if (norm(RangeType(1.0)-sum) > TOL)
    {
      std::cout << "Finite element type " << Dune::className(fe)
                << " cannot represent constant functions!" << std::endl;
      std::cout << "    At position: " << testPoint << "," << std::endl;
      std::cout << "    discrete approximation of the '1' function has value " << sum << std::endl;
      std::cout << std::endl;
      success = false;
    }

  } // Loop over all quadrature points

  return success;
}

// check whether Jacobian agrees with FD approximation
template <class FE>
bool testJacobian(
    const FE &fe, unsigned order = 2,
    const std::function<bool(const typename FE::Traits::LocalBasisType::Traits::DomainType &)>
        derivativePointSkip
    = nullptr)
{
  typedef typename FE::Traits::LocalBasisType LB;

  bool success = true;

  // ////////////////////////////////////////////////////////////
  //   Check the partial derivatives by comparing them
  //   to finite difference approximations
  // ////////////////////////////////////////////////////////////

  // A set of test points
  const Dune::QuadratureRule<double, LB::Traits::dimDomain> quad
      = Dune::QuadratureRules<double, LB::Traits::dimDomain>::rule(fe.type(), order);

  // Loop over all quadrature points
  for (size_t i = 0; i < quad.size(); i++)
  {

    // Get a test point
    const Dune::FieldVector<double, LB::Traits::dimDomain> &testPoint = quad[i].position();

    // Get the shape function derivatives there
    std::vector<typename LB::Traits::JacobianType> jacobians;
    fe.localBasis().evaluateJacobian(testPoint, jacobians);
    if (jacobians.size() != fe.localBasis().size())
    {
      std::cout << "Bug in evaluateJacobian() for finite element type " << Dune::className(fe)
                << std::endl;
      std::cout << "    Jacobian vector has size " << jacobians.size() << std::endl;
      std::cout << "    Basis has size " << fe.localBasis().size() << std::endl;
      std::cout << std::endl;
      return false;
    }

    // Skip this test point if we are supposed to
    if (derivativePointSkip && derivativePointSkip(quad[i].position()))
      continue;

    // Loop over all directions
    for (int k = 0; k < LB::Traits::dimDomain; k++)
    {

      // Compute an approximation to the derivative by finite differences
      Dune::FieldVector<double, LB::Traits::dimDomain> upPos = testPoint;
      Dune::FieldVector<double, LB::Traits::dimDomain> downPos = testPoint;

      upPos[k] += jacobianTOL;
      downPos[k] -= jacobianTOL;

      std::vector<typename LB::Traits::RangeType> upValues, downValues;

      fe.localBasis().evaluateFunction(upPos, upValues);
      fe.localBasis().evaluateFunction(downPos, downValues);

      // Loop over all shape functions in this set
      for (unsigned int j = 0; j < fe.localBasis().size(); ++j)
      {
          // The current partial derivative, just for ease of notation
        auto derivative = getDerivativeComponent(jacobians[j],k);

        auto finiteDiff = (upValues[j] - downValues[j])
                            / (2*jacobianTOL);

          // Check
        if ( norm(derivative-finiteDiff) >
              TOL/jacobianTOL*((norm(finiteDiff)>1) ? norm(finiteDiff) : 1.) )
        {
          std::cout << std::setprecision(16);
          std::cout << "Bug in evaluateJacobian() for finite element type "
                    << Dune::className(fe) << std::endl;
          std::cout << "    Shape function derivative does not agree with "
                    << "FD approximation" << std::endl;
          std::cout << "    Shape function " << j
                    << " at position " << testPoint << ": derivative in "
                    << "direction " << k << " is " << derivative << ", but "
                    << finiteDiff << " is expected." << std::endl;
          std::cout << std::endl;
          success = false;
        }
      } // Loop over all shape functions in this set
    } // Loop over all directions
  } // Loop over all quadrature points

  return success;
}

/** \brief Helper class to test the 'partial' method
 *
 * It implements a static loop over the available diff orders
 */
struct TestPartial
{
  template <class FE, unsigned int diffOrder>
  static bool test(const FE& fe,
                   double eps, double delta, std::size_t order = 2,
                  const std::function<bool(const typename FE::Traits::LocalBasisType::Traits::DomainType&)> derivativePointSkip = nullptr)
  {
    bool success = true;

    if constexpr(diffOrder > 2)
      std::cout << "No test for differentiability orders larger than 2!" << std::endl;

    if constexpr (diffOrder >= 2)
      success = success and testOrder2(fe, eps, delta, order, derivativePointSkip);

    if constexpr (diffOrder >= 1)
      success = success and testOrder1(fe, eps, delta, order, derivativePointSkip);

    success = success and testOrder0(fe, eps, delta, order);

    return success;
  }

  /** \brief Test the 'partial' method for zero-order partial derivatives, i.e., values */
  template <class FE>
  static bool testOrder0(const FE& fe,
                   double eps,
                   double delta,
                   std::size_t order = 2)
  {
    typedef typename FE::Traits::LocalBasisType::Traits::RangeType RangeType;
    constexpr auto dimDomain = FE::Traits::LocalBasisType::Traits::dimDomain;

    bool success = true;

    //////////////////////////////////////////////////////////////
    //   Check the partial derivatives by comparing them
    //   to finite difference approximations
    //////////////////////////////////////////////////////////////

    // A set of test points
    const auto& quad = Dune::QuadratureRules<double, dimDomain>::rule(fe.type(),
                                                                      order);

    // Loop over all quadrature points
    for (size_t i = 0; i < quad.size(); i++)
    {
      // Get a test point
      const Dune::FieldVector<double, dimDomain>& testPoint = quad[i].position();

      // Get the shape function values there using the 'partial' method
      std::vector<RangeType> partialValues;
      std::array<unsigned int, dimDomain> multiIndex;
      std::fill(multiIndex.begin(), multiIndex.end(), 0);

      fe.localBasis().partial(multiIndex, testPoint, partialValues);

      if (partialValues.size() != fe.localBasis().size())
      {
        std::cout << "Bug in partial() for finite element type "
                  << Dune::className(fe) << std::endl;
        std::cout << "    values vector has size "
                  << partialValues.size() << std::endl;
        std::cout << "    Basis has size " << fe.localBasis().size()
                  << std::endl;
        std::cout << std::endl;
        return false;
      }

      // Get reference values
      std::vector<RangeType> referenceValues;
      fe.localBasis().evaluateFunction(testPoint, referenceValues);

      // Loop over all shape functions in this set
      for (unsigned int j = 0; j < fe.localBasis().size(); ++j)
      {
        // Check the 'partial' method
        if (norm(partialValues[j] - referenceValues[j])
            > TOL / jacobianTOL
              * ((norm(referenceValues[j]) > 1) ? norm(referenceValues[j]) : 1.))
        {
          std::cout << std::setprecision(16);
          std::cout << "Bug in partial() for finite element type "
                    << Dune::className(fe) << std::endl;
          std::cout << "    Shape function value does not agree with "
                    << "output of method evaluateFunction." << std::endl;
          std::cout << "    Shape function " << j
                    << " at position " << testPoint << ": value is " << partialValues[j]
                    << ", but " << referenceValues[j] << " is expected." << std::endl;
          std::cout << std::endl;
          success = false;
        }

      } // Loop over all shape functions in this set
    } // Loop over all quadrature points

    return success;
  }

  /** \brief Test the 'partial' method for first-order partial derivatives */
  template <class FE>
  static bool testOrder1(const FE& fe,
                   double eps,
                   double delta,
                   std::size_t order = 2,
                  const std::function<bool(const typename FE::Traits::LocalBasisType::Traits::DomainType&)> derivativePointSkip = nullptr)
  {
    typedef typename FE::Traits::LocalBasisType LB;
    typedef typename LB::Traits::RangeFieldType RangeField;

    bool success = true;

    //////////////////////////////////////////////////////////////
    //   Check the partial derivatives by comparing them
    //   to finite difference approximations
    //////////////////////////////////////////////////////////////

    // A set of test points
    const Dune::QuadratureRule<double, LB::Traits::dimDomain> quad =
          Dune::QuadratureRules<double, LB::Traits::dimDomain>::rule(fe.type(),
                                                                     order);

    // Loop over all quadrature points
    for (size_t i = 0; i < quad.size(); i++)
    {
      // Get a test point
      const Dune::FieldVector<double, LB::Traits::dimDomain>& testPoint = quad[i].position();

      // Skip the test points we are supposed to skip
      if (derivativePointSkip && derivativePointSkip(testPoint))
        continue;

      // Loop over all directions
      for (int k = 0; k < LB::Traits::dimDomain; k++)
      {
        // Get the shape function derivatives there using the 'partial' method
        std::vector<typename LB::Traits::RangeType> firstPartialDerivatives;
        std::array<unsigned int, LB::Traits::dimDomain> multiIndex;
        std::fill(multiIndex.begin(), multiIndex.end(), 0);
        multiIndex[k]++;
        fe.localBasis().partial(multiIndex, testPoint, firstPartialDerivatives);
        if (firstPartialDerivatives.size() != fe.localBasis().size())
        {
          std::cout << "Bug in partial() for finite element type "
                    << Dune::className(fe) << std::endl;
          std::cout << "    firstPartialDerivatives vector has size "
                    << firstPartialDerivatives.size() << std::endl;
          std::cout << "    Basis has size " << fe.localBasis().size()
                    << std::endl;
          std::cout << std::endl;
          return false;
        }

        // Compute an approximation to the derivative by finite differences
        Dune::FieldVector<double, LB::Traits::dimDomain> upPos = testPoint;
        Dune::FieldVector<double, LB::Traits::dimDomain> downPos = testPoint;

        upPos[k] += jacobianTOL;
        downPos[k] -= jacobianTOL;

        std::vector<typename LB::Traits::RangeType> upValues, downValues;

        fe.localBasis().evaluateFunction(upPos, upValues);
        fe.localBasis().evaluateFunction(downPos, downValues);

        // Loop over all shape functions in this set
        for (unsigned int j = 0; j < fe.localBasis().size(); ++j)
        {

            auto finiteDiff = (upValues[j] - downValues[j])
                              / (2 * jacobianTOL);

            // Check the 'partial' method
            auto partialDerivative = firstPartialDerivatives[j];
            if (norm(partialDerivative - finiteDiff)
                > TOL / jacobianTOL
                  * ((norm(finiteDiff) > 1) ? norm(finiteDiff) : 1.))
            {
              std::cout << std::setprecision(16);
              std::cout << "Bug in partial() for finite element type "
                        << Dune::className(fe) << std::endl;
              std::cout << "    Shape function derivative does not agree with "
                        << "FD approximation" << std::endl;
              std::cout << "    Shape function " << j
                        << " at position " << testPoint << ": derivative in "
                        << "direction " << k << " is " << partialDerivative << ", but "
                        << finiteDiff << " is expected." << std::endl;
              std::cout << std::endl;
              success = false;
            }

        } // Loop over all shape functions in this set
      } // Loop over all directions
    } // Loop over all quadrature points

    return success;
  }

  /** \brief Test second-order partial derivatives */
  template <class FE>
  static bool testOrder2(const FE& fe,
                   double eps,
                   double delta,
                   std::size_t order = 2,
                   const std::function<bool(const typename FE::Traits::LocalBasisType::Traits::DomainType&)> derivativePointSkip = nullptr)
  {
    typedef typename FE::Traits::LocalBasisType LocalBasis;
    typedef typename LocalBasis::Traits::DomainFieldType DF;
    typedef typename LocalBasis::Traits::DomainType Domain;
    static const int dimDomain = LocalBasis::Traits::dimDomain;

    static const std::size_t dimR = LocalBasis::Traits::dimRange;
    typedef typename LocalBasis::Traits::RangeType Range;
    typedef typename LocalBasis::Traits::RangeFieldType RangeField;

    bool success = true;

    //////////////////////////////////////////////////////////////
    //   Check the partial derivatives by comparing them
    //   to finite difference approximations
    //////////////////////////////////////////////////////////////

    // A set of test points
    const Dune::QuadratureRule<DF, dimDomain> quad
       = Dune::QuadratureRules<DF,dimDomain>::rule(fe.type(), order);

    // Loop over all quadrature points
    for (std::size_t i = 0; i < quad.size(); i++)
    {
      // Get a test point
      const Domain& testPoint = quad[i].position();

      // Skip the test points we are supposed to skip
      if (derivativePointSkip && derivativePointSkip(testPoint))
        continue;

      // For testing the 'partial' method
      std::array<std::vector<Dune::FieldMatrix<RangeField, dimDomain, dimDomain> >, dimR> partialHessians;
      for (size_t k = 0; k < dimR; k++)
        partialHessians[k].resize(fe.size());

      //loop over all local directions
      for (int dir0 = 0; dir0 < dimDomain; dir0++)
      {
        for (int dir1 = 0; dir1 < dimDomain; dir1++)
        {
          // Get the shape function derivatives there using the 'partial' method
          std::vector<Range> secondPartialDerivative;
          std::array<unsigned int,dimDomain> multiIndex;
          std::fill(multiIndex.begin(), multiIndex.end(), 0);
          multiIndex[dir0]++;
          multiIndex[dir1]++;
          secondPartialDerivative.resize(fe.size());
          fe.localBasis().partial(multiIndex, testPoint, secondPartialDerivative);

          if (secondPartialDerivative.size() != fe.localBasis().size())
          {
            std::cout << "Bug in partial() for finite element type "
                      << Dune::className<FE>() << ":" << std::endl;
            std::cout << "    return vector has size " << secondPartialDerivative.size()
                      << std::endl;
            std::cout << "    Basis has size " << fe.localBasis().size()
                      << std::endl;
            std::cout << std::endl;
            return false;
          }

          //combine to Hesse matrices
          for (size_t k = 0; k < dimR; k++)
            for (std::size_t j = 0; j < fe.localBasis().size(); ++j)
              partialHessians[k][j][dir0][dir1] = secondPartialDerivative[j][k];

        }
      }  //loop over all directions

      // Loop over all local directions
      for (std::size_t dir0 = 0; dir0 < dimDomain; ++dir0)
      {
        for (unsigned int dir1 = 0; dir1 < dimDomain; dir1++)
        {
          // Compute an approximation to the derivative by finite differences
          std::array<Domain,4> neighbourPos;
          std::fill(neighbourPos.begin(), neighbourPos.end(), testPoint);

          neighbourPos[0][dir0] += delta;
          neighbourPos[0][dir1] += delta;
          neighbourPos[1][dir0] -= delta;
          neighbourPos[1][dir1] += delta;
          neighbourPos[2][dir0] += delta;
          neighbourPos[2][dir1] -= delta;
          neighbourPos[3][dir0] -= delta;
          neighbourPos[3][dir1] -= delta;

          std::array<std::vector<Range>, 4> neighbourValues;
          for (int k = 0; k < 4; k++)
            fe.localBasis().evaluateFunction(neighbourPos[k],
                                             neighbourValues[k]);

          // Loop over all shape functions in this set
          for (std::size_t j = 0; j < fe.localBasis().size(); ++j)
          {
            //Loop over all components
            for (std::size_t k = 0; k < dimR; ++k)
            {
              RangeField finiteDiff = (neighbourValues[0][j][k]
                  - neighbourValues[1][j][k] - neighbourValues[2][j][k]
                  + neighbourValues[3][j][k]) / (4 * delta * delta);

              // The current partial derivative, just for ease of notation, evaluated by the 'partial' method
              RangeField partialDerivative = partialHessians[k][j][dir0][dir1];

              // Check
              if (std::abs(partialDerivative - finiteDiff)
                  > eps / delta * (std::max(std::abs(finiteDiff), 1.0)))
              {
                std::cout << std::setprecision(16);
                std::cout << "Bug in partial() for finite element type "
                          << Dune::className<FE>() << ":" << std::endl;
                std::cout << "    Second shape function derivative does not agree with "
                          << "FD approximation" << std::endl;
                std::cout << "    Shape function " << j << " component " << k
                          << " at position " << testPoint << ": derivative in "
                          << "local direction (" << dir0 << ", " << dir1 << ") is "
                          << partialDerivative << ", but " << finiteDiff
                          << " is expected." << std::endl;
                std::cout << std::endl;
                success = false;
              }

            } //Loop over all components
          }
        } // Loop over all local directions
      } // Loop over all shape functions in this set
    } // Loop over all quadrature points

    return success;
  }

};

// Flags for disabling parts of testFE
enum {
  DisableNone = 0,
  DisableLocalInterpolation = 1,
  DisableVirtualInterface = 2,
  DisableJacobian = 4,
  DisableEvaluate = 8,
  DisableRepresentConstants = 16
};

/** \brief Call tests for given finite element
 *
 * \param derivativePointSkip This is a small predicate class that allows to skip certain
 *   points when testing the derivative implementations.  It exists because some
 *   finite elements are not everywhere differentiable, but we still want to run
 *   the tests for derivatives.  Rather than constructing special sets of test
 *   points that avoid the problematic parts of the domain, we simply skip
 *   all test points that happen to be somewhere where the shape functions are
 *   not differentiable.
 */
template <unsigned int diffOrder = 0, class FE, class Element>
bool testFE(
    const FE &fe, Element const &element, char disabledTests = DisableNone,
    const std::function<bool(const typename FE::Traits::LocalBasisType::Traits::DomainType &)>
        derivativePointSkip
    = nullptr)
{
  // Order of the quadrature rule used to generate test points
  unsigned int quadOrder = 2;

  bool success = true;

  if (FE::Traits::LocalBasisType::Traits::dimDomain != fe.type().dim())
  {
    std::cout << "Bug in type() for finite element type " << Dune::className(fe) << std::endl;
    std::cout << "    Coordinate dimension is " << FE::Traits::LocalBasisType::Traits::dimDomain
              << std::endl;
    std::cout << "    but GeometryType is " << fe.type() << " with dimension " << fe.type().dim()
              << std::endl;
    success = false;
  }

  if (fe.size() != fe.localBasis().size())
  {
    std::cout << "Bug in finite element type " << Dune::className(fe) << std::endl;
    std::cout << "    Size reported by LocalFiniteElement is " << fe.size() << std::endl;
    std::cout << "    but size reported by LocalBasis is " << fe.localBasis().size() << std::endl;
    success = false;
  }

  // Make sure evaluateFunction returns the correct number of values
  std::vector<typename FE::Traits::LocalBasisType::Traits::RangeType> values;
  fe.localBasis().evaluateFunction(
      Dune::ReferenceElements<double, FE::Traits::LocalBasisType::Traits::dimDomain>::general(
          fe.type())
          .position(0, 0),
      values);

  if (values.size() != fe.size())
  {
    std::cout << "Bug in finite element type " << Dune::className(fe) << std::endl;
    std::cout << "    LocalFiniteElement.size() returns " << fe.size() << "," << std::endl;
    std::cout << "    but LocalBasis::evaluateFunction returns " << values.size() << " values!"
              << std::endl;
    success = false;
  }

  if (fe.size() != fe.localCoefficients().size())
  {
    std::cout << "Bug in finite element type " << Dune::className(fe) << std::endl;
    std::cout << "    Size reported by LocalFiniteElement is " << fe.size() << std::endl;
    std::cout << "    but size reported by LocalCoefficients is " << fe.localCoefficients().size()
              << std::endl;
    success = false;
  }

  const auto &lc = fe.localCoefficients();
  for (size_t i = 0; i < lc.size(); i++)
  {
    const auto &lk = lc.localKey(i);
    if (lk.codim() > fe.type().dim())
    {
      std::cout << "Bug in finite element type " << Dune::className(fe) << std::endl;
      std::cout << "    Codimension reported by localKey(" << i << ") is " << lk.codim()
                << std::endl;
      std::cout << "    but geometry is " << fe.type() << " with dimension " << fe.type().dim()
                << std::endl;
      success = false;
    }
  }

  if (not(disabledTests & DisableLocalInterpolation))
  {
    success = testLocalInterpolation<FE>(fe, element) and success;
  }

  if (not (disabledTests & DisableRepresentConstants))
  {
    success = testCanRepresentConstants<FE>(fe) and success;
  }

  if (not(disabledTests & DisableJacobian))
  {
    success = testJacobian<FE>(fe, quadOrder, derivativePointSkip) and success;
  }
  else
  {
    // make sure diffOrder is 0
    success = (diffOrder == 0) and success;
  }

  if (not(disabledTests & DisableEvaluate))
  {
    success = TestPartial::test<FE,diffOrder>(fe, TOL, jacobianTOL, quadOrder, derivativePointSkip) and success;
  }

#if 0 // virtual interface currently interpolates into a lambda functions
  if (not(disabledTests & DisableVirtualInterface))
  {
    typedef typename FE::Traits::LocalBasisType::Traits ImplementationLBTraits;
    typedef typename Dune::LocalFiniteElementVirtualInterface<ImplementationLBTraits>
        VirtualFEInterface;
    typedef typename Dune::LocalFiniteElementVirtualImp<FE> VirtualFEImp;

    const VirtualFEImp virtualFE(fe);
    if (not(disabledTests & DisableLocalInterpolation))
      success = testLocalInterpolation<VirtualFEInterface, Element>(virtualFE, element) and success;
    if (not(disabledTests & DisableJacobian))
    {
      success
          = testJacobian<VirtualFEInterface>(virtualFE, quadOrder, derivativePointSkip) and success;
    }
    else
    {
      // make sure diffOrder is 0
      success = (diffOrder == 0) and success;
    }
  }
#endif
  return success;
}

#define TEST_FE(A, element)                                                                        \
  {                                                                                                \
    bool b = testFE(A, element);                                                                   \
    std::cout << "testFE(" #A ") " << (b ? "succeeded\n" : "failed\n");                            \
    success &= b;                                                                                  \
  }
#define TEST_FE2(A, element, B)                                                                    \
  {                                                                                                \
    bool b = testFE(A, element, B);                                                                \
    std::cout << "testFE(" #A ", " #B ") " << (b ? "succeeded\n" : "failed\n");                    \
    success &= b;                                                                                  \
  }
#define TEST_FE3(A, element, B, C)                                                                 \
  {                                                                                                \
    bool b = testFE<C>(A, element, B);                                                             \
    std::cout << "testFE(" #A ", " #B ", " #C ") " << (b ? "succeeded\n" : "failed\n");            \
    success &= b;                                                                                  \
  }
#define TEST_FE4(A, element, B, C, D)                                                              \
  {                                                                                                \
    bool b = testFE<C>(A, element, B, D);                                                          \
    std::cout << "testFE(" #A ", " #B ", " #C ", " #D ") " << (b ? "succeeded\n" : "failed\n");    \
    success &= b;                                                                                  \
  }
#endif // DUNE_LOCALFUNCTIONS_TEST_TEST_BOUND_LOCALFE_HH
