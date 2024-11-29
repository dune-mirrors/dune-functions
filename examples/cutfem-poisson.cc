// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#include <config.h>

#include <vector>
#include <cmath>

#include <dune/common/bitsetvector.hh>
#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/multidomain.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/tpmc/tpmcrefinement.hh>

using namespace Dune;

// material properties...
std::vector<double> D{0.01,100};

// this is a quick and dirty implementation, using many
// shortcuts as we know we work on a YaspGrid in 2D
class CutCellInfo
{
public:
  using Coord = FieldVector<double,2>;
  using QP = QuadraturePoint<double,2>;
  using LevelSet = std::function<double(Coord)>;

  CutCellInfo(LevelSet l0, LevelSet l1) :
    levelset0(l0), levelset1(l1) {}

  template<typename Element>
  void bind(const Element & e)
  {
    std::array<double,4> values0;
    std::array<double,4> values1;
    bool outside = false;
    bool inbetween = false;
    bool innermost = false;
    for (int i=0; i<e.geometry().corners(); i++)
    {
      auto v0 = levelset0(e.geometry().corner(i));
      auto v1 = levelset1(e.geometry().corner(i));
      innermost = innermost || (v0 < 0);
      inbetween = inbetween || (v1 < 0 && v0 >= 0);
      outside   = outside   || (v1 >= 0);
      values0[i] = v0;
      values1[i] = v1;
    }
    // which domain have support?
    dom0 = innermost;
    dom1 = inbetween;
    // cut cell geometry
    interface0.emplace(values0);
    interface1.emplace(values1);
    // normal vectors
    interfaceNormal[0] = 0.5*((values0[1]-values0[0])+(values0[3]-values0[2]));
    interfaceNormal[1] = 0.5*((values0[2]+values0[3])-(values0[0]+values0[1]));
    interfaceNormal = interfaceNormal / interfaceNormal.two_norm();
    boundaryNormal[0] = 0.5*((values1[1]-values1[0])+(values1[3]-values1[2]));
    boundaryNormal[1] = 0.5*((values1[2]+values1[3])-(values1[0]+values1[1]));
    boundaryNormal = boundaryNormal / boundaryNormal.two_norm();
  }

  std::vector<QP> quadrature(int order, int domain) const
  {
    using Snippets = IteratorRange<typename TpmcRefinement<double,2>::const_volume_iterator>;
    Snippets snippets;
    if (!dom0 && !dom1)
      return {};
    if (domain == 0)
      snippets = interface0->volume(tpmc::ReconstructionType::InteriorDomain);
    if (domain == 1 && dom0)
      snippets = interface0->volume(tpmc::ReconstructionType::ExteriorDomain);
    if (domain == 1 && ! dom0)
      snippets = interface1->volume(tpmc::ReconstructionType::InteriorDomain);
    // create quadrature
    std::vector<QP> newquad;
    for (const auto& geo : snippets)
    {
      const auto& quad = QuadratureRules<double,2>::rule(geo.type(), order);
      for (const auto& qp : quad)
      {
        // transform quadrature point
        newquad.push_back(
          QuadraturePoint<double,2>(
            geo.global(qp.position()),
            geo.integrationElement(qp.position())*qp.weight()));
      }
    }
    return newquad;
  }
  FieldVector<double,2> interfaceNormal;
  FieldVector<double,2> boundaryNormal;
private:
  std::optional<TpmcRefinement<double,2>> interface0;
  std::optional<TpmcRefinement<double,2>> interface1;
  LevelSet levelset0;
  LevelSet levelset1;
  bool dom0;
  bool dom1;
};

// Compute the stiffness matrix for a single element
template <class Element, class FENode, class MatrixType>
void assemblePoisson(const Element& element, const FENode& feNode, MatrixType& elementMatrix, const CutCellInfo& ccinfo, int subDomain)
{

  const int dim = Element::dimension;
  auto geometry = element.geometry();

  const auto& localFiniteElement = feNode.finiteElement();
  if (! localFiniteElement.active()) return;

  // Get a quadrature rule
  int order = 2*(dim*localFiniteElement.localBasis().order()-1);
  // const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(element.type(), order);
  auto quad = ccinfo.quadrature(order, subDomain);

  // Loop over all quadrature points
  for (auto qp : quad) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = qp.position();

    // The inverse Jacobian of the map from the reference element to the element
    const auto& jacobianInverse = geometry.jacobianInverse(quadPos);

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = geometry.integrationElement(quadPos);

    // The gradients of the shape functions on the reference element
    std::vector<FieldMatrix<double,1,dim> > referenceJacobians;
    localFiniteElement.localBasis().evaluateJacobian(quadPos, referenceJacobians);

    // Compute the shape function gradients on the real element
    std::vector<FieldMatrix<double,1,dim> > jacobians(referenceJacobians.size());
    for (size_t i=0; i<jacobians.size(); i++)
      jacobians[i] = referenceJacobians[i] * jacobianInverse;

    // Compute the actual matrix entries
    for (size_t i=0; i<feNode.size(); i++)
      for (size_t j=0; j<feNode.size(); j++ )
      {
        size_t row = feNode.localIndex(i);
        size_t col = feNode.localIndex(j);
        elementMatrix[row][col] +=
          (D[subDomain] * jacobians[i] * transpose(jacobians[j]))
          * qp.weight() * integrationElement;
      }

  }
}

// Compute the stiffness matrix for a single element
template <class Element, class FENode, class MatrixType>
void assembleCouplings(const Element& element, const FENode& node, MatrixType& elementMatrix, const CutCellInfo& ccinfo)
{
  using namespace Indices;

  const auto& node0 = node.child(_0);
  const auto& node1 = node.child(_1);
  const auto& lfe0 = node0.finiteElement();
  const auto& lfe1 = node1.finiteElement();

  // only do something if both function spaces are active
  if (node0.size() == 0 || node1.size() == 0) return;

  // assemble nitsche term

  const int dim = Element::dimension;
  auto geometry = element.geometry();

  // pseudo intersection (we know there is currently only one!)
  double weight = std::sqrt(2.0);
  Dune::FieldVector<double,2> normal {weight,weight};

#if 0
  for (auto&& intersection : interface)
  {
    // Get a quadrature rule
    int order = 3;
    const auto& quad = QuadratureRules<double, dim-1>::rule(intersection.type(), order);

    // Loop over all quadrature points
    for (size_t pt=0; pt < quad.size(); pt++) {

      auto positionInElement = intersection.geometryInInside().global(quadPoint.position());

      // Position of the current quadrature point in the reference element
      const FieldVector<double,dim>& quadPos = quad[pt].position();

      // The inverse Jacobian of the map from the reference element to the element
      const auto& jacobianInverse = geometry.jacobianInverse(quadPos);

      // The values of the shape functions at the quadrature point
      // --- we know both subdomains use the same shape functions
      std::vector<FieldVector<double,1> > values;
      lfe0.localBasis().evaluateFunction(positionInElement, values);

      // The gradients of the shape functions on the reference element
      // --- we know both subdomains use the same shape functions
      std::vector<FieldMatrix<double,1,dim> > referenceJacobians;
      lfe0.localBasis().evaluateJacobian(quadPos, referenceJacobians);

      // Compute the shape function gradients on the real element
      std::vector<FieldMatrix<double,1,dim> > jacobians(referenceJacobians.size());
      for (size_t i=0; i<jacobians.size(); i++)
        jacobians[i] = referenceJacobians[i] * jacobianInverse;

      // The multiplicative factor in the integral transformation formula
      const double integrationElement = geometry.integrationElement(quadPos);

      // auto factor = intersection.integrationOuterNormal(quadPoint.position());
      auto n = intersection.integrationOuterNormal(quadPoint.position());

      // Compute the actual matrix entries
      for (size_t i=0; i<values.size(); i++)
        for (size_t j=0; j<values.size(); j++ )
        {
          // row: test / col: trial
          // inside/inside
          {
            size_t row = lfe0.localIndex(i);
            size_t col = lfe0.localIndex(j);
            elementMatrix[row][col] += (jacobians[i] * transpose(jacobians[j])) * quad[pt].weight() * integrationElement;
          }
          // inside/outside
          {
            size_t row = lfe0.localIndex(i);
            size_t col = lfe1.localIndex(j);
            elementMatrix[row][col] += (jacobians[i] * transpose(jacobians[j])) * quad[pt].weight() * integrationElement;
          }
          // outside/inside
          {
            size_t row = lfe1.localIndex(i);
            size_t col = lfe0.localIndex(j);
            elementMatrix[row][col] += (jacobians[i] * transpose(jacobians[j])) * quad[pt].weight() * integrationElement;
          }
          // outside/outside
          {
            size_t row = lfe1.localIndex(i);
            size_t col = lfe1.localIndex(j);
            elementMatrix[row][col] += (jacobians[i] * transpose(jacobians[j])) * quad[pt].weight() * integrationElement;
          }
        }
    }
  }
#endif
}

// Compute the stiffness matrix for a single element
template <class LocalView, class MatrixType>
void getLocalMatrix(const LocalView& localView, const CutCellInfo& ccinfo, MatrixType& elementMatrix)
{
  // Get the grid element from the local FE basis view
  typedef typename LocalView::Element Element;
  const Element& element = localView.element();

  // Set all matrix entries to zero
  elementMatrix.setSize(localView.size(), localView.size());
  elementMatrix = 0;      // fills the entire matrix with zeros

  // Get FE tree
  const auto& feTree = localView.tree();

  using namespace Indices;

  if (feTree.child(_0).active())
    assemblePoisson(element, feTree.child(_0), elementMatrix, ccinfo, 0); // inner
  if (feTree.child(_1).active())
    assemblePoisson(element, feTree.child(_1), elementMatrix, ccinfo, 1); // outer
  if (feTree.child(_0).active() && feTree.child(_1).active())
    assembleCouplings(element, feTree, elementMatrix, ccinfo);
}


// Compute the source term for a single element
template <class LocalView, class LocalVolumeTerm>
void getSourceTerm( const LocalView& localView,
                    BlockVector<FieldVector<double,1> >& localRhs,
                    LocalVolumeTerm&& localVolumeTerm)
{
  using namespace Indices;

  // Get the grid element from the local FE basis view
  typedef typename LocalView::Element Element;
  const Element& element = localView.element();

  const int dim = Element::dimension;

  // Get set of shape functions for this element
  const auto& node0 = localView.tree().child(_0);
  const auto& node1 = localView.tree().child(_1);
  const auto& lfe0 = node0.finiteElement();
  const auto& lfe1 = node1.finiteElement();

  // Set all entries to zero
  localRhs.resize(localView.tree().size());
  localRhs = 0;

  // A quadrature rule
  int order = dim*1;
  const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(element.type(), order);

  // storage for shape function values
  std::vector<FieldVector<double,1> > shapeFunctionValues;

  // Loop over all quadrature points
  for ( size_t pt=0; pt < quad.size(); pt++ ) {

    // Position of the current quadrature point in the reference element
    const FieldVector<double,dim>& quadPos = quad[pt].position();

    // The multiplicative factor in the integral transformation formula
    const double integrationElement = element.geometry().integrationElement(quadPos);

    double functionValue = localVolumeTerm(quadPos);

    {
      // Evaluate all shape function values at this point
      lfe0.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

      // Actually compute the vector entries
      for (size_t i=0; i<shapeFunctionValues.size(); i++)
      {
        size_t row = node0.localIndex(i);
        localRhs[i] += D[0] * shapeFunctionValues[i] * functionValue * quad[pt].weight() * integrationElement;
      }
    }

    {
      // Evaluate all shape function values at this point
      lfe1.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

      // Actually compute the vector entries
      // Somehow this size differs from node1.size()
      for (size_t i=0; i<shapeFunctionValues.size(); i++)
      {
        size_t row = node1.localIndex(i);
        localRhs[i] += D[1] * shapeFunctionValues[i] * functionValue * quad[pt].weight() * integrationElement;
      }
    }
  }

}

// Get the occupation pattern of the stiffness matrix
template <class FEBasis>
void getOccupationPattern(const FEBasis& feBasis, MatrixIndexSet& nb)
{
  // Total number of grid vertices
  auto n = feBasis.size();

  nb.resize(n, n);

  // A view on the FE basis on a single element
  auto localView = feBasis.localView();

  const auto & gridView = feBasis.gridView();

  // Loop over all leaf elements
  for(const auto& e : elements(gridView))
  {
    // Bind the local FE basis view to the current element
    localView.bind(e);

    // Cell Pattern, including Nitsche term

    // There is a matrix entry a_ij if the i-th and j-th vertex are connected in the grid
    for (size_t i=0; i<localView.tree().size(); i++) {

      for (size_t j=0; j<localView.tree().size(); j++) {

        auto iIdx = localView.index(i);
        auto jIdx = localView.index(j);

        // Add a nonzero entry to the matrix
        nb.add(iIdx, jIdx);

      }

    }

    // Skeleton pattern for ghost penalty
    if (gridView.domains(e).size() <= 1)
      continue;

    // Now let's get the off-diagonal element stiffness matrix
    for (auto&& is : intersections(localView.globalBasis().gridView(), e))
    {
      if (!is.neighbor())
        continue;

      if (gridView.domains(is.outside()).size() <= 1)
        continue;

      // Get a local view and local index set for the element on the other side of the intersection
      auto outsideLocalView = feBasis.localView();
      outsideLocalView.bind(is.outside());

      // Add element stiffness matrix onto the global stiffness matrix
      for (size_t i=0; i<localView.tree().size(); i++)
      {
        // The global index of the i-th local degree of freedom of the element 'e'
        auto row = localView.index(i)[0];

        for (size_t j=0; j<outsideLocalView.tree().size(); j++ )
        {
          // The global index of the j-th local degree of freedom
          // of the element on the other side of the intersection
          auto col = outsideLocalView.index(j)[0];
          nb.add(row,col);
        }

      }

    }

  }

}


/** \brief Assemble the Laplace stiffness matrix on the given grid view */
template <class FEBasis, class VolumeTerm>
void assemble(const FEBasis& feBasis,
              BCRSMatrix<double>& matrix,
              BlockVector<double>& rhs,
              VolumeTerm&& volumeTerm,
              CutCellInfo& ccinfo)
{
  // Get the grid view from the finite element basis
  typedef typename FEBasis::GridView GridView;
  GridView gridView = feBasis.gridView();

  auto localVolumeTerm = localFunction(Functions::makeGridViewFunction(volumeTerm, gridView));

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  MatrixIndexSet occupationPattern;
  getOccupationPattern(feBasis, occupationPattern);

  // ... and give it the occupation pattern we want.
  occupationPattern.exportIdx(matrix);

  // set rhs to correct length -- the total number of basis vectors in the basis
  rhs.resize(feBasis.size());

  // Set all entries to zero
  matrix = 0;
  rhs = 0;

  // A view on the FE basis on a single element
  auto localView = feBasis.localView();

  // A loop over all elements of the grid
  for(const auto& e : elements(gridView))
  {

    // Bind the local FE basis view to the current element
    localView.bind(e);

    // bind cutcell infos
    ccinfo.bind(e);

    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    Matrix<FieldMatrix<double,1,1> > elementMatrix;
    getLocalMatrix(localView, ccinfo, elementMatrix);

    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i=0; i<elementMatrix.N(); i++) {

      // The global index of the i-th local degree of freedom of the element 'e'
      auto row = localView.index(i);

      for (size_t j=0; j<elementMatrix.M(); j++ ) {

        // The global index of the j-th local degree of freedom of the element 'e'
        auto col = localView.index(j);
        matrix[row][col] += elementMatrix[i][j];

      }

    }

    // Skeleton pattern for ghost penalty
    if (gridView.domains(e).size() <= 1)
      continue;

    // Now let's get the off-diagonal element stiffness matrix
    for (auto&& is : intersections(gridView, e))
    {
      if (!is.neighbor())
        continue;

      if (gridView.domains(is.outside()).size() <= 1)
        continue;

      // Get a local view and local index set for the element on the other side of the intersection
      auto outsideLocalView = feBasis.localView();
      outsideLocalView.bind(is.outside());

      // assembleStabilization(...);
      // getOffDiagonalLocalMatrix(is, localView, outsideLocalView, elementMatrix, localVelocityField);

      // Add element stiffness matrix onto the global stiffness matrix
      for (size_t i=0; i<elementMatrix.N(); i++) {

        // The global index of the i-th local degree of freedom of the element 'e'
        auto row = localView.index(i)[0];

        for (size_t j=0; j<elementMatrix.M(); j++ ) {

          // The global index of the j-th local degree of freedom
          // of the element on the other side of the intersection
          auto col = outsideLocalView.index(j)[0];
          matrix[row][col] += elementMatrix[i][j];

        }
      }
    }

    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> > localRhs;
    localVolumeTerm.bind(e);
    getSourceTerm(localView, localRhs, localVolumeTerm);

    for (size_t i=0; i<localRhs.size(); i++) {

      // The global index of the i-th vertex of the element 'e'
      auto row = localView.index(i);
      rhs[row] += localRhs[i];

    }

  }

}


// This method marks all Lagrange nodes on the boundary of the grid.
// In our problem we want to use them as the Dirichlet nodes.
// The result can be found in the 'dirichletNodes' variable.  There, a bit
// is set precisely when the corresponding Lagrange node is on the grid boundary.
//
// Since interpolating into a vector<bool> is currently not supported,
// we use a vector<char> which, in contrast to vector<bool>
// is a real container.
template <class FEBasis>
void boundaryTreatment (const FEBasis& feBasis, std::vector<char>& dirichletNodes )
{
  dirichletNodes.clear();
  dirichletNodes.resize(feBasis.size(), false);

  Functions::forEachBoundaryDOF(feBasis, [&] (auto&& index) {
    dirichletNodes[index] = true;
  });
}

///
template<typename GridView, typename F>
auto createDomainInfo(const GridView& gridView, const F& levelSet0, const F& levelSet1)
{
  // auto lv0 = localFunction(levelSet0);
  // auto lv1 = localFunction(levelSet1);

  // FLAGS:
  // innermost = 1
  // inbetween = 2
  // outside   = 4
  //
  // partition = sum(FLAGS)
  // innermost = 1 -> 0
  // Interface = 3    1
  // inbetween = 2    2
  // Interface = 6    3
  // outside   = 4    4
  std::vector<int> partitions(gridView.size(0), 0);
  auto const& indexSet = gridView.indexSet();
  for (auto const& e : elements(gridView)) {
    // lv0.bind(e);
    // lv1.bind(e);
    bool outside = false;
    bool inbetween = false;
    bool innermost = false;
    for (int i = 0; i < e.geometry().corners(); i++)
    {
      auto v0 = levelSet0(e.geometry().corner(i));
      auto v1 = levelSet1(e.geometry().corner(i));
      innermost = innermost || (v0 < 0);
      inbetween = inbetween || (v1 < 0 && v0 >= 0);
      outside   = outside   || (v1 >= 0);
    }
    //                  0, 1, 2, 3, 4, 5, 6
    int partition[] = {-1, 0, 2, 1, 4,-1, 3};
    int key = (outside << 2) + (inbetween << 1) + innermost;
    partitions[indexSet.index(e)] = partition[key];
    if (partition[key] < 0)
      std::cout << "strange cell : "
                << outside << ", "
                << inbetween << ", "
                << innermost << " -> "
                << key << " -> " << partition[key] << std::endl;
    assert(partition[key] >= 0);
  }

  // subdomains:
  // domain0 (inside) = {0,1}, domain1 (outside) = {1,2,3}
  // domain2 (interface) = {1}, domain3 (boundary) = {3}
  using namespace Dune::Functions::MultiDomain;
  auto domainInfo = createPartitionedDomainInfo(
    std::move(partitions),
    {{0,1}, {1,2,3}, {1}, {3}});

  return domainInfo;
}

int main (int argc, char *argv[]) try
{
  using namespace Dune::Functions::BasisFactory;
  using namespace Dune::Functions::MultiDomain;

  // Set up MPI, if available
  MPIHelper::instance(argc, argv);

  ///////////////////////////////////
  //   Generate the grid
  ///////////////////////////////////

  constexpr int dim = 2;
  double center = 0.5;
  Dune::YaspGrid<dim> grid{ {1.0,1.0}, {3,3} };
  grid.globalRefine(4);
  auto gridView = grid.leafGridView();
  using GridView = decltype(gridView);

  ///////////////////////////////////
  //   Define Geometry and partition
  //   based on level-sets
  ///////////////////////////////////

  double R0 = 0.3;
  double R1 = 0.45;
  double freq = 10;
  double ampl = 0.3;
  auto l0 = [&](auto x) {
    x -= center;
    return 10*(x*x - R0*R0) + ampl * std::sin(freq * std::atan2(x[0],x[1]));
  };
  auto l1 = [&](auto x) {
    x -= center;
    return x*x - R1*R1;
  };

  auto levelSetBasis = makeBasis(gridView,lagrange<1>());

  using LSVectorType = BlockVector<double>;
  LSVectorType lx0;
  LSVectorType lx1;
  interpolate(levelSetBasis, lx0, l0);
  interpolate(levelSetBasis, lx1, l1);

  auto levelset0 = Functions::makeDiscreteGlobalBasisFunction<double>(levelSetBasis, lx0);
  auto levelset1 = Functions::makeDiscreteGlobalBasisFunction<double>(levelSetBasis, lx1);

  CutCellInfo ccinfo(levelset0, levelset1);

  /////////////////////////////////////////////////////////
  //   Choose a finite element space
  /////////////////////////////////////////////////////////
  auto domainInfo = createDomainInfo(gridView, levelset0, levelset1);

  auto multiDomainPreBasis =
    multiDomainComposite(domainInfo,
      restrict(lagrange<1>(), subdomain(0)),
      restrict(lagrange<1>(), subdomain(1)),
      flatLexicographic());

  auto basis = makeBasis(gridView,multiDomainPreBasis);

  /////////////////////////////////////////////////////////
  //   Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////

  using VectorType = BlockVector<double>;
  using MatrixType = BCRSMatrix<double>;

  VectorType rhs;
  MatrixType stiffnessMatrix;

  auto pi = StandardMathematicalConstants<double>::pi();

  /////////////////////////////////////////////////////////
  //  Assemble the system
  /////////////////////////////////////////////////////////

  std::cout << "Number of DOFs is " << basis.dimension() << std::endl;

  auto rightHandSide = [pi] (const auto& x) { return 4*pi*std::sin(2*pi*x[1]); };

  Dune::Timer timer;
  assemble(basis, stiffnessMatrix, rhs, rightHandSide, ccinfo);
  std::cout << "Assembling the problem took " << timer.elapsed() << "s" << std::endl;

  /////////////////////////////////////////////////
  //   Choose an initial iterate
  /////////////////////////////////////////////////
  VectorType x(basis.size());
  x = 0;

  // Determine Dirichlet dofs
  std::vector<char> dirichletNodes;
  boundaryTreatment(basis, dirichletNodes);

  // Don't trust on non-standard M_PI.
  auto dirichletValueFunction = [pi](const auto& x){ return std::sin(2*pi*x[1]); };

  // Interpolate dirichlet values at the boundary nodes
  interpolate(basis, x, dirichletValueFunction, dirichletNodes);

  //////////////////////////////////////////////////////
  //   Incorporate Dirichlet values in a symmetric way
  //////////////////////////////////////////////////////

  // Compute residual for non-homogeneous Dirichlet values
  // stored in x. Since x is zero for non-Dirichlet DOFs,
  // we can simply multiply by the matrix.
  stiffnessMatrix.mmv(x, rhs);

  // Change Dirichlet matrix rows and columns to the identity.
  for (size_t i=0; i<stiffnessMatrix.N(); i++) {
    if (dirichletNodes[i]) {
      rhs[i] = x[i];
      // loop over nonzero matrix entries in current row
      for (auto&& [entry, idx] : sparseRange(stiffnessMatrix[i]))
        entry = (i==idx) ? 1.0 : 0.0;
    }
    else {
      // loop over nonzero matrix entries in current row
      for (auto&& [entry, idx] : sparseRange(stiffnessMatrix[i]))
        if (dirichletNodes[idx])
          entry = 0.0;
    }
  }

  ////////////////////////////
  //   Compute solution
  ////////////////////////////

  // Technicality:  turn the matrix into a linear operator
  MatrixAdapter<MatrixType,VectorType,VectorType> op(stiffnessMatrix);

  // Sequential incomplete LU decomposition as the preconditioner
  SeqILDL<MatrixType,VectorType,VectorType> ildl(stiffnessMatrix,1.0);

  // Preconditioned conjugate-gradient solver
  CGSolver<VectorType> cg(op,
                          ildl, // preconditioner
                          1e-20, // desired residual reduction factor
                          1000,   // maximum number of iterations
                          2);   // verbosity of the solver

  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;

  // Solve!
  cg.apply(x, rhs, statistics);

  ////////////////////////////////////////////////////////////////////////////
  //  Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////
  using namespace Indices;

  auto x0Function
          = Functions::makeDiscreteGlobalBasisFunction<double>(
            Functions::subspaceBasis(basis, _0), x);
  auto x1Function
          = Functions::makeDiscreteGlobalBasisFunction<double>(
            Functions::subspaceBasis(basis, _1), x);

  //////////////////////////////////////////////////////////////////////////////////////////////
  //  Write result to VTK file
  //  We need to subsample, because VTK cannot natively display real second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////
  SubsamplingVTKWriter<GridView> vtkWriter(gridView, Dune::refinementLevels(0));
  vtkWriter.addVertexData(
          x0Function,
          VTK::FieldInfo("x0", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.addVertexData(
          x1Function,
          VTK::FieldInfo("x1", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.addVertexData(
          levelset0,
          VTK::FieldInfo("levelset0", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.addVertexData(
          levelset1,
          VTK::FieldInfo("levelset1", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("cutfem-poisson");

 }
// Error handling
 catch (Exception& e) {
    std::cout << e.what() << std::endl;
 }
