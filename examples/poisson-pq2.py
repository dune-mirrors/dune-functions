import numpy as np
import scipy.sparse.linalg
from scipy.sparse import lil_matrix
from io import StringIO

import dune.geometry
import dune.grid as grid
import dune.functions as functions

# Compute element stiffness matrix and element load vector
#
# TODO: This assembler loop is very inefficient in terms of run time and should be improved using Python vectorization.
# See discussion at https://gitlab.dune-project.org/staging/dune-functions/-/merge_requests/295 for hints and code pointers.
def localAssembler(localView, volumeTerm):

    # Number of degrees of freedom on this element
    n = len(localView)

    # The grid element
    element = localView.element()

    # The set of shape functions
    localBasis = localView.tree().finiteElement().localBasis

    # Initialize element matrix and load vector
    localA = np.zeros((n,n))
    localB = np.zeros(n)

    # choose a high enough quadrature order
    quadOrder = 4

    # create a quadrature rule and integrate
    quadRule = dune.geometry.quadratureRule(element.type, quadOrder)
    for pt in quadRule:

        # Position of the quadrature point
        quadPos = pt.position

        # The determinant that appears in the integral transformation formula
        integrationElement = element.geometry.integrationElement(quadPos)

        # Evaluate all shape functions (The return value is an array!)
        values = localBasis.evaluateFunction(quadPos)

        # Evaluate the shape function Jacobians on the reference element (array of arrays)
        referenceJacobians = localBasis.evaluateJacobian(quadPos)

        # Transform the reference Jacobians to the actual element
        geometryJacobianInverse = element.geometry.jacobianInverse(quadPos)
        jacobians = [ np.dot(np.array(g)[0], geometryJacobianInverse) for g in referenceJacobians ]

        quadPosGlobal = element.geometry.toGlobal(quadPos)

        for i in range( n ):
            for j in range( n ):
                localA[i,j] += pt.weight * integrationElement * np.dot(jacobians[i], jacobians[j])

            localB[i] += pt.weight * integrationElement * values[i] * volumeTerm(quadPosGlobal)

    return localA, localB


# The assembler for the global stiffness matrix
def assembleLaplaceMatrix(basis, volumeTerm):

    # Total number of degrees of freedom
    n = len(basis)

    # Make empty sparse matrix
    A = lil_matrix((n,n))

    # Make empty rhs vector
    b = np.zeros(n)

    # View on the finite element basis on a single element
    localView = basis.localView()

    # Loop over all grid elements
    grid = basis.gridView
    for element in grid.elements:

        # Bind the localView to the current element
        localView.bind(element)

        # Number of degrees of freedom on the current element
        localN = len(localView)

        # Assemble the local stiffness matrix and load vector
        localA, localb = localAssembler(localView, volumeTerm)

        # Copy the local entries into the global matrix and vector
        for i in range(localN):

            gi = localView.index(i)[0]

            b[gi] += localb[i]

            for j in range(localN):
                gj = localView.index(j)[0]
                A[gi, gj] += localA[i, j]

    # Convert matrix to CSR format
    return A.tocsr(), b

# Mark all degrees of freedom on the grid boundary
#
# This method simply calls the corresponding C++ code.  A more Pythonic solution
# is planned to appear eventually...
def markBoundaryDOFs(basis, vector):
    code="""
    #include<utility>
    #include<functional>
    #include<dune/functions/functionspacebases/boundarydofs.hh>
    #include<dune/common/rangeutilities.hh>

    template<class Basis, class Vector>
    void run(const Basis& basis, Vector& vector)
    {
      {
        std::cout << basis.gridView().indexSet().size(0) << std::endl;
        auto localView = basis.localView();
        std::size_t k=0;
        for(const auto& e : elements(basis.gridView()))
        {
          localView.bind(e);
          std::cout << "Element " << k;
          std::cout << " index=" << basis.gridView().indexSet().index(e);
          std::cout << " dof indices [ ";
          for(auto i: Dune::range(localView.tree().size()))
            std::cout << localView.index(i) << " ";
          std::cout << " ]" << std::endl;
          ++k;
        }
      }

      auto vectorBackend = vector.mutable_unchecked();
      Dune::Functions::forEachBoundaryDOF(basis, [&] (auto&& index) {
          vectorBackend[index] = true;
      });
    }
    """
    dune.generator.algorithm.run("run",StringIO(code), basis, vector)


############################  main program  ###################################

# Number of grid elements in one direction
#gridSize = 32

# Create a grid of the unit square
#grid = grid.structuredGrid([0,0],[1,1],[gridSize,gridSize])
grid = dune.grid.ugGrid( (dune.grid.reader.gmsh, "square_triangles.msh"), dimgrid=2 )
grid.hierarchicalGrid.globalRefine(1)

# Create a second-order Lagrange FE basis
basis = functions.defaultGlobalBasis(grid, functions.Lagrange(order=2))

# Load term
rightHandSide = lambda x : 10

# Compute stiffness matrix and rhs vector
A,b = assembleLaplaceMatrix(basis, rightHandSide)

# Determine all coefficients that are on the boundary
isDirichlet = np.zeros(len(basis))

#def isNear(a,b):
#  return np.abs(a-b) <= 1e-10
#isDirichletIndicator = lambda x : isNear(x[0], 0) or isNear(x[0], 1) or isNear(x[1], 0) or isNear(x[1], 1)
#basis.interpolate(isDirichlet,isDirichletIndicator)

markBoundaryDOFs(basis, isDirichlet)

# The function that implements the Dirichlet values
dirichletValueFunction = lambda x : np.sin(2*np.pi*x[0])

# Get coefficients of a Lagrange-FE approximation of the Dirichlet values
dirichletCoeffs = np.zeros(len(basis))
basis.interpolate(dirichletCoeffs, dirichletValueFunction)

# Integrate Dirichlet conditions into the matrix and rhs vector
rows, cols = A.nonzero()

for i,j in zip(rows, cols):
    if isDirichlet[i]:
        if i==j:
            A[i,j] = 1.0
        else:
            A[i,j] = 0
        b[i] = dirichletCoeffs[i]


# Solve linear system!
x = scipy.sparse.linalg.spsolve(A, b)

# Write result as vtu file
u = basis.asFunction(x)

vtk = grid.vtkWriter(2)
u.addToVTKWriter("sol", vtk, dune.grid.DataType.PointData)
vtk.write("poisson-pq2-result")
