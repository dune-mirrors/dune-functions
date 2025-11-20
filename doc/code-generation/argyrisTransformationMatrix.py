# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

##### This file was used to generate the code for the transformation in dune/functions/functionspacebases/argyrisbasis.hh

from sympy import *
init_printing()

### Sympy generates c code, that treats a matrix as a flat vector. 
### Here we restore this to double bracket notation and treat pow(h,2) and o[j] for efficiancy
### dic is a dictionary with all submatrix names and their shapes
def restoreArrayNotation(code, dic):
  for name, shape in dic.items():
    for i in range(shape[0]):
      for j in range(shape[1]):
        code = code.replace(name+"["+str(i*shape[1]+j)+"]",
                            name+"["+str(i)+"fixed]["+str(j)+"]")
  code = code.replace("fixed", "")
  # Here 'h' should be replaced by a regex for generalization
  for j in range(3):
    code = code.replace("*pow(h["+str(j)+"], 2)",
                          "*h["+str(j)+"] * h["+str(j)+"]")
    code= code.replace("/pow(h["+str(j)+"], 2)",
                          "/h["+str(j)+"] /h["+str(j)+"]")
  for j in range(3):
    code= code.replace("o["+str(j)+"]",
                          "(o["+str(j)+"] ? -1 : 1)")
  return code

### Take a matrix, write it as bcrs matrix. Includes pattern setup
def writeBCRS(mat, replacementDict, path="argyrisTransformationMethods.txt"):
  with open(path, 'w') as f:
    for row in range(mat.shape[0]):
      print(mat.row(row).RL)
      f.write("mat_.setrowsize("+str(row)+","+str(len(mat.row(row).RL))+");\n")
    f.write("mat_.endrowsizes();\n")
    for row in range(mat.shape[0]):
      for element in mat.row(row).RL:
        f.write("mat_.addindex("+str(row)+","+str(element[1])+");\n")
    f.write("mat_.endindices();\n\n")

    for element in mat.RL:
      code = ccode(element[2], "mat_["+str(element[0]) +
                   "]["+str(element[1])+"]", standard="c99")
      codeRestored = restoreArrayNotation(code, replacementDict)
      f.write(codeRestored+"\n")

def writeMatVec(mat, replacementDict, path="argyrisTransformationMatrixFreeMethods.txt"):
  with open(path, 'w') as f:

    V = SparseMatrix(21,1,{(0,0):Matrix(MatrixSymbol('inValues', 21,1))})

    result = simplify( mat*V)

    for element in result.RL:
      code = ccode(element[2], "outValues["+str(element[0]) + "]", standard="c99")
      codeRestored = restoreArrayNotation(code, replacementDict)

      f.write(codeRestored+"\n")

if __name__="__main__":

  ### Main routine: create sympy matrix 
  ### We follow Kirbys notation in https://www.numdam.org/item/10.5802/smai-jcm.33.pdf

  ### E reduces the degrees of freedom to the original by removing tangential derivatives
  E = SparseMatrix(21, 24, {(0, 0): eye(19), (19, 20): 1, (20, 22): 1})
  #pprint(E)
  ## J transforms gradients from global xy to reference xy
  J_0 = Matrix(MatrixSymbol('J_0', 2, 2))
  J_1 = Matrix(MatrixSymbol('J_1', 2, 2))
  J_2 = Matrix(MatrixSymbol('J_2', 2, 2))
  ### Include Dir to transform gradients from general global directions
  Dir_0 = Matrix(MatrixSymbol('dir_0',2,2))
  Dir_1 = Matrix(MatrixSymbol('dir_1',2,2))
  Dir_2 = Matrix(MatrixSymbol('dir_2',2,2))

  ### Theta transforms Hessians 
  Theta_0 = Matrix(MatrixSymbol('theta_0', 3, 3))
  Theta_1 = Matrix(MatrixSymbol('theta_1', 3, 3))
  Theta_2 = Matrix(MatrixSymbol('theta_2', 3, 3))
  ### Again a direction modification
  ThetaDir_0 = Matrix(MatrixSymbol('thetaDir_0',3,3))
  ThetaDir_1 = Matrix(MatrixSymbol('thetaDir_1',3,3))
  ThetaDir_2 = Matrix(MatrixSymbol('thetaDir_2',3,3))
  ###  B transforms normal and tangential derivatives from global to local
  B1 = Matrix(MatrixSymbol('b_0', 2, 2))
  B2 = Matrix(MatrixSymbol('b_1', 2, 2))
  B3 = Matrix(MatrixSymbol('b_2', 2, 2))

  # Dir = SparseMatrix(21, 21, {(0, 0): 1, (1, 1): Dir_0.T, (3, 3): ThetaDir_0, (6, 6): 1, (7, 7): Dir_1.T, (
      # 9, 9): ThetaDir_1, (12, 12): 1, (13, 13): Dir_2.T, (15, 15): ThetaDir_2, (18, 18): eye(3)})

  ### Full transformation matrix for extended dofs
  VC = SparseMatrix(24, 24, {(0, 0): 1, (1, 1): J_0, (3, 3): Theta_0, (6, 6): 1, (7, 7): J_1, (
      9, 9): Theta_1, (12, 12): 1, (13, 13): J_2, (15, 15): Theta_2, (18, 18): B1, (20, 20): B2, (22, 22): B3})
  pprint(VC)
  
  ## edge lengths
  l = Matrix(MatrixSymbol('l', 3, 1))
  ## edge tangents
  t = Matrix(MatrixSymbol('t', 3, 2))
  ## once vector per edge
  tau = Matrix(MatrixSymbol('tau', 3, 3))
  ## boolean vector indicating wrongly oriented edges
  orientation = Matrix(MatrixSymbol('o', 3, 1))

  ### Matrix that extends the dofs by adding tangential derivatives at edge midpoints
  D = SparseMatrix(24, 21, {(0, 0): eye(19), (20, 19): 1, (22, 20): 1})
  shift = [[0, 6], [0, 12], [6, 12]]
  for i in range(3):
    D[19+2*i, shift[i][0]] = -15/(8*l[i])
    D[19+2*i, shift[i][0]+1] = -7/16 * t[i, :]
    D[19+2*i, shift[i][0]+3] = -l[i]/32 * tau[i, :]
    D[19+2*i, shift[i][1]] = 15/(8*l[i])
    D[19+2*i, shift[i][1]+1] = -7/16 * t[i, :]
    D[19+2*i, shift[i][1]+3] = l[i]/32 * tau[i, :]
  #pprint(D)
  # pprint(Dir)

  ### Which degrees of freedom to flip on a wrongly oriented edge
  OrientationMatrix = SparseMatrix(21,21,{(0,0):eye(18),(18,18):orientation[0],(19,19):orientation[1],(20,20):orientation[2]})

  ### matrix to scale dofs by inverted h-average of vertex patch
  h = Matrix(MatrixSymbol('h', 3,1))
  ## Not yet inverted
  ScalingMatrix = banded({0 : (1,h[0],h[0],h[0]*h[0], h[0]*h[0], h[0]*h[0],1,h[1],h[1],h[1]*h[1], h[1]*h[1], h[1]*h[1],1,h[2],h[2],h[2]*h[2], h[2]*h[2], h[2]*h[2],1,1,1)})

  ###### Now collect the matrix that relates global dofs and pushed forward local dofs, as well as reference basis and pulled back global basis 
  ## Pure V as in Kirbys Paper
  # V = E*VC*D
  # Also handling orientation
  # V = OrientationMatrix*D.T*VC*E.T
  # Also handling vertex scaling
  V = OrientationMatrix * ScalingMatrix.inv()*D.T*VC*E.T

  # Also handling Directions of Derivatives
  # V = OrientationMatrix*Dir*D.T*VC*E.T

  replacementDict = {"theta_0": (3, 3), "theta_1": (3, 3), "theta_2": (3, 3), "J_0": (2, 2),"J_1": (2, 2), "J_2": (2, 2),"thetaDir_0": (3, 3), "thetaDir_1": (3, 3), "thetaDir_2": (3, 3), "dir_0": (2, 2),"dir_1": (2, 2), "dir_2": (2, 2),  "b_1": (
          2, 2), "b_2": (2, 2), "b_0": (2, 2), 't': (3, 2), 'tau': (3, 3)}
  print("==============================================================")
  print("start writing file")
  writeBCRS(V, replacementDict)
  writeMatVec(V, replacementDict)
  # with open("argyrisMatrix.tex", 'w') as f:
  #  f.write(latex(V))
