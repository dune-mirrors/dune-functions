# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

"""Generate the Arnold-Winther reference basis from Symfem."""

import symfem
from sympy import Matrix, cxxcode, diff, init_printing, shape
from sympy.abc import x, y
from sympy.polys.polyfuncs import horner


def createGenericReferenceElement(refName, feName, order, **kwargs):
  fe = symfem.create_element(refName, feName, order, **kwargs)
  assert fe.reference.name == "triangle"
  assert fe.reference.edges == ((0, 1), (0, 2), (1, 2))
  return fe

def verifyDuneDOFOrdering(fe):
  """Verify the local-key order assumed by ArnoldWintherLocalCoefficients."""
  def sympy_value(value):
    return value.as_sympy() if hasattr(value, "as_sympy") else value

  assert fe.order == 2
  assert fe.variant == "dune"
  assert len(fe.dofs) == 24

  expected_entities = (
      [(0, vertex) for vertex in range(3) for _ in range(3)]
      + [(1, edge) for edge in range(3) for _ in range(4)]
      + [(2, 0)] * 3
  )
  assert [dof.entity for dof in fe.dofs] == expected_entities

  expected_vertex_components = (
      ((1, 0), (1, 0)),
      ((1, 0), (0, 1)),
      ((0, 1), (0, 1)),
  )
  for vertex in range(3):
    for component, (left, right) in enumerate(expected_vertex_components):
      dof = fe.dofs[3 * vertex + component]
      assert sympy_value(dof.lvec) == left
      assert sympy_value(dof.rvec) == right

  for edge in range(3):
    edge_reference = fe.reference.sub_entity(1, edge)
    normal = sympy_value(edge_reference.normal())
    tangent = sympy_value(edge_reference.tangent())
    for moment in range(2):
      normal_dof = fe.dofs[9 + 4 * edge + 2 * moment]
      tangent_dof = fe.dofs[10 + 4 * edge + 2 * moment]
      assert normal_dof.dof.dof_point() == (moment,)
      assert tangent_dof.dof.dof_point() == (moment,)
      assert sympy_value(normal_dof.inner_with_left) == normal
      assert sympy_value(normal_dof.inner_with_right) == normal
      assert sympy_value(tangent_dof.inner_with_left) == tangent
      assert sympy_value(tangent_dof.inner_with_right) == normal

  expected_cell_components = (
      Matrix(((1, 0), (0, 0))),
      Matrix(((0, 1), (0, 0))),
      Matrix(((0, 0), (0, 1))),
  )
  assert tuple(dof.f.as_sympy() for dof in fe.dofs[21:]) == expected_cell_components

# Apply a Horner scheme to a tensor-valued function.
def hornerScheme(f, derivative=None, **kwargs):
  s = shape(f)
  assert(len(s) == 2)

  if derivative == "Divdiv":
    result = horner(diff(f[0,0].diff(x) + f[0,1].diff(y), x) +  diff(f[1,0].diff(x) + f[1,1].diff(y), y))
  else:
    result =  [] # result is at least tensor order 1
    for i in range(s[0]):
      if derivative == "Div":
        result.append(horner(f[i,0].diff(x) + f[i,1].diff(y)))
      else:
        result.append([]) # result is at least tensor order 2
        for j in range(s[1]):
          if derivative is None: #values
            result[i].append(horner(f[i,j]))
          elif isinstance(derivative, list):  #multiple derivatives, typically jacobian
            result[i].append([])  ## tensor order 3
            for k,direction in enumerate(derivative):
              result[i][j].append(horner(diff(f[i,j], direction), **kwargs))
          else:   # single derivative, i.e. partial
            result[i].append(horner(diff(f[i,j], derivative)))

  return result

## get the code for a list (or a single function)
def getCodeForList(f, **kwargs):
  code = ""
  if isinstance(f, list):
    code+= "{" + cxxcode(f[0], **kwargs)
    for ff in f[1:]:
      code += ", " + getCodeForList(ff, **kwargs)
    code += "}"
  else:
    code += cxxcode(f, **kwargs)
  return code

def getCodeForScalarorVector(f, **kwargs):
  code = ""
  code += "\n*(iter++) = " + getCodeForList(f , **kwargs) + ";\n"
  return code

def getCodeForMatrix(tensor, **kwargs):
  symmetric = kwargs.pop("symmetric", False)
  code = ""
  if symmetric:
    code += "\n*(iter++) = sym<Range>(" + getCodeForList(tensor[0][0] , **kwargs) + ", " + getCodeForList( tensor[0][1], **kwargs) + ", " + getCodeForList( tensor[1][1], **kwargs) + ");"

  else :
    code += "\n*(iter++) = {{" + getCodeForList(tensor[0][0] , **kwargs)

    code += ",\t" + getCodeForList( tensor[0][1], **kwargs)+ "},\n\t{"  + getCodeForList( tensor[1][0], **kwargs)+  ",\t"
    code += getCodeForList( tensor[1][1], **kwargs)+"}};\n"

  return code

### function body for evaluateFunction
def getCodeForEvaluation(basis, **kwargs):
  symmetric = kwargs.pop("symmetric")
  derivative = kwargs.pop("derivative", None)
  powsubs={'Pow': [(lambda b,e: e == 2, lambda b, e: ("{0}*{0}".format(b))),
  (lambda b,e: e == 3, lambda b, e: ("{0}*{0}*{0}".format(b))),
  (lambda b,e: e == 4, lambda b, e: ("{0}*{0}*{0}*{0}".format(b))),
  (lambda b, e: not b.is_integer, 'pow')]}

  code = ""
  kwargs.update( {"user_functions": powsubs})
  name = "val" if kwargs.get("assign_to") is None else kwargs.get("assign_to").name
  # kwargs.pop("assign_to")
  for i,f in enumerate(basis):
    code += "\n//{}th basis function".format(i)

    mat = hornerScheme(f, derivative)
    if derivative is None:
      code += getCodeForMatrix(mat,symmetric = symmetric, **kwargs)+"\n"
    else:
      code += getCodeForScalarorVector(mat, **kwargs) + "\n"
  return code

## Generate an include file for evaluation methods
def printEvaluationCode(name, reference, feType,minOrder  = 0, maxOrder  = 3, symmetric = False, **kwargs):
  variant = kwargs.pop("variant", "dune")
  assert(isinstance(name, str))
  guard_name = name.upper().replace("REFERENCE", "BASIS")
  # Assemble the marker at runtime so REUSE does not mistake these generated
  # annotations for additional annotations of this Python source file.
  spdxMarker = "SPDX"
  code = f"// {spdxMarker}-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md\n"
  code += f"// {spdxMarker}-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later\n\n"
  code += "#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_" + guard_name + "_INC_HH\n"
  code += "#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_" + guard_name + "_INC_HH\n"
  code += "namespace Dune::Functions{\n  namespace Impl{ \n    "
  code += "template<class D, class R,int dim, unsigned int k>\n"
  code +='     void ' + name + 'LocalBasis<D,R, dim,k>::evaluateFunction(const typename Traits::DomainType &in,std::vector<typename Traits::RangeType> &out) const\n{\nout.resize(size());\n auto iter = out.begin();'
  code += "\n\n// generated with sympy from symfem library\n"
  code += "auto const&x = in[0], y = in[1];"
  for i in range(minOrder, maxOrder +1):
    fe = createGenericReferenceElement(reference, feType, i, variant=variant, **kwargs)
    verifyDuneDOFOrdering(fe)
    basis = fe.get_basis_functions()
    code += "\n if constexpr (k =="+str(i)+"){\n"
    code += getCodeForEvaluation(basis, symmetric = symmetric)

    code += "\n}"

  code +="\n}"
  code += 'template<class D, class R, int dim, unsigned int k>\n      void ' + name + 'LocalBasis<D,R,dim,k>::evaluateDivergence(const typename Traits::DomainType &in,std::vector<typename Traits::DivergenceType> &out) const\n{\nout.resize(size());\nauto iter = out.begin();'
  code += "\n\n// generated with sympy from symfem library\n"
  code += "auto const&x = in[0], y = in[1];"
  for i in range(minOrder, maxOrder +1):
    fe = createGenericReferenceElement(reference, feType, i, variant=variant, **kwargs)
    verifyDuneDOFOrdering(fe)
    basis = fe.get_basis_functions()
    code += "if constexpr (k =="+str(i)+"){"

    code += getCodeForEvaluation(basis, derivative = "Div", symmetric = symmetric)

    code += "\n}"
  code += "\n}"

  print(code)

  print( "  }//namespace Impl\n}//namespace Dune\n#endif" ) # namespace braces

if __name__== "__main__":
  init_printing()

  printEvaluationCode("ArnoldWintherReference", reference="triangle", feType="AW",
                      minOrder=2, maxOrder=2, symmetric=True, variant="dune")
