# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
from sympy import *
from sympy.polys.polyfuncs import horner
from sympy.abc import x, y
import symfem
import numpy as np

from symfem.functions import parse_function_input as parse, _to_sympy_format

def adaptReferenceElementToDune(fe):
  newRef = fe.reference
  assert(newRef.name == "triangle")
  newRef.edges = ((0,1),(0,2),(1,2))
  if hasattr(fe, "variant"):
    newFe = type(fe)(newRef, fe.order, fe.variant)
  else:
    newFe = type(fe)(newRef, fe.order)
  return newFe

def adaptReferenceToPhysicalElement(fe, vertices):
  physRef = symfem.create_reference("triangle", vertices= vertices)
  assert(physRef.name == "triangle")
  # physRef.edges = ((0,1),(0,2),(1,2))
  newFe = type(fe)(physRef, fe.order)
  return newFe

def createGenericReferenceElement(refName, feName, order, **kwargs):
  return(adaptReferenceElementToDune(symfem.create_element(refName, feName,order, **kwargs)))

def createPhysicalElement(refName, feName, order, vertices):
  return(adaptReferenceToPhysicalElement(symfem.create_element(refName, feName, order), vertices))

## apply a horner scheme on a function f
def hornerScheme(f, derivative = [x,y], **kwargs):
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
    code += "\n*(iter++) = sym<Range>(" + getCodeForList(tensor[0][0] , **kwargs) \
    + ", " + getCodeForList( tensor[0][1], **kwargs)\
    + ", " + getCodeForList( tensor[1][1], **kwargs)\
    + ");"

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

  assert(isinstance(name, str))
  code = "// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-\n \
// vi: set et ts=4 sw=2 sts=2: \n\n \
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md\n\
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later\n\
// NOTE: This is an auto generated file from CodeGeneration_HHJ.py\n"
  code += "#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_" + name.upper() + "_INC_HH\n#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_" + name.upper() + "_INC_HH\n namespace Dune::Functions{\n  namespace Impl{ \n    "
  code += "template<class D, class R,int dim, unsigned int k>\n"
  code +='     void ' + name + 'LocalBasis<D,R, dim,k>::evaluateFunction(const typename Traits::DomainType &in,std::vector<typename Traits::RangeType> &out) const\n{\nout.resize(size());\n auto iter = out.begin();'
  code += "\n\n// generated with sympy from symfem library\n"
  code += "auto const&x = in[0], y = in[1];"
  for i in range(minOrder, maxOrder +1):
    fe = createGenericReferenceElement(reference, feType, i, **kwargs)
    basis = fe.get_basis_functions()
    code += "\n if constexpr (k =="+str(i)+"){\n"
    code += getCodeForEvaluation(basis, symmetric = symmetric)

    code += "\n}"

  code +="\n}"
  code += 'template<class D, class R, int dim, unsigned int k>\n      void ' + name + 'LocalBasis<D,R,dim,k>::evaluateDivDiv(const typename Traits::DomainType &in,std::vector<typename Traits::DivDivType> &out) const\n{\nout.resize(size());\nauto iter = out.begin();'
  code += "\n\n// generated with sympy from symfem library\n"
  code += "auto const&x = in[0], y = in[1];"
  for i in range(minOrder, maxOrder +1):
    fe = createGenericReferenceElement(reference, feType, i)
    basis = fe.get_basis_functions()
    code += "if constexpr (k =="+str(i)+"){"

    code += getCodeForEvaluation(basis, derivative = "Divdiv", symmetric = symmetric)

    code += "\n}"
  code += "\n}"

  print(code)

  print( "  }//namespace Impl\n}//namespace Dune\n#endif" ) # namespace braces

if __name__== "__main__":
  init_printing()

  ## Note: in order to match the indices of the Lagrange from symfem to Dune we needed to introduce a "Dune-variant" of the Lagrange
  ### The code for 1d and 2d added to symfem/elements/lagrange.py line 63 is
  #          elif variant == "Dune" and reference.name in ("interval", "triangle"):
    # if reference.name == "interval":
    #     for i in range(order +1):
    #         dim  = 0 if i == 0 or i == order else reference.tdim
    #         subEntityCount = 1 if i == order else 0
    #         dofs.append(PointEvaluation(reference, (sympy.Rational(i,order),), entity=(dim, subEntityCount)))
  # hhj = createGenericReferenceElement("triangle", "HHJ", 3, variant="Dune")
  # print(hhj.get_basis_function(0))
  # print(diff(hhj.get_basis_function(0)[0][1], x))
  # print(hhj.get_basis_function(0)[0][1].diff(x))


  printEvaluationCode("HellanHerrmannJohnsonReference", reference = "triangle", feType = "HHJ",minOrder  = 0, maxOrder= 6, symmetric = True, variant = "Dune")