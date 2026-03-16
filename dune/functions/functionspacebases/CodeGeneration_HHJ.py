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
  newFe = type(fe)(newRef, fe.order)
  return newFe

def adaptReferenceToPhysicalElement(fe, vertices):
  physRef = symfem.create_reference("triangle", vertices= vertices)
  assert(physRef.name == "triangle")
  # physRef.edges = ((0,1),(0,2),(1,2))
  newFe = type(fe)(physRef, fe.order)
  return newFe

def createGenericReferenceElement(refName, feName, order):
  return(adaptReferenceElementToDune(symfem.create_element(refName, feName, order)))

def createPhysicalElement(refName, feName, order, vertices):
  return(adaptReferenceToPhysicalElement(symfem.create_element(refName, feName, order), vertices))

## apply a horner scheme on a function f
def hornerScheme(f, derivative = [x,y], **kwargs):
  s = shape(f)
  assert(len(s) == 2)

  result =  [] #
  for i in range(s[0]):
    result.append([])
    for j in range(s[1]):
      if derivative is None: #values
        result[i].append(horner(f[i,j]))
      elif derivative == "Divdiv":
        result[i].append(horner(diff(f[i,j].diff(x), x)+ diff(f[i,j].diff(y), y)))
      elif isinstance(derivative, list):  #multiple derivatives, typically jacobian
        result[i].append([])
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
  symmetric = kwargs.get("symmetric")
  derivative = kwargs.pop("derivative", None)
  powsubs={'Pow': [(lambda b,e: e == 2, lambda b, e: ("{0}*{0}".format(b))),
  (lambda b,e: e == 3, lambda b, e: ("{0}*{0}*{0}".format(b))),
  (lambda b,e: e == 4, lambda b, e: ("{0}*{0}*{0}*{0}".format(b))),
  (lambda b, e: not b.is_integer, 'pow')]}
  code = "// generated with sympy from symfem library\n"
  code += "auto const&x = in[0], y = in[1];"

  kwargs.update( {"user_functions": powsubs})
  name = "val" if kwargs.get("assign_to") is None else kwargs.get("assign_to").name
  # kwargs.pop("assign_to")
  for i,f in enumerate(basis):
    code += "\n//{}th basis function\n".format(i+1)

    mat = hornerScheme(f, derivative)
    code += getCodeForMatrix(mat,**kwargs)
  return code

## Generate an include file for evaluation methods
def printEvaluationCode(name, reference, feType, maxOrder, symmetric = False):

  assert(isinstance(name, str))
  code = "#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_" + name.upper() + "_INC_HH\n#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_" + name.upper() + "_INC_HH\n namespace Dune{\n  namespace Impl{ \n    "
  code += "// generated with sympy from symfem library\ntemplate<class D, class R>\n"
  code +='     void ' + name + 'LocalBasis<D,R>::evaluateFunction(const typename Traits::DomainType &in,std::vector<typename Traits::RangeType> &out) const\n{\nout.resize(size());\n auto iter = out.begin();'

  for i in range(maxOrder):
    fe = createGenericReferenceElement(reference, feType, i)
    basis = fe.get_basis_functions()
    code += "if constexpr (k =="+str(i)+"){"
    code += getCodeForEvaluation(basis, symmetric = symmetric)

    code += "\n}"

  code +="\n}"
  code += '// generated with sympy from symfem library\ntemplate<class D, class R>\n      void ' + name + 'LocalBasis<D,R>::evaluateDivDiv(const typename Traits::DomainType &in,std::vector<typename Traits::DivDivType> &out) const\n{\nout.resize(size());\nauto iter = out.begin();'
  for i in range(maxOrder):
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

  printEvaluationCode("HellanHermannJohnsonReference", reference =  "triangle", feType = "HHJ",maxOrder= 8, symmetric = True)