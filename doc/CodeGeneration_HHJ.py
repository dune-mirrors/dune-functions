# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

"""Generate the one- and two-dimensional HHJ reference bases from Symfem."""

import symfem
from symfem.functions import _to_sympy_format, parse_function_input
from sympy import Matrix, cxxcode, diff, init_printing, shape
from sympy.abc import x, y
from sympy.polys.polyfuncs import horner


coordinates = (x, y)


def create_reference_element(reference, order, **kwargs):
  """Create the scalar 1D analogue or the HHJ element with Dune ordering."""
  if reference == "interval":
    element = symfem.create_element(reference, "Lagrange", order, **kwargs)
  else:
    element = symfem.create_element(reference, "HHJ", order, **kwargs)
    assert element.reference.name == "triangle"
    assert element.reference.edges == ((0, 1), (0, 2), (1, 2))
  return element


def matrix_basis(element):
  """Represent the scalar interval basis as one-by-one matrices."""
  basis = element.get_basis_functions()
  if element.reference.name == "interval":
    return [Matrix([[parse_function_input(_to_sympy_format(function))]]) for function in basis]
  return basis


def apply_horner_scheme(function, derivative=None):
  """Apply a Horner scheme to a matrix function or its div-div."""
  function_shape = shape(function)
  assert len(function_shape) == 2 and function_shape[0] == function_shape[1]
  dim = function_shape[0]

  if derivative == "Divdiv":
    result = sum(
        diff(diff(function[i, j], coordinates[j]), coordinates[i])
        for i in range(dim)
        for j in range(dim)
    )
    return horner(result)

  return [
      [horner(function[i, j]) for j in range(dim)]
      for i in range(dim)
  ]


def code_for_list(value, **kwargs):
  if isinstance(value, list):
    return "{" + ", ".join(code_for_list(entry, **kwargs) for entry in value) + "}"
  return cxxcode(value, **kwargs)


def code_for_matrix(tensor, **kwargs):
  entries = [
      code_for_list(tensor[i][j], **kwargs)
      for i in range(len(tensor))
      for j in range(i, len(tensor))
  ]
  return "\n*(iter++) = sym<Range>(" + ", ".join(entries) + ");\n"


def code_for_basis(basis, derivative=None):
  power_substitutions = {
      "Pow": [
          (lambda base, exponent: exponent == 2, lambda base, exponent: f"{base}*{base}"),
          (lambda base, exponent: exponent == 3, lambda base, exponent: f"{base}*{base}*{base}"),
          (lambda base, exponent: exponent == 4, lambda base, exponent: f"{base}*{base}*{base}*{base}"),
          (lambda base, exponent: not base.is_integer, "pow"),
      ]
  }

  code = ""
  for index, function in enumerate(basis):
    code += f"\n// {index}th basis function"
    value = apply_horner_scheme(function, derivative)
    if derivative is None:
      code += code_for_matrix(value, user_functions=power_substitutions)
    else:
      code += f"\n*(iter++) = {code_for_list(value, user_functions=power_substitutions)};\n"
  return code


def generate_method(name, references, min_order, max_order, derivative=None, **kwargs):
  if derivative is None:
    method = "evaluateFunction"
    result_type = "RangeType"
  else:
    method = "evaluateDivDiv"
    result_type = "DivDivType"

  code = "template<class D, class R, int dim, unsigned int k>\n"
  code += f"void {name}LocalBasis<D, R, dim, k>::{method}(\n"
  code += "    const typename Traits::DomainType& in,\n"
  code += f"    std::vector<typename Traits::{result_type}>& out) const\n{{\n"
  code += "  out.resize(size());\n  auto iter = out.begin();\n"
  if derivative is None:
    code += "  using Range = typename Traits::RangeType;\n"
  code += "\n  static_assert(dim == 1 || dim == 2);\n"
  code += "\n  // Generated with SymPy from the Symfem library.\n"

  for reference in references:
    ref_dim = symfem.create_reference(reference).tdim
    coordinate_definitions = ", ".join(
        f"{coordinates[i]} = in[{i}]" for i in range(ref_dim)
    )
    code += f"  if constexpr (dim == {ref_dim}) {{\n"
    code += f"    auto const& {coordinate_definitions};\n"
    for order in range(min_order, max_order + 1):
      element = create_reference_element(reference, order, **kwargs)
      code += f"    if constexpr (k == {order}) {{\n"
      code += code_for_basis(matrix_basis(element), derivative)
      code += "    }\n"
    code += "  }\n"
  code += "}\n"
  return code


def print_evaluation_code(name, references, min_order=0, max_order=3, **kwargs):
  """Print the generated include containing value and div-div evaluation."""
  guard_name = name.upper().replace("REFERENCE", "BASIS")

  # Assemble the marker at runtime so REUSE does not treat the generated-file
  # annotations as additional annotations of this Python source file.
  spdx_marker = "SPDX"
  code = f"// {spdx_marker}-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md\n"
  code += f"// {spdx_marker}-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later\n\n"
  code += "// Generated by doc/CodeGeneration_HHJ.py; do not edit manually.\n\n"
  code += f"#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_{guard_name}_INC_HH\n"
  code += f"#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_{guard_name}_INC_HH\n\n"
  code += "namespace Dune::Functions\n{\n  namespace Impl\n  {\n"
  code += generate_method(name, references, min_order, max_order, **kwargs)
  code += generate_method(name, references, min_order, max_order, derivative="Divdiv", **kwargs)
  code += "  } // namespace Impl\n} // namespace Dune::Functions\n\n#endif\n"
  print(code)


if __name__ == "__main__":
  init_printing()
  print_evaluation_code(
      "HellanHerrmannJohnsonReference",
      references=("interval", "triangle"),
      min_order=0,
      max_order=6,
      variant="dune",
  )
