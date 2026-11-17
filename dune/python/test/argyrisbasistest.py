# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

import dune.grid
import dune.functions

import basistest   # The general test suite for function space bases

mshfile = "../../dune-grid/doc/grids/gmsh/circle1storder.msh"
grid = dune.grid.ugGrid( (dune.grid.reader.gmsh, mshfile), dimgrid=2 )

basis = dune.functions.defaultGlobalBasis(grid, dune.functions.Argyris())
basistest.checkBasis(basis)
