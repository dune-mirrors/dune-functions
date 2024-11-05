# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

from .tree import Composite, DG, Lagrange, Nedelec, Power, RaviartThomas, Tree

duneFunctionsLayouts = {"lexicographic": "Lexicographic", "interleaved": "Interleaved"}

def indexMergingStrategy(blocked, layout):
    return "Dune::Functions::BasisFactory::" + ("Blocked" if blocked else "Flat") + duneFunctionsLayouts[layout]


def preBasisTypeName(tree, gridViewTypeName):
    assert isinstance(tree, Tree)
    if isinstance(tree, Lagrange):
        return "Dune::Functions::LagrangePreBasis< " + gridViewTypeName + " , " + str(tree.order) + " >"
    elif isinstance(tree, DG):
        raise Exception(repr(tree) + " not supported by dune-functions.")
    elif isinstance(tree, Nedelec):
        return "Dune::Functions::NedelecPreBasis< " + gridViewTypeName + ", double, " + str(tree.kind) + ", " + str(tree.order) + " >"
    elif isinstance(tree, RaviartThomas):
        return "Dune::Functions::RaviartThomasPreBasis< " + gridViewTypeName + ", " + str(tree.order) + " >"
    elif isinstance(tree, Composite):
        IMS = indexMergingStrategy(tree.blocked, tree.layout)
        ChildPreBases = " , ".join(preBasisTypeName(c, gridViewTypeName) for c in tree.children)
        return "Dune::Functions::CompositePreBasis< " + gridViewTypeName + ", " + IMS + " , " + ChildPreBases + " >"
    elif isinstance(tree, Power):
        IMS = indexMergingStrategy(tree.blocked, tree.layout)
        ChildPreBasis = preBasisTypeName(tree.children[0], gridViewTypeName)
        return "Dune::Functions::PowerPreBasis< " + IMS + " , " + ChildPreBasis + " , " + str(tree.exponent) + " >"
    else:
        raise Exception("Unknown type of tree: " + repr(tree))


def defaultGlobalBasis(gridView, tree):
    from dune.functions import load

    headers = ["powerbasis", "compositebasis", "lagrangebasis", "nedelecbasis", "raviartthomasbasis", "subspacebasis", "defaultglobalbasis"]

    includes = []
    includes += list(gridView.cppIncludes)
    includes += ["dune/functions/functionspacebases/" + h + ".hh" for h in headers]

    typeName = "Dune::Functions::DefaultGlobalBasis< " + preBasisTypeName(tree, gridView.cppTypeName) + " >"

    return load(includes, typeName).GlobalBasis(gridView)
