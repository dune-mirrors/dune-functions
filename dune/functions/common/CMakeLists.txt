# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

add_subdirectory("test")

install(FILES
        childextraction.hh
        compositenode.hh
        defaultderivativetraits.hh
        densevectorview.hh
        differentiablefunction.hh
        differentiablefunction_imp.hh
        differentiablefunctionfromcallables.hh
        dynamicpowernode.hh
        functionconcepts.hh
        geometryinancestor.hh
        indexaccess.hh
        interfaces.hh
        leafnode.hh
        localfunction.hh
        localfunction_imp.hh
        mapperutilities.hh
        multiindex.hh
        nodeinterface.hh
        nodetags.hh
        overflowarray.hh
        polymorphicsmallobject.hh
        powernode.hh
        reserveddeque.hh
        signature.hh
        staticforloop.hh
        traversal.hh
        treecontainer.hh
        treepath.hh
        type_traits.hh
        typetraits.hh
        typeerasure.hh
        typetreeutility.hh
        utility.hh
        visitor.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/common)
