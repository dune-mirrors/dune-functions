# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

add_subdirectory("test")

install(FILES
        analyticgridviewfunction.hh
        composedgridfunction.hh
        coarsefunctiononfinegridview.hh
        discreteglobalbasisfunction.hh
        finefunctiononcoarsegridview.hh
        gridfunction.hh
        gridfunction_imp.hh
        gridviewentityset.hh
        gridviewfunction.hh
        localderivativetraits.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/gridfunctions)
