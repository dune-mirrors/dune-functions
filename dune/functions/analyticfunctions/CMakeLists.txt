# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

if(DUNE_ENABLE_TESTING)
  add_subdirectory("test")
endif()

install(FILES
        monomialset.hh
        polynomial.hh
        trigonometricfunction.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/analyticfunctions)
