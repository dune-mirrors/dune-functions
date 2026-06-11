# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

install(FILES
        basisset.hh
        bindcontext.hh
        concepts.hh
        derivative.hh
        geometryderivative.hh
        localfiniteelement.hh
        pipeline.hh
        piola.hh
        referenceevaluation.hh
        simplexcontext.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases/transformed)
