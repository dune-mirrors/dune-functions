# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

install(FILES
  prebasis.hh
  finiteelement.hh
  localfiniteelement.hh
  utility.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases/basix)
