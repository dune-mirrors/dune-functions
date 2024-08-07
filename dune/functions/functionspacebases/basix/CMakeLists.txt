# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

install(FILES
  prebasis.hh
  finiteelement.hh
  localfiniteelement.hh
  utility.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/functionspacebases/basix)
