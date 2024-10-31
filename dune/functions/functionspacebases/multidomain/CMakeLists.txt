# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#install headers
install(FILES
  subdomainbasis.hh
  subdomaingridview.hh
  subdomainindexset.hh
  domaininfo.hh
  multidomainnode.hh
  #surfacebasis.hh
  #surfacegridview.hh
  #surfaceindexset.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/subdomain)
