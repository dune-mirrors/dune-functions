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
