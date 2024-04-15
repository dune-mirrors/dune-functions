add_subdirectory("test")

install(FILES
        monomialset.hh
        polynomial.hh
        trigonometricfunction.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/analyticfunctions)
