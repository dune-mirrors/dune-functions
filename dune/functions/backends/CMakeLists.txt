add_subdirectory("test")

install(FILES
        concepts.hh
        istlmatrixfactory.hh
        istlvectorbackend.hh
        singlerowcolmatrix.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/backends)
