add_subdirectory("test")

install(FILES
        concepts.hh
        istlcommunication.hh
        istlvectorbackend.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/backends)
