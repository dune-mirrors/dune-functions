add_subdirectory("test")

install(FILES
        entityowner.hh
        globalgridview.hh
        globalindexset.hh
        subpartitionprovider.hh
        uniqueborderpartition.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/functions/parallel)
