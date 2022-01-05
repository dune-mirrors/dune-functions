add_subdirectory(functions)
add_subdirectory(typetree2)

if( DUNE_ENABLE_PYTHONBINDINGS )
  add_subdirectory("python")
endif()
