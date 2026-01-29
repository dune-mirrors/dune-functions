# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

# By default we first search for a local version of MathJax and use it if found.
# If this is not found, we use the version provided by the cdn. Both versions
# can individually be disabled using the following two variables. Setting both
# will disable MathJax in Doxygen.
set(DUNE_FUNCTIONS_MATHJAX_DISABLE_LOCAL OFF CACHE BOOL "Flag to disable usage of local MathJax")
set(DUNE_FUNCTIONS_MATHJAX_DISABLE_CDN OFF CACHE BOOL "Flag to disable usage of cdn MathJax")

# These variables are set automatically and can be used in Doxylocal
# to set the respective Doxygen options.
set(DUNE_FUNCTIONS_MATHJAX_RELPATH "")
set(DUNE_FUNCTIONS_USE_MATHJAX "")

# Check for local MathJax and use it if found unless disabled.
if(NOT DUNE_FUNCTIONS_MATHJAX_DISABLE_LOCAL)
  # This currently searches for MathJax2 only which is the default in Doxygen.
  # Newer versions do not provide MathJax.js and have to be enabled manually
  # in Doxygen.
  find_path(DUNE_FUNCTIONS_MATHJAX_PATH
    MathJax.js
    PATHS
      "${CMAKE_PREFIX_PATH}/share/javascript/mathjax"
      "/usr/share/javascript/mathjax"
      "/usr/local/share/javascript/mathjax/"
      "/usr/share/yelp/mathjax"
  )
  if(DUNE_FUNCTIONS_MATHJAX_PATH MATCHES "DUNE_FUNCTIONS_MATHJAX_PATH-NOTFOUND")
    message(STATUS "Local MathJax not found")
  else()
    message(STATUS "Local MathJax found in ${DUNE_FUNCTIONS_MATHJAX_PATH}")
    set(DUNE_FUNCTIONS_USE_MATHJAX "YES")
    set(DUNE_FUNCTIONS_MATHJAX_RELPATH "${DUNE_FUNCTIONS_MATHJAX_PATH}")
  endif()
else()
  message(STATUS "MathJax was disabled")
endif()

# Fall back to cdn version of MathJax unless disabled.
# Note, that this is independent of the result of the local version check,
# because Doxygen uses the cdn version automatically if no path is given.
if(NOT DUNE_FUNCTIONS_MATHJAX_DISABLE_CDN)
  set(DUNE_FUNCTIONS_USE_MATHJAX "YES")
endif()



# Since the macros in dune-common do not allow to math custom variables
# the following is a copy of the macro from dune-common with added support
# for forwarding DUNE_FUNCTIONS_USE_MATHJAX and DUNE_FUNCTIONS_MATHJAX_RELPATH.
macro(prepare_doxyfile)
  cmake_parse_arguments(DOXYFILE "" "TARGET" "" ${ARGN})

  # default target name is the module name
  if(NOT DOXYFILE_TARGET)
    set(DOXYFILE_TARGET ${PROJECT_NAME})
  endif()

  message(STATUS "using ${DOXYSTYLE_FILE} to create doxystyle file")
  message(STATUS "using C macro definitions from ${DOXYGENMACROS_FILE} for Doxygen")

  # check whether module has a Doxylocal file
  find_file(_DOXYLOCAL Doxylocal PATHS ${CMAKE_CURRENT_SOURCE_DIR} NO_DEFAULT_PATH)

  if(_DOXYLOCAL)
    set(make_doxyfile_command ${CMAKE_COMMAND} -D DOT_TRUE=${DOT_TRUE} -D DUNE_MOD_NAME=${ProjectName} -D DUNE_MOD_VERSION=${ProjectVersion} -D DOXYSTYLE=${DOXYSTYLE_FILE} -D DOXYGENMACROS=${DOXYGENMACROS_FILE}  -D DOXYLOCAL=${CMAKE_CURRENT_SOURCE_DIR}/Doxylocal -D abs_top_srcdir=${CMAKE_SOURCE_DIR} -D srcdir=${CMAKE_CURRENT_SOURCE_DIR} -D top_srcdir=${${PROJECT_NAME}_SOURCE_DIR} -D DUNE_FUNCTIONS_USE_MATHJAX="${DUNE_FUNCTIONS_USE_MATHJAX}" -D DUNE_FUNCTIONS_MATHJAX_RELPATH="${DUNE_FUNCTIONS_MATHJAX_RELPATH}" -P ${scriptdir}/CreateDoxyFile.cmake)
    add_custom_command(OUTPUT Doxyfile.in Doxyfile
      COMMAND ${make_doxyfile_command}
      COMMENT "Creating Doxyfile.in"
      DEPENDS ${DOXYSTYLE_FILE} ${DOXYGENMACROS_FILE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxylocal)
  else()
    set(make_doxyfile_command ${CMAKE_COMMAND} -D DOT_TRUE=${DOT_TRUE} -D DUNE_MOD_NAME=${ProjectName} -D DUNE_MOD_VERSION=${DUNE_MOD_VERSION} -D DOXYSTYLE=${DOXYSTYLE_FILE} -D DOXYGENMACROS=${DOXYGENMACROS_FILE} -D abs_top_srcdir=${CMAKE_SOURCE_DIR} -D top_srcdir=${${PROJECT_NAME}_SOURCE_DIR} -P ${scriptdir}/CreateDoxyFile.cmake)
    add_custom_command(OUTPUT Doxyfile.in Doxyfile
      COMMAND ${make_doxyfile_command}
      COMMENT "Creating Doxyfile.in"
      DEPENDS ${DOXYSTYLE_FILE} ${DOXYGENMACROS_FILE})
  endif()
  add_custom_target(doxyfile_${DOXYFILE_TARGET} DEPENDS Doxyfile.in Doxyfile)
endmacro(prepare_doxyfile)

# shortcut for creating the Doxyfile.in and Doxyfile
add_doxygen_target()
