cmake_minimum_required(VERSION 2.4) 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/tools/cmake/")

# Version of library
SET(VERSION_MAJOR 1 CACHE INTERNAL "")
SET(VERSION_MINOR 1 CACHE INTERNAL "")
SET(VERSION_PATCH 0 CACHE INTERNAL "")

PROJECT(meteoio)
ADD_SUBDIRECTORY(meteoio)
ADD_SUBDIRECTORY(meteoio/plugins)
ADD_SUBDIRECTORY(meteoio/plugins/gsn)

# require proper c++
set(DEBUG_ARITHM ON CACHE BOOL "Use the debug flags -DDEBUG_ARITHM to catch arithmetic exceptions")
set(WARNINGS "-Wall")
set(EXTRA_WARNINGS "-Wextra -ansi -pedantic")

IF(DEBUG_ARITHM)
  LIST(APPEND CFLAGS "-DDEBUG_ARITHM")
ENDIF(DEBUG_ARITHM)

SET(POPCC popcc)

###########################################################
#platform specific
IF(UNIX)
  IF(APPLE)
    set(PLATFORM "-DAPPLE")
  ELSE(APPLE)
    set(PLATFORM "-DLINUX")
  ENDIF(APPLE)
ELSE(UNIX)
  IF(WIN32)
    set(PLATFORM "-DWIN32 -U__STRICT_ANSI__")
  ENDIF(WIN32)
ENDIF(UNIX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "release" CACHE STRING
      "Choose the type of build, options are: debug release "
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

###########################################################
#targets providing sets of compiler options
if(NOT DEST)
  set(DEST "safe" CACHE STRING "Choose safe, optimized, hera, or par, hera_par (with popc)" FORCE)
endif(NOT DEST)

if (DEST STREQUAL "par")
  set(EXTRA_WARNINGS "-Wextra ")
  set(OPTIM "-g -O3 -DNDEBUG")
  set(ARCH  "-mtune=native")
  set(DEBUG "-g -O0 -D__DEBUG")
endif(DEST STREQUAL "par")

if (DEST STREQUAL "hera_par")
  set(EXTRA_WARNINGS "-Wextra ")
  set(OPTIM "-g -O3 -DNDEBUG")
  #required for AMD64: amdfam10 or fPIC, even for static libs
  set(ARCH "-march=amdfam10")
  set(DEBUG "-g -O0 -D__DEBUG")
endif(DEST STREQUAL "hera_par")

if (DEST STREQUAL "safe")
  set(OPTIM "-g -O3 -DNDEBUG")
  set(ARCH  "")
  set(DEBUG "-g -O0 -D__DEBUG")
endif(DEST STREQUAL "safe")

IF(DEST STREQUAL "optimized")
  SET(OPTIM "-g -O3 -DNDEBUG")
  SET(ARCH  "-march=pentium4 -mtune=native")
  SET(DEBUG "-g -O0 -D__DEBUG")
ENDIF(DEST STREQUAL "optimized")

if (DEST STREQUAL "hera")
  set(OPTIM "-g -O3 -DNDEBUG")
  set(ARCH  "-march=amdfam10")
  set(DEBUG "-g -O0 -D__DEBUG")
endif(DEST STREQUAL "hera")

set(CMAKE_CXX_FLAGS_RELEASE "${ARCH} ${OPTIM} ${CFLAGS} ${PLATFORM}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "${DEBUG} ${ARCH} ${CFLAGS} ${PLATFORM} ${WARNINGS} ${EXTRA_WARNINGS}" CACHE STRING "" FORCE)

message(STATUS "DEBUG CXXFLAGS  : ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "RELEASE CXXFLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")

###########################################################
#plugins
set(PLUGIN_ARCIO ON CACHE BOOL "Compilation ARCIO ON or OFF")
set(PLUGIN_ARPSIO ON CACHE BOOL "Compilation ARPSIO ON or OFF")
set(PLUGIN_PGMIO ON CACHE BOOL "Compilation PGMIO ON or OFF")
set(PLUGIN_BORMAIO OFF CACHE BOOL "Compilation BORMAIO ON or OFF")
set(PLUGIN_GEOTOPIO ON CACHE BOOL "Compilation GEOTOPIO ON or OFF")
set(PLUGIN_SNIO ON CACHE BOOL "Compilation SNIO ON or OFF")
set(PLUGIN_SMETIO ON CACHE BOOL "Compilation SMETIO ON or OFF")
set(PLUGIN_GRASSIO ON CACHE BOOL "Compilation GRASSIO ON or OFF")
set(PLUGIN_GSNIO OFF CACHE BOOL "Compilation GSNIO ON or OFF")
set(PLUGIN_IMISIO OFF CACHE BOOL "Compilation IMISIO ON or OFF")
set(PROJ4 OFF CACHE BOOL "Use PROJ4 for the class MapProj ON or OFF")

set(LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

###########################################################
#for the install target
FILE(GLOB hfiles "${CMAKE_CURRENT_SOURCE_DIR}/meteoio/*.h")
if (DEST MATCHES "par")
  FILE(GLOB phfiles "${CMAKE_CURRENT_SOURCE_DIR}/meteoio/*.ph")
endif(DEST MATCHES "par")
INSTALL(FILES ${hfiles} ${phfiles} DESTINATION include/meteoio)

FILE(GLOB solibs "${PROJECT_SOURCE_DIR}/lib/libmeteoio*.so")
FILE(GLOB plugins "${PROJECT_SOURCE_DIR}/lib/lib*.so")
FILE(GLOB alibs "${PROJECT_SOURCE_DIR}/lib/*.a")
if (DEST MATCHES "par")
  FILE(GLOB modules "${PROJECT_SOURCE_DIR}/lib/*.module")
endif(DEST MATCHES "par")
INSTALL(FILES ${solibs} ${alibs} ${modules} DESTINATION lib)
INSTALL(FILES ${plugins} DESTINATION lib/meteoio/plugins)


add_custom_target(distclean make clean
  COMMAND cmake -E remove CMakeCache.txt cmake_install.cmake meteoio/cmake_install.cmake meteoio/plugins/cmake_install.cmake meteoio/plugins/gsn/cmake_install.cmake
  COMMAND cmake -E remove_directory CMakeFiles
  COMMAND cmake -E remove_directory meteoio/CMakeFiles
  COMMAND cmake -E remove_directory meteoio/plugins/gsn/CMakeFiles
  COMMAND cmake -E remove_directory meteoio/plugins/CMakeFiles)

###########################################################
#for the uninstall target
CONFIGURE_FILE(
  "${PROJECT_SOURCE_DIR}/tools/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

#FindDoxygen
include(FindDoxygen)
FIND_PACKAGE(Doxygen)
if (DOXYGEN_FOUND)
  ADD_CUSTOM_TARGET(documentation "doxygen" ${PROJECT_SOURCE_DIR}/meteoio/config.dox)
  get_target_property(DOC_TARGET doc TYPE)
  if(NOT DOC_TARGET)
	add_custom_target(doc)
  endif(NOT DOC_TARGET)
  add_dependencies(doc documentation)
  set(DOC_FILES License.txt doc/SMET_specifications.pdf)
  set(DOC_PATH "share/doc/${CMAKE_PROJECT_NAME}")
  install(FILES ${DOC_FILES} DESTINATION ${DOC_PATH})
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html DESTINATION ${DOC_PATH})
endif (DOXYGEN_FOUND)

#some cmake specific things
if (GSOAPPP AND NOT EXISTS ${GSOAPPP})
  set(GSOAPPP NOTFOUND CACHE FORCE)
endif(GSOAPPP AND NOT EXISTS ${GSOAPPP})

###########################################################
#creation of packages NOT FULLY WORKING YET
set(CPACK_PACKAGE_DESCRIPTION "Meteorological preprocessing library")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Meteorological Input/Output preprocessing library: it covers data retrieval, data filtering, data transformations and spatial interpolations.")
set(CPACK_PACKAGE_NAME "MeteoIO")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6")

if (DEST MATCHES "par")
  set(PACKAGE_EXTRA "-par")
endif(DEST MATCHES "par")
set(CPACK_PACKAGE_CONTACT "Mathias Bavay ")
set(CPACK_PACKAGE_VENDOR "WSL/SLF")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

set(CPACK_GENERATOR "DEB;RPM;")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}${PACKAGE_EXTRA}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}-${CMAKE_SYSTEM_PROCESSOR}")

include(CPack)
