cmake_minimum_required(VERSION 2.4) 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/tools/cmake/")

# Version of library
SET(VERSION_MAJOR 1 CACHE INTERNAL "")
SET(VERSION_MINOR 0 CACHE INTERNAL "")
SET(VERSION_PATCH 0 CACHE INTERNAL "")

PROJECT(meteoio)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(src/plugins)
ADD_SUBDIRECTORY(src/plugins/gsn)

FIND_PROGRAM(COLORGCC_EXECUTABLE colorgcc)
#IF(COLORGCC_EXECUTABLE)
#  set(CMAKE_C_COMPILER "colorgcc" CACHE STRING "Compiler set to colorgcc" FORCE)
#  set(CMAKE_CXX_COMPILER "colorgcc" CACHE STRING "Compiler set to colorgcc" FORCE)
#ENDIF(COLORGCC_EXECUTABLE)

# require proper c++
set(DEBUG_ARITHM ON CACHE BOOL "Use the debug flags -DDEBUG_ARITHM to catch arithmetic exceptions")
IF(DEBUG_ARITHM)
  set(CFLAGS "-Wall -Wextra -ansi -pedantic -Wno-long-long -DDEBUG_ARITHM")
ELSE(DEBUG_ARITHM)
  set(CFLAGS "-Wall -Wextra -ansi -pedantic -Wno-long-long")
ENDIF(DEBUG_ARITHM)

SET(POPCC popcc)

#platform specific
IF(UNIX)
  IF(APPLE)
    set(PLATFORM "-DAPPLE")
    set(CFLAGS "-Wall -Wextra -Wno-long-long")
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

if(NOT DEST)
  set(DEST "safe" CACHE STRING "Choose safe, par (with popcc),zeus or grid" FORCE)
endif(NOT DEST)

if (DEST STREQUAL "par")
  set(CFLAGS "-Wall -Wextra -Wno-long-long")
  set(OPTIM "-O3 -fomit-frame-pointer -DNDEBUG" CACHE STRING "" FORCE)
  set(DEBUG "-g -O0 -D__DEBUG" CACHE STRING "" FORCE)
endif(DEST STREQUAL "par")

if (DEST STREQUAL "safe")
  set(OPTIM "-O3 -fomit-frame-pointer -DNDEBUG")
  set(ARCH  "")
  set(DEBUG "-g -O0 -D__DEBUG")
endif(DEST STREQUAL "safe")

IF(DEST STREQUAL "grid")
  SET(OPTIM "-O3 -fomit-frame-pointer -DNDEBUG")
  SET(ARCH  "-march=native")
  SET(DEBUG "-g -O0 -D__DEBUG")
ENDIF(DEST STREQUAL "grid")

if (DEST STREQUAL "zeus")
  set(OPTIM "-DNDEBUG")
  set(ARCH  "-march=x86-64")
  set(DEBUG "-g -O0 -D__DEBUG")
endif(DEST STREQUAL "zeus")

set(CMAKE_CXX_FLAGS_RELEASE "${ARCH} ${OPTIM} ${CFLAGS} ${PLATFORM}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "${DEBUG} ${CFLAGS} ${PLATFORM}" CACHE STRING "" FORCE)

message(STATUS "DEBUG CXXFLAGS  : ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "RELEASE CXXFLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")

set(ARCIO ON CACHE BOOL "Compilation ARCIO ON or OFF")
set(ARPSIO ON CACHE BOOL "Compilation ARPSIO ON or OFF")
set(BORMAIO OFF CACHE BOOL "Compilation BORMAIO ON or OFF")
set(GEOTOPIO OFF CACHE BOOL "Compilation GEOTOPIO ON or OFF")
set(SNIO ON CACHE BOOL "Compilation SNIO ON or OFF")
set(GRASSIO ON CACHE BOOL "Compilation GRASSIO ON or OFF")
set(GSNIO OFF CACHE BOOL "Compilation GSNIO ON or OFF")
set(IMISIO OFF CACHE BOOL "Compilation IMISIO ON or OFF")
set(PROJ4 OFF CACHE BOOL "Use PROJ4 for the class MapProj ON or OFF")

set(LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

#for the install target
FILE(GLOB hfiles "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB phfiles "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${hfiles} ${phfiles} DESTINATION include)
FILE(GLOB solibs "${PROJECT_SOURCE_DIR}/lib/*.so")
FILE(GLOB alibs "${PROJECT_SOURCE_DIR}/lib/*.a")
FILE(GLOB modules "${PROJECT_SOURCE_DIR}/lib/*.module")
INSTALL(FILES ${solibs} ${alibs} ${modules} DESTINATION lib)

add_custom_target(distclean make clean
  COMMAND cmake -E remove CMakeCache.txt cmake_install.cmake src/cmake_install.cmake src/plugins/cmake_install.cmake src/plugins/gsn/cmake_install.cmake
  COMMAND cmake -E remove_directory CMakeFiles
  COMMAND cmake -E remove_directory src/CMakeFiles
  COMMAND cmake -E remove_directory src/plugins/gsn/CMakeFiles
  COMMAND cmake -E remove_directory src/plugins/CMakeFiles)

#for the uninstall target
CONFIGURE_FILE(
  "${PROJECT_SOURCE_DIR}/tools/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

ADD_CUSTOM_TARGET(documentation "doxygen" ${PROJECT_SOURCE_DIR}/src/config.dox)

#some cmake specific things
if (GSOAPPP AND NOT EXISTS ${GSOAPPP})
  set(GSOAPPP NOTFOUND CACHE FORCE)
endif(GSOAPPP AND NOT EXISTS ${GSOAPPP})
