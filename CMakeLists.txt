cmake_minimum_required(VERSION 2.6) 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/tools/cmake/")

set(CMAKE_C_COMPILER "colorgcc")
set(CMAKE_CXX_COMPILER "colorgcc")

PROJECT(meteoio CXX C)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(src/plugins)
ADD_SUBDIRECTORY(src/plugins/gsn)
ADD_SUBDIRECTORY(src/filter)

# require proper c++
set(CFLAGS "-Wall -ansi -pedantic -Wno-long-long")
set(POPCC popcc)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "release" CACHE STRING
      "Choose the type of build, options are: debug release "
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(NOT DEST)
  set(DEST "safe" CACHE STRING "Choose safe, zeus or grid" FORCE)
endif(NOT DEST)

if (DEST STREQUAL "safe")
  set(OPTIM "-O3 -fomit-frame-pointer -DNDEBUG")
  set(ARCH  "-march=pentium3")
  set(DEBUG "-g -O0 -D__DEBUG")
endif(DEST STREQUAL "safe")

if (DEST STREQUAL "grid")
  set(CMAKE_CXX_COMPILER "colorgcc" CACHE STRING 
    "Compiler set to colorgcc for the grid" FORCE)
  set(CMAKE_CC_COMPILER "colorgcc" CACHE STRING 
    "Compiler set to colorgcc for the grid" FORCE)
  set(OPTIM "-O3 -fomit-frame-pointer -DNDEBUG")
  set(ARCH  "-march=native")
  set(DEBUG "-g -O0 -D__DEBUG")
endif(DEST STREQUAL "grid")

if (DEST STREQUAL "zeus")
  set(CMAKE_CXX_COMPILER "g++" CACHE STRING 
    "Compiler set to g++ for zeus" FORCE)
  set(CMAKE_CC_COMPILER "gcc" CACHE STRING 
    "Compiler set to gcc for zeus" FORCE)
  set(OPTIM "-DNDEBUG")
  set(ARCH  "-march=x86-64")
  set(DEBUG "-g -O0 -D__DEBUG")
endif(DEST STREQUAL "zeus")

set(CMAKE_CXX_FLAGS_RELEASE "${ARCH} ${OPTIM} ${CFLAGS}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "${DEBUG} ${CFLAGS}" CACHE STRING "" FORCE)

message(STATUS "DEBUG CXXFLAGS  : ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "RELEASE CXXFLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")

set(ARCIO ON CACHE BOOL "Compilation ARCIO ON or OFF")
set(BOSCHUNGIO OFF CACHE BOOL "Compilation BOSCHUNGIO ON or OFF")
set(GEOTOPIO OFF CACHE BOOL "Compilation GEOTOPIO ON or OFF")
set(GRASSIO ON CACHE BOOL "Compilation GRASSIO ON or OFF")
set(GSNIO OFF CACHE BOOL "Compilation GSNIO ON or OFF")
set(IMISIO OFF CACHE BOOL "Compilation IMISIO ON or OFF")

set(LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

add_custom_target(distclean make clean
  COMMAND cmake -E remove CMakeCache.txt cmake_install.cmake src/cmake_install.cmake src/filter/cmake_install.cmake src/plugins/cmake_install.cmake src/plugins/gsn/cmake_install.cmake
  COMMAND cmake -E remove_directory CMakeFiles
  COMMAND cmake -E remove_directory src/CMakeFiles
  COMMAND cmake -E remove_directory src/filter/CMakeFiles
  COMMAND cmake -E remove_directory src/plugins/gsn/CMakeFiles
  COMMAND cmake -E remove_directory src/plugins/CMakeFiles)


#for the uninstall target
CONFIGURE_FILE(
  "${PROJECT_SOURCE_DIR}/tools/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")


#some cmake specific things
if (GSOAPPP AND NOT EXISTS ${GSOAPPP})
  set(GSOAPPP NOTFOUND CACHE FORCE)
endif(GSOAPPP AND NOT EXISTS ${GSOAPPP})

