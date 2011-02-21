cmake_minimum_required(VERSION 2.4)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/tools/cmake/")

# Version of library
SET(VERSION_MAJOR 2 CACHE INTERNAL "")
SET(VERSION_MINOR 0 CACHE INTERNAL "")
SET(VERSION_PATCH 0 CACHE INTERNAL "")

PROJECT(meteoio)
ADD_SUBDIRECTORY(meteoio)
ADD_SUBDIRECTORY(meteoio/plugins)
ADD_SUBDIRECTORY(meteoio/plugins/gsn)
ADD_SUBDIRECTORY(meteoio/meteolaws)
ADD_SUBDIRECTORY(meteoio/meteofilters)

# require proper c++
set(DEBUG_ARITHM ON CACHE BOOL "Use the debug flags -DDEBUG_ARITHM to catch arithmetic exceptions")
set(WARNINGS "-Wall")
set(PROFILING "-pg -fprofile-arcs")
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
  set(ARCH "-march=amdfam10 -fPIC")
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
set(PLUGIN_COSMOXMLIO OFF CACHE BOOL "Compilation COSMOXMLIO ON or OFF")
set(PROJ4 OFF CACHE BOOL "Use PROJ4 for the class MapProj ON or OFF")

set(LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

###########################################################
#for the install target
FILE(GLOB hfiles "${CMAKE_CURRENT_SOURCE_DIR}/meteoio/*.h")
FILE(GLOB hlawsfiles "${CMAKE_CURRENT_SOURCE_DIR}/meteoio/meteolaws/*.h")
FILE(GLOB hfilterfiles "${CMAKE_CURRENT_SOURCE_DIR}/meteoio/meteofilters/*.h")

if (DEST MATCHES "par")
  FILE(GLOB phfiles "${CMAKE_CURRENT_SOURCE_DIR}/meteoio/*.ph")
endif(DEST MATCHES "par")
INSTALL(FILES ${hfiles} ${phfiles} DESTINATION include/meteoio COMPONENT headers)
INSTALL(FILES ${hlawsfiles} DESTINATION include/meteoio/meteolaws COMPONENT headers)
INSTALL(FILES ${hfilterfiles} DESTINATION include/meteoio/meteofilters COMPONENT headers)

FILE(GLOB solibs "${PROJECT_SOURCE_DIR}/lib/libmeteoio*.so*")
FILE(GLOB plugins "${PROJECT_SOURCE_DIR}/lib/plugins/lib*.so*")
FILE(GLOB alibs "${PROJECT_SOURCE_DIR}/lib/*.a*")
if (DEST MATCHES "par")
  FILE(GLOB modules "${PROJECT_SOURCE_DIR}/lib/*.module")
  INSTALL(FILES ${modules} DESTINATION lib/meteoio COMPONENT libraries)
endif(DEST MATCHES "par")
INSTALL(FILES ${solibs} ${alibs} DESTINATION lib COMPONENT libraries)
INSTALL(FILES ${plugins} DESTINATION lib/meteoio/plugins COMPONENT plugins)


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

#Build the documentation if doxygen is found
include(FindDoxygen)
FIND_PACKAGE(Doxygen)
if (DOXYGEN_FOUND)
	ADD_CUSTOM_TARGET(documentation "doxygen" ${PROJECT_SOURCE_DIR}/meteoio/config.dox)
	get_target_property(DOC_TARGET doc TYPE)
	if(NOT DOC_TARGET)
		add_custom_target(doc)
	endif(NOT DOC_TARGET)
	add_dependencies(doc documentation)
endif (DOXYGEN_FOUND)
#Install the documentation if the user desires it
set(INSTALL_DOC OFF CACHE BOOL "Install documentation ON or OFF")
if(INSTALL_DOC)
	set(DOC_PATH "share/doc/${CMAKE_PROJECT_NAME}")
	if (DOXYGEN_FOUND)
		install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html DESTINATION ${DOC_PATH} COMPONENT documentation)
	else (DOXYGEN_FOUND)
		message("Doxygen has not been found and is needed to build online documentation. Please either install doxygen or download a documentation-only package!")
	endif (DOXYGEN_FOUND)
	set(DOC_FILES License.txt doc/SMET_specifications.pdf)
	install(FILES ${DOC_FILES} DESTINATION ${DOC_PATH} COMPONENT documentation)
	install(DIRECTORY "${PROJECT_SOURCE_DIR}/doc/examples"
		DESTINATION "share/doc/${CMAKE_PROJECT_NAME}"
		PATTERN ".svn" EXCLUDE
		COMPONENT examples)
endif(INSTALL_DOC)

#some cmake specific things
if (GSOAPPP AND NOT EXISTS ${GSOAPPP})
  set(GSOAPPP NOTFOUND CACHE FORCE)
endif(GSOAPPP AND NOT EXISTS ${GSOAPPP})

###########################################################
set(CPACK_PACKAGE_DESCRIPTION "Meteorological preprocessing library")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Meteorological Input/Output preprocessing library: it covers data retrieval, data filtering, data transformations and spatial interpolations. Please refer to the included documentation and the included code examples!")
set(CPACK_PACKAGE_NAME "MeteoIO")
set(CPACK_PACKAGE_CONTACT "Mathias Bavay ")
set(CPACK_PACKAGE_VENDOR "WSL Institute for Snow and Avalanche Research")
set(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/meteoio/icons\\\\icon.ico")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/License.txt")

if (DEST MATCHES "par")
  set(PACKAGE_EXTRA "-par")
endif(DEST MATCHES "par")
set(CPACK_PACKAGE_VERSION_MAJOR "${VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${VERSION_PATCH}")
set(VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
#set(CPACK_ARCH "${CMAKE_SYSTEM_PROCESSOR}")
set(CPACK_ARCH "i386")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}${PACKAGE_EXTRA}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-${CPACK_ARCH}")

#installation by components for formats that support it
#set(CPACK_COMPONENTS_ALL_GROUPS_IN_ONE_PACKAGE 1)
set(CPACK_COMPONENTS_ALL libraries headers plugins documentation examples)
set(CPACK_COMPONENT_libraries_DISPLAY_NAME "MeteoIO library")
set(CPACK_COMPONENT_libraries_DESCRIPTION "The MeteoIO library itself")
set(CPACK_COMPONENT_libraries_GROUP Runtime)
set(CPACK_COMPONENT_headers_DISPLAY_NAME "Header files")
set(CPACK_COMPONENT_headers_DESCRIPTION "Necessary C++ header files for developing applications using MeteoIO")
set(CPACK_COMPONENT_headers_GROUP Development)
set(CPACK_COMPONENT_headers_DEPENDS libraries)
set(CPACK_COMPONENT_plugins_DISPLAY_NAME "Input/Output Plugins")
set(CPACK_COMPONENT_plugins_DESCRIPTION "Various plugins for supporting various data formats/protocols")
set(CPACK_COMPONENT_plugins_GROUP Runtime)
set(CPACK_COMPONENT_plugins_DEPENDS libraries)
set(CPACK_COMPONENT_documentation_DISPLAY_NAME "Documentation")
set(CPACK_COMPONENT_documentation_DESCRIPTION "Documentation of the MeteoIO usage and API both for the end user and the developer using MeteoIO in his/her application")
set(CPACK_COMPONENT_examples_DISPLAY_NAME "Examples")
set(CPACK_COMPONENT_examples_DESCRIPTION "C++ Source code of various simple examples of applications with their example data")
set(CPACK_COMPONENT_examples_GROUP Development)
set(CPACK_COMPONENT_examples_DEPENDS libraries headers)

set(CPACK_COMPONENT_Runtime_GROUP "Runtime")
set(CPACK_COMPONENT_GROUP_Runtime_DESCRIPTION "What you need to run an application based on MeteoIO")
set(CPACK_COMPONENT_Development_GROUP "Development")
set(CPACK_COMPONENT_GROUP_Development_DESCRIPTION "What you need to develop an application based on MeteoIO")

#cpack_add_component(libraries DISPLAY_NAME "MeteoIO library" DESCRIPTION "The MeteoIO library itself" GROUP Runtime)
# cpack_add_component(headers DISPLAY_NAME "Header files" DESCRIPTION "Necessary C++ header files for developing applications using MeteoIO" DEPENDS libraries GROUP Development)
# cpack_add_component(plugins DISPLAY_NAME "Input/Output Plugins" DESCRIPTION "Various plugins for supporting various data formats/protocols" DEPENDS libraries GROUP Runtime)
# cpack_add_component(documentation DISPLAY_NAME "Documentation" DESCRIPTION "Documentation of the MeteoIO usage and API both for the end user and the developer using MeteoIO in his/her application")
# cpack_add_component(examples DISPLAY_NAME "Examples" DESCRIPTION "C++ Source code of various simple examples of applications with their example data" DEPENDS "libraries;headers" GROUP Development)

# cpack_add_component_group(Runtime EXPANDED DESCRIPTION "What you need to run an application based on MeteoIO")
# cpack_add_component_group(Development EXPANDED DESCRIPTION "What you need to develop an application based on MeteoIO")

IF(UNIX)
	IF(APPLE)
		set(CPACK_GENERATOR "PackageMaker")
		#set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/Welcome.txt")
	ELSE(APPLE)
		set(CPACK_GENERATOR "DEB;RPM;TGZ")
		set(CPACK_DEBIAN_PACKAGE_NAME "meteoio")
		set(CPACK_DEBIAN_PACKAGE_SECTION "devel")
		set(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6, libproj0")
		set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")
		set(CPACK_RPM_PACKAGE_NAME "meteoio")
		set(CPACK_RPM_PACKAGE_LICENSE "LPGLv3")
		set(CPACK_RPM_PACKAGE_REQUIRES "libstdc++6, libproj4")
		set(CPACK_RPM_PACKAGE_ARCHITECTURE "i386")
		set(CPACK_RPM_PACKAGE_GROUP "Development/Libraries")
	ENDIF(APPLE)
ELSE(UNIX)
	IF(WIN32)
		set(CPACK_GENERATOR "NSIS")
		set(CPACK_NSIS_MUI_ICON "${CPACK_PACKAGE_ICON}")
		set(CPACK_PACKAGE_ICON "${CPACK_PACKAGE_ICON}")
		set(CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}-${VERSION}")
		set(CPACK_NSIS_MENU_LINKS
			"doc/html/index.html" "MeteoIO documentation"
			"doc/examples" "MeteoIO examples"
			"https://slfsmm.indefero.net/p/meteoio/" "Home page")
		set(CPACK_NSIS_URL_INFO_ABOUT "https://slfsmm.indefero.net/p/meteoio/")
		set(CPACK_NSIS_HELP_LINK "doc/html/index.html")
		#set(CPACK_NSIS_PACKAGE_NAME "${EWS_APP_NAME}")
		#set(CPACK_NSIS_DISPLAY_NAME "${EWS_APP_NAME}")
		#set(CPACK_NSIS_CONTACT "${EWS_EMAIL}")
	ENDIF(WIN32)
ENDIF(UNIX)


set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${VERSION}" CACHE INTERNAL "tarball basename")
set(CPACK_SOURCE_GENERATOR TGZ)
# The following components are regex's to match anywhere (unless anchored)
# in absolute path + filename to find files or directories to be excluded
# from source tarball.
set(CPACK_SOURCE_IGNORE_FILES
	"~$"
	"\\\\.svn/"
	"/CMakeFiles/"
	"/CMakeCache*"
	"/cmake_install\\\\.cmake$"
	"/cmake_uninstall\\\\.cmake$"
	"^_CPack_.*/"
	"/CPack*"
	"\\\\.o$"
	"/m4/"
)

INCLUDE(CPack)
