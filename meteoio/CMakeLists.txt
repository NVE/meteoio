INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")

SET(LIBPREFIX "lib") #make sure ALL plugins are libXXX for ALL plateforms
SET(CMAKE_IMPORT_LIBRARY_PREFIX "lib") #for Visual C++ import libraries
OPTION( BUILD_SHARED_LIBS "Set to ON to build dynamic libraries" ON )
IF(WIN32)
	#otherwise the lib files are not generated
	SET(BUILD_STATIC_LIBS ON)
ELSE(WIN32)
	OPTION( BUILD_STATIC_LIBS "Set to ON to build static libraries" OFF )
ENDIF(WIN32)

IF(PROJ4)
	FIND_PATH(LIBPROJ_INCLUDE_DIR "proj_api.h")
	FIND_LIBRARY(LIBPROJ proj)
	INCLUDE_DIRECTORIES(${LIBPROJ_INCLUDE_DIR})
	ADD_DEFINITIONS(-DPROJ4) #it looks like some VC++ versions don't support -D syntax
ENDIF(PROJ4)

INCLUDE("${PROJECT_SOURCE_DIR}/meteoio/meteolaws/CMakeLists.txt")
INCLUDE("${PROJECT_SOURCE_DIR}/meteoio/meteostats/CMakeLists.txt")
INCLUDE("${PROJECT_SOURCE_DIR}/meteoio/meteofilters/CMakeLists.txt")

IF(POPC)
	FIND_PACKAGE(Popc REQUIRED)
	INCLUDE_DIRECTORIES(${Popc_INCLUDE_DIRS})
	#TARGET_LINK_LIBRARIES(${STATICNAME} ${Popc_LIBRARIES} ${SHAREDNAME})
	SET(popc_sources marshal_meteoio.cc)
	#ADD_DEFINITIONS(-D_POPC_) #necessary if NOT forcing c++ compiler to popcc
	SET(CMAKE_CXX_COMPILER popcc)
	SET(CMAKE_C_COMPILER popcc)
	#INCLUDE(CMakeForceCompiler)
	#CMAKE_FORCE_C_COMPILER(popcc POPCC)
	#CMAKE_FORCE_CXX_COMPILER(popcc POPCC)
	SET(POPC_EXT "popc")
ENDIF(POPC)

#	FilterAlgorithms.cc obsoleted
SET(meteoio_sources
	IOManager.cc
	libsmet.cc
	Matrix.cc
	IOPlugin.cc
	MeteoProcessor.cc
	ResamplingAlgorithms.cc
	ResamplingAlgorithms2D.cc
	Meteo1DInterpolator.cc
	InterpolationAlgorithms.cc
	A3DIO.cc
	Date.cc
	Timer.cc
	Grid2DObject.cc
	IOHandler.cc
	Coords.cc
	Graphics.cc
	Meteo2DInterpolator.cc
	BufferedIOHandler.cc
	DEMObject.cc
	Grid3DObject.cc
	IOInterface.cc
	StationData.cc
	Config.cc
	DynamicLibrary.cc
	IOExceptions.cc
	IOUtils.cc
	MeteoData.cc
	${meteolaws_sources}
	${meteostats_sources}
	${meteofilters_sources}
	${popc_sources}
)

IF(BUILD_STATIC_LIBS)
	SET(STATICNAME ${PROJECT_NAME}_STATIC)
	SET(STATICLIBNAME ${PROJECT_NAME}${POPC_EXT})
	ADD_LIBRARY(${STATICNAME} STATIC ${meteoio_sources})
	TARGET_LINK_LIBRARIES(${STATICNAME} ${LIBPROJ} ${Popc_LIBRARIES})
ENDIF(BUILD_STATIC_LIBS)

SET(SHAREDNAME ${PROJECT_NAME}${POPC_EXT})
ADD_LIBRARY(${SHAREDNAME} ${meteoio_sources})
TARGET_LINK_LIBRARIES(${SHAREDNAME} ${LIBPROJ} ${Popc_LIBRARIES})
SET_TARGET_PROPERTIES(${SHAREDNAME} PROPERTIES
	PREFIX "${LIBPREFIX}"
	LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
	ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
	RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
	CLEAN_DIRECT_OUTPUT 1
	VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}"
	SOVERSION "${VERSION_MAJOR}"
)
INSTALL(TARGETS ${SHAREDNAME}
	RUNTIME DESTINATION lib
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	COMPONENT libraries
)

#static library. Please keep in mind that for Windows, a dll is a runtime
IF(BUILD_STATIC_LIBS)
	SET_TARGET_PROPERTIES(${STATICNAME} PROPERTIES
		PREFIX "${LIBPREFIX}"
		LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
		ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
		RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
		CLEAN_DIRECT_OUTPUT 1
		OUTPUT_NAME ${STATICLIBNAME}
	)
	INSTALL(TARGETS ${STATICNAME}
		RUNTIME DESTINATION lib
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
		COMPONENT libraries
	)
ENDIF(BUILD_STATIC_LIBS)
