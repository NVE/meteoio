#OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")

IF(PROJ4)
  FIND_PATH(LIBPROJ_INCLUDE_DIR "proj_api.h")
  FIND_LIBRARY(LIBPROJ proj)
  INCLUDE_DIRECTORIES(${LIBPROJ_INCLUDE_DIR})
  ADD_DEFINITIONS(-DPROJ4)
ENDIF(PROJ4)

INCLUDE("${PROJECT_SOURCE_DIR}/meteoio/meteolaws/CMakeLists.txt")

SET(meteoio_sources
	IOManager.cc
	WindowedFilter.cc
	ProcessingStack.cc
	FilterMeanAvg.cc
	FilterMedianAvg.cc
	FilterBlock.cc
	FilterMinMax.cc
	ProcessingBlock.cc
	Matrix.cc
	IOPlugin.cc
	MeteoProcessor.cc
	ResamplingAlgorithms.cc
	Meteo1DInterpolator.cc
	InterpolationAlgorithms.cc
	MeteoFilter.cc
	FilterAlgorithms.cc
	A3DIO.cc
	Date.cc
	Timer.cc
	Grid2DObject.cc
	IOHandler.cc
	Coords.cc
	Meteo2DInterpolator.cc
	BufferedIOHandler.cc
	DEMObject.cc
	Grid3DObject.cc
	IOInterface.cc
	libinterpol1D.cc
	libfit1D.cc
	StationData.cc
	Config.cc
	DynamicLibrary.cc
	IOExceptions.cc
	IOUtils.cc
	libinterpol2D.cc
	MeteoData.cc
	${meteolaws_sources}
)

if (DEST MATCHES "par")
  #set_source_files_properties(SMETIO.cc PROPERTIES COMPILE_FLAGS "-fPIC")
  SET(CMAKE_CXX_SOURCE_FILE_EXTENSIONS cc;ph;)

  SET(meteoio_sources 
    marshal_meteoio.cc
    IOHandler.ph
    ${meteoio_sources}
    )
  SET_SOURCE_FILES_PROPERTIES( ${meteoio_sources} PROPERTIES LANGUAGE CXX )
  SET_SOURCE_FILES_PROPERTIES( IOHandler.ph PROPERTIES KEEP_EXTENSION TRUE)
  SET(STATICNAME ${PROJECT_NAME}_STATIC_par)
  SET(STATICLIBNAME ${PROJECT_NAME}popc)
  ADD_LIBRARY(${STATICNAME} STATIC ${meteoio_sources})
  TARGET_LINK_LIBRARIES(${STATICNAME} ${LIBPROJ})
  #SET(SHAREDNAME ${PROJECT_NAME}_par)
  INCLUDE(CMakeForceCompiler)
  #CMAKE_FORCE_C_COMPILER(popcc POPCC)
  CMAKE_FORCE_CXX_COMPILER(popcc POPCC)

  ADD_CUSTOM_COMMAND(
    OUTPUT "${PROJECT_SOURCE_DIR}/lib/PackMeteoIO.o"
    COMMAND popcc -o ${PROJECT_SOURCE_DIR}/lib/PackMeteoIO.o -c ${PROJECT_SOURCE_DIR}/meteoio/PackMeteoIO.cc -I${PROJECT_SOURCE_DIR} -combine -pipe "${ARCH} ${OPTIM} ${CFLAGS} ${PLATFORM} ${WARNINGS} ${EXTRA_WARNINGS}"
    DEPENDS ${STATICNAME}
    )
  ADD_CUSTOM_COMMAND(
    OUTPUT ${PROJECT_SOURCE_DIR}/lib/meteoio.module
    COMMAND popcc -object -parocld=c++ -o ${PROJECT_SOURCE_DIR}/lib/meteoio.module ${PROJECT_SOURCE_DIR}/lib/PackMeteoIO.o -L${PROJECT_SOURCE_DIR}/lib -lmeteoiopopc -rdynamic -lc -ldl -lm -lstdc++ ${LIBPROJ}
    DEPENDS ${STATICNAME} ${PROJECT_SOURCE_DIR}/lib/PackMeteoIO.o
    )

  ADD_CUSTOM_TARGET(weirdstuff ALL DEPENDS ${PROJECT_SOURCE_DIR}/lib/meteoio.module)

ELSE(DEST)
  set(STATICNAME ${PROJECT_NAME}_STATIC)
  SET(STATICLIBNAME ${PROJECT_NAME})
  set(SHAREDNAME ${PROJECT_NAME})
  ADD_LIBRARY(${SHAREDNAME} ${meteoio_sources})
  TARGET_LINK_LIBRARIES(${SHAREDNAME} ${LIBPROJ})
  ADD_LIBRARY(${STATICNAME} STATIC ${meteoio_sources})
  TARGET_LINK_LIBRARIES(${STATICNAME} ${LIBPROJ})
  SET_TARGET_PROPERTIES(${SHAREDNAME} PROPERTIES 
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
    CLEAN_DIRECT_OUTPUT 1
    VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}"
    SOVERSION "${VERSION_MAJOR}"
    )

  INSTALL(TARGETS ${SHAREDNAME}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
endif(DEST MATCHES "par")

SET_TARGET_PROPERTIES(${STATICNAME} PROPERTIES 
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
  CLEAN_DIRECT_OUTPUT 1
  OUTPUT_NAME ${STATICLIBNAME}
  )

INSTALL(TARGETS ${STATICNAME}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )
