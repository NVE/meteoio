#OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/src")

SET(meteoio_sources 
	InterpolationAlgorithms.cc
	MeteoFilter.cc
	FilterAlgorithms.cc
  	A3DIO.cc
	Date_IO.cc         
	Grid2DObject.cc  
	IOHandler.cc    
	LegacyIO.cc       
	Coords.cc
	Meteo2DInterpolator.cc  
	BufferedIOHandler.cc  
	DEMObject.cc       
	Grid3DObject.cc  
	IOInterface.cc  
	libinterpol1D.cc  
	StationData.cc
	ConfigReader.cc       
	DynamicLibrary.cc  
	IOExceptions.cc  
	IOUtils.cc
	libinterpol2D.cc  
	MeteoData.cc
)

if (DEST STREQUAL "par")
  SET(CMAKE_CXX_SOURCE_FILE_EXTENSIONS cc;ph;)

  SET(meteoio_sources 
    marshal_meteoio.cc
    LegacyIO.ph
    IOHandler.ph
    ${meteoio_sources}
    )
  SET_SOURCE_FILES_PROPERTIES( ${meteoio_sources} PROPERTIES LANGUAGE CXX )
  SET_SOURCE_FILES_PROPERTIES( LegacyIO.ph PROPERTIES KEEP_EXTENSION TRUE)
  SET_SOURCE_FILES_PROPERTIES( IOHandler.ph PROPERTIES KEEP_EXTENSION TRUE)
  SET(STATICNAME ${PROJECT_NAME}_STATIC_par)
  SET(STATICLIBNAME ${PROJECT_NAME}popc)
  ADD_LIBRARY(${STATICNAME} STATIC ${meteoio_sources})
  #SET(SHAREDNAME ${PROJECT_NAME}_par)
  INCLUDE(CMakeForceCompiler)
  CMAKE_FORCE_C_COMPILER(popcc POPCC)
  CMAKE_FORCE_CXX_COMPILER(popcc POPCC)


  ADD_CUSTOM_COMMAND(
    OUTPUT "${PROJECT_SOURCE_DIR}/lib/PackMeteoIO.o"
    COMMAND popcc -o ${PROJECT_SOURCE_DIR}/lib/PackMeteoIO.o -c ${PROJECT_SOURCE_DIR}/src/PackMeteoIO.cc
    DEPENDS ${STATICNAME}
    )
  ADD_CUSTOM_COMMAND(
    OUTPUT ${PROJECT_SOURCE_DIR}/lib/meteoio.module
    COMMAND popcc "${CMAKE_CXX_FLAGS_RELEASE}" -object -parocld=g++ -o "${PROJECT_SOURCE_DIR}/lib/meteoio.module" -lmeteoiopopc -rdynamic -lc -ldl -lm -lstdc++
    DEPENDS ${STATICNAME} ${PROJECT_SOURCE_DIR}/lib/PackMeteoIO.o
    )

  ADD_CUSTOM_TARGET(weidstuff ALL DEPENDS ${PROJECT_SOURCE_DIR}/lib/meteoio.module)

  #$(LIBDIR)/meteoio.module: $(LIBDIR)/libmeteoiopopc.a $(SRCDIR)/PackMeteoIO_par.o
  #$(POPCC) $(CCFLAGS) -object -parocld=$(LINKER) -o $@ $(SRCDIR)/PackMeteoIO_par.o $(LDFLAGS) $(LDFLAGS_POPC)


ELSE(DEST)
  set(STATICNAME ${PROJECT_NAME}_STATIC)
  SET(STATICLIBNAME ${PROJECT_NAME})
  set(SHAREDNAME ${PROJECT_NAME})
  ADD_LIBRARY(${SHAREDNAME} ${meteoio_sources})
  ADD_LIBRARY(${STATICNAME} STATIC ${meteoio_sources})
  SET_TARGET_PROPERTIES(${SHAREDNAME} PROPERTIES 
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
    CLEAN_DIRECT_OUTPUT 1
    VERSION 0.1.0
    SOVERSION 0.1.0
    )

  INSTALL(TARGETS ${SHAREDNAME}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
endif(DEST STREQUAL "par")


SET_TARGET_PROPERTIES(${STATICNAME} PROPERTIES 
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
  CLEAN_DIRECT_OUTPUT 1
  OUTPUT_NAME ${STATICLIBNAME}
  VERSION 0.1.0
  SOVERSION 0.1.0
  )

#install targets:
INSTALL(TARGETS ${STATICNAME}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )

IF(PROJ4)
  FIND_PATH(LIBPROJ_INCLUDE_DIR "proj_api.h")
  FIND_LIBRARY(LIBPROJ proj)
  INCLUDE_DIRECTORIES(${LIBPROJ_INCLUDE_DIR})

  if (DEST STREQUAL "par")
    TARGET_LINK_LIBRARIES(${STATICNAME} ${LIBPROJ})
  else(DEST)
    TARGET_LINK_LIBRARIES(${SHAREDNAME} ${LIBPROJ})
  endif(DEST STREQUAL "par")
  ADD_DEFINITIONS(-DPROJ4)
ENDIF(PROJ4)

if (DEST STREQUAL "par")
  FILE(GLOB files "${PROJECT_SOURCE_DIR}/src/*.h ${PROJECT_SOURCE_DIR}/src/*.ph")
else(DEST)
  FILE(GLOB files "${PROJECT_SOURCE_DIR}/src/*.h")
endif(DEST STREQUAL "par")

INSTALL(FILES ${files} DESTINATION include)
