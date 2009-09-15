#OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/src")

SET(meteoio_sources 
	Laws.c
  	A3DIO.cc
	Date_IO.cc         
	Grid2DObject.cc  
	IOHandler.cc    
	LegacyIO.cc       
	MapProj.cc           
	Meteo2DInterpolator.cc  
	BufferedIOHandler.cc  
	DEMObject.cc       
	Grid3DObject.cc  
	IOInterface.cc  
	libinterpol1D.cc  
	MeteoBuffer.cc          
	StationData.cc
	ConfigReader.cc       
	DynamicLibrary.cc  
	IOExceptions.cc  
	IOUtils.cc
	libinterpol2D.cc  
	Meteo1DResampler.cc  
	MeteoData.cc
)

ADD_LIBRARY(${PROJECT_NAME} STATIC ${meteoio_sources})
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES 
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
  )

#install targets:
INSTALL(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

FILE(GLOB files "${PROJECT_SOURCE_DIR}/src/*.h")
INSTALL(FILES ${files} DESTINATION include)
