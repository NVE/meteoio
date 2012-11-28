INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")

IF(PLUGIN_ARCIO)
	SET(plugins_sources ${plugins_sources} plugins/ARCIO.cc)
ENDIF(PLUGIN_ARCIO)

IF(PLUGIN_A3DIO)
	SET(plugins_sources ${plugins_sources} plugins/A3DIO.cc)
ENDIF(PLUGIN_A3DIO)

IF(PLUGIN_ARPSIO)
	SET(plugins_sources ${plugins_sources} plugins/ARPSIO.cc)
ENDIF(PLUGIN_ARPSIO)

IF(PLUGIN_BORMAIO)
	FIND_PACKAGE(LibXML++ REQUIRED)
	INCLUDE_DIRECTORIES(SYSTEM ${LibXML++_INCLUDE_DIRS})
	SET(plugin_libs ${plugin_libs} ${LibXML++_LIBRARIES})
	SET(plugins_sources ${plugins_sources} plugins/BormaIO.cc)
ENDIF(PLUGIN_BORMAIO)

IF(PLUGIN_COSMOXMLIO)
	FIND_PACKAGE(LibXML++ REQUIRED)
	INCLUDE_DIRECTORIES(SYSTEM ${LibXML++_INCLUDE_DIRS})
	SET(plugin_libs ${plugin_libs} ${LibXML++_LIBRARIES})
	SET(plugins_sources ${plugins_sources} plugins/CosmoXMLIO.cc)
ENDIF(PLUGIN_COSMOXMLIO)

IF(PLUGIN_GRASSIO)
	SET(plugins_sources ${plugins_sources} plugins/GrassIO.cc)
ENDIF(PLUGIN_GRASSIO)

IF(PLUGIN_GRIBIO)
	FIND_PACKAGE(GribApi REQUIRED)
	INCLUDE_DIRECTORIES(SYSTEM ${GribApi_INCLUDE_DIRS})
	SET(plugin_libs ${plugin_libs} ${GribApi_LIBRARIES})
	SET(plugins_sources ${plugins_sources} plugins/GRIBIO.cc)
ENDIF(PLUGIN_GRIBIO)

IF(PLUGIN_GEOTOPIO)
	SET(plugins_sources ${plugins_sources} plugins/GeotopIO.cc)
ENDIF(PLUGIN_GEOTOPIO)

IF(PLUGIN_GSNIO)
	#HACK
	#gsoap generates files that do lots of warnings, there is nothing we can do about it -> off
	IF(NOT MSVC)
		SET_SOURCE_FILES_PROPERTIES(plugins/gsn/soapC.cpp PROPERTIES COMPILE_FLAGS "-w")
		SET_SOURCE_FILES_PROPERTIES(plugins/gsn/stdsoap2.cpp PROPERTIES COMPILE_FLAGS "-w")
		SET_SOURCE_FILES_PROPERTIES(plugins/gsn/soapGSNWebServiceSoap12BindingProxy.cpp PROPERTIES COMPILE_FLAGS "-w")
		SET_SOURCE_FILES_PROPERTIES(plugins/GSNIO.cc PROPERTIES COMPILE_FLAGS "-Wno-long-long")
	ENDIF(NOT MSVC)
	IF(MSVC)
		SET_SOURCE_FILES_PROPERTIES(plugins/gsn/soapC.cpp PROPERTIES COMPILE_FLAGS "/W0")
		SET_SOURCE_FILES_PROPERTIES(plugins/gsn/stdsoap2.cpp PROPERTIES COMPILE_FLAGS "/W0")
		SET_SOURCE_FILES_PROPERTIES(plugins/gsn/soapGSNWebServiceSoap12BindingProxy.cpp PROPERTIES COMPILE_FLAGS "/W0")
	ENDIF(MSVC)

	IF(WIN32)
		IF(MINGW)
			ADD_DEFINITIONS(-D_isnan=isnan)
		ENDIF(MINGW)
		SET(EXTRA_LINK_LIBS "ws2_32")
	ENDIF(WIN32)

	SET(plugins_sources ${plugins_sources} plugins/GSNIO.cc plugins/gsn/soapGSNWebServiceSoap12BindingProxy.cpp plugins/gsn/soapC.cpp plugins/gsn/stdsoap2.cpp)
ENDIF(PLUGIN_GSNIO)

IF(PLUGIN_IMISIO)
	FIND_PACKAGE(OCCI REQUIRED)
	INCLUDE_DIRECTORIES(SYSTEM ${OCCI_INCLUDE_DIR})
	SET(plugin_libs ${plugin_libs} ${OCCI_LIBRARIES})
	SET(plugins_sources ${plugins_sources} plugins/ImisIO.cc)
ENDIF(PLUGIN_IMISIO)

IF(PLUGIN_PGMIO)
	SET(plugins_sources ${plugins_sources} plugins/PGMIO.cc)
ENDIF(PLUGIN_PGMIO)

IF(PLUGIN_PNGIO)
	FIND_PACKAGE(PNG REQUIRED)
	INCLUDE_DIRECTORIES(SYSTEM ${PNG_INCLUDE_DIR})
	SET(plugin_libs ${plugin_libs} ${PNG_LIBRARY})
	SET(plugins_sources ${plugins_sources} plugins/PNGIO.cc)
ENDIF(PLUGIN_PNGIO)

IF(PLUGIN_SMETIO)
	SET(plugins_sources ${plugins_sources} plugins/SMETIO.cc)
ENDIF(PLUGIN_SMETIO)

IF(PLUGIN_SNIO)
	SET(plugins_sources ${plugins_sources} plugins/SNIO.cc)
ENDIF(PLUGIN_SNIO)

CONFIGURE_FILE( ${PROJECT_SOURCE_DIR}/meteoio/IOHandler.cmake.cc ${PROJECT_SOURCE_DIR}/meteoio/IOHandler.cc )

