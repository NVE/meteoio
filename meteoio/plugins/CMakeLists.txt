INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")

IF(PLUGIN_SMETIO)
	IF(MSVC)
		ADD_DEFINITIONS(/DPLUGIN_SMETIO) #it looks like some VC++ versions don't support -D syntax
	ELSE(MSVC)
		ADD_DEFINITIONS(-DPLUGIN_SMETIO)
	ENDIF(MSVC)
	SET(plugins_sources ${plugins_sources} plugins/SMETIO.cc)
ENDIF(PLUGIN_SMETIO)

IF(PLUGIN_ARCIO)
	IF(MSVC)
		ADD_DEFINITIONS(/DPLUGIN_ARCIO) #it looks like some VC++ versions don't support -D syntax
	ELSE(MSVC)
		ADD_DEFINITIONS(-DPLUGIN_ARCIO)
	ENDIF(MSVC)
	SET(plugins_sources ${plugins_sources} plugins/ARCIO.cc)
ENDIF(PLUGIN_ARCIO)

IF(PLUGIN_A3DIO)
	IF(MSVC)
		ADD_DEFINITIONS(/DPLUGIN_A3DIO) #it looks like some VC++ versions don't support -D syntax
	ELSE(MSVC)
		ADD_DEFINITIONS(-DPLUGIN_A3DIO)
	ENDIF(MSVC)
	SET(plugins_sources ${plugins_sources} plugins/A3DIO.cc)
ENDIF(PLUGIN_A3DIO)

IF(PLUGIN_ARPSIO)
	IF(MSVC)
		ADD_DEFINITIONS(/DPLUGIN_ARPSIO) #it looks like some VC++ versions don't support -D syntax
	ELSE(MSVC)
		ADD_DEFINITIONS(-DPLUGIN_ARPSIO)
	ENDIF(MSVC)
	SET(plugins_sources ${plugins_sources} plugins/ARPSIO.cc)
ENDIF(PLUGIN_ARPSIO)

IF(PLUGIN_BORMAIO)
	FIND_PACKAGE(LibXML++ REQUIRED)
	INCLUDE_DIRECTORIES(${LibXML++_INCLUDE_DIRS})
	SET(plugin_libs ${plugin_libs} ${LibXML++_LIBRARIES})
	IF(MSVC)
		ADD_DEFINITIONS(/DPLUGIN_BORMAIO) #it looks like some VC++ versions don't support -D syntax
	ELSE(MSVC)
		ADD_DEFINITIONS(-DPLUGIN_BORMAIO)
	ENDIF(MSVC)
	SET(plugins_sources ${plugins_sources} plugins/BormaIO.cc)
ENDIF(PLUGIN_BORMAIO)

IF(PLUGIN_COSMOXMLIO)
	FIND_PACKAGE(LibXML++ REQUIRED)
	INCLUDE_DIRECTORIES(${LibXML++_INCLUDE_DIRS})
	SET(plugin_libs ${plugin_libs} ${LibXML++_LIBRARIES})
	IF(MSVC)
		ADD_DEFINITIONS(/DPLUGIN_COSMOXMLIO) #it looks like some VC++ versions don't support -D syntax
	ELSE(MSVC)
		ADD_DEFINITIONS(-DPLUGIN_COSMOXMLIO)
	ENDIF(MSVC)
	SET(plugins_sources ${plugins_sources} plugins/CosmoXMLIO.cc)
ENDIF(PLUGIN_COSMOXMLIO)

IF(PLUGIN_GRASSIO)
	IF(MSVC)
		ADD_DEFINITIONS(/DPLUGIN_GRASSIO) #it looks like some VC++ versions don't support -D syntax
	ELSE(MSVC)
		ADD_DEFINITIONS(-DPLUGIN_GRASSIO)
	ENDIF(MSVC)
	SET(plugins_sources ${plugins_sources} plugins/GrassIO.cc)
ENDIF(PLUGIN_GRASSIO)

IF(PLUGIN_GRIBIO)
	FIND_PACKAGE(GribApi REQUIRED)
	INCLUDE_DIRECTORIES(${GribApi_INCLUDE_DIRS})
	SET(plugin_libs ${plugin_libs} ${GribApi_LIBRARIES})
	IF(MSVC)
		ADD_DEFINITIONS(/DPLUGIN_GRIBIO) #it looks like some VC++ versions don't support -D syntax
	ELSE(MSVC)
		ADD_DEFINITIONS(-DPLUGIN_GRIBIO)
	ENDIF(MSVC)
	SET(plugins_sources ${plugins_sources} plugins/GRIBIO.cc)
ENDIF(PLUGIN_GRIBIO)

IF(PLUGIN_GEOTOPIO)
	IF(MSVC)
		ADD_DEFINITIONS(/DPLUGIN_GEOTOPIO) #it looks like some VC++ versions don't support -D syntax
	ELSE(MSVC)
		ADD_DEFINITIONS(-DPLUGIN_GEOTOPIO)
	ENDIF(MSVC)
	SET(plugins_sources ${plugins_sources} plugins/GeotopIO.cc)
ENDIF(PLUGIN_GEOTOPIO)

IF(PLUGIN_GSNIO)
	#HACK
	#gsoap generates files that do lots of warnings, there is nothing we can do about it -> off
	IF(NOT MSVC)
		SET_SOURCE_FILES_PROPERTIES(plugins/gsn/soapC.cpp PROPERTIES COMPILE_FLAGS "-w")
		SET_SOURCE_FILES_PROPERTIES(plugins/gsn/stdsoap2.cpp PROPERTIES COMPILE_FLAGS "-w")
	ENDIF(NOT MSVC)
	IF(MSVC)
		SET_SOURCE_FILES_PROPERTIES(plugins/gsn/soapC.cpp PROPERTIES COMPILE_FLAGS "/W0")
		SET_SOURCE_FILES_PROPERTIES(plugins/gsn/stdsoap2.cpp PROPERTIES COMPILE_FLAGS "/W0")
	ENDIF(MSVC)

	IF(WIN32)
		IF(MINGW)
			ADD_DEFINITIONS(-D_isnan=isnan)
		ENDIF(MINGW)
		SET(EXTRA_LINK_LIBS "ws2_32")
	ENDIF(WIN32)

	IF(MSVC)
		ADD_DEFINITIONS(/DPLUGIN_GSNIO) #it looks like some VC++ versions don't support -D syntax
	ELSE(MSVC)
		ADD_DEFINITIONS(-DPLUGIN_GSNIO)
	ENDIF(MSVC)

	SET(plugins_sources ${plugins_sources} plugins/GSNIO.cc plugins/gsn/soapGSNWebServiceSoap12BindingProxy.cpp plugins/gsn/soapC.cpp plugins/gsn/stdsoap2.cpp)
ENDIF(PLUGIN_GSNIO)

IF(PLUGIN_IMISIO)
	FIND_PACKAGE(OCCI REQUIRED)
	INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${OCCI_INCLUDE_DIR})
	SET(plugin_libs ${plugin_libs} ${OCCI_LIBRARIES})
	IF(MSVC)
		ADD_DEFINITIONS(/DPLUGIN_IMISIO) #it looks like some VC++ versions don't support -D syntax
	ELSE(MSVC)
		ADD_DEFINITIONS(-DPLUGIN_IMISIO)
	ENDIF(MSVC)
	SET(plugins_sources ${plugins_sources} plugins/ImisIO.cc)
ENDIF(PLUGIN_IMISIO)

IF(PLUGIN_PGMIO)
	IF(MSVC)
		ADD_DEFINITIONS(/DPLUGIN_PGMIO) #it looks like some VC++ versions don't support -D syntax
	ELSE(MSVC)
		ADD_DEFINITIONS(-DPLUGIN_PGMIO)
	ENDIF(MSVC)
	SET(plugins_sources ${plugins_sources} plugins/PGMIO.cc)
ENDIF(PLUGIN_PGMIO)

IF(PLUGIN_PNGIO)
	FIND_PACKAGE(PNG REQUIRED)
	INCLUDE_DIRECTORIES(${PNG_INCLUDE_DIR})
	SET(plugin_libs ${plugin_libs} ${PNG_LIBRARY})
	IF(MSVC)
		ADD_DEFINITIONS(/DPLUGIN_PNGIO) #it looks like some VC++ versions don't support -D syntax
	ELSE(MSVC)
		ADD_DEFINITIONS(-DPLUGIN_PNGIO)
	ENDIF(MSVC)
	SET(plugins_sources ${plugins_sources} plugins/PNGIO.cc)
ENDIF(PLUGIN_PNGIO)

IF(PLUGIN_SNIO)
	IF(MSVC)
		ADD_DEFINITIONS(/DPLUGIN_SNIO) #it looks like some VC++ versions don't support -D syntax
	ELSE(MSVC)
		ADD_DEFINITIONS(-DPLUGIN_SNIO)
	ENDIF(MSVC)
	SET(plugins_sources ${plugins_sources} plugins/SNIO.cc)
ENDIF(PLUGIN_SNIO)

# IF(PLUGIN_GSNIO)
# 	INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} "/usr/include")
# 	#gsoap generates files that do lots of warnings, there is nothing we can do about it -> off
# 	IF(NOT MSVC)
# 		SET_SOURCE_FILES_PROPERTIES(gsn/soapC.cpp PROPERTIES COMPILE_FLAGS "-w")
# 		SET_SOURCE_FILES_PROPERTIES(gsn/stdsoap2.cpp PROPERTIES COMPILE_FLAGS "-w")
# 	ENDIF(NOT MSVC)
# 	IF(MSVC)
# 		SET_SOURCE_FILES_PROPERTIES(gsn/soapC.cpp PROPERTIES COMPILE_FLAGS "/W0")
# 		SET_SOURCE_FILES_PROPERTIES(gsn/stdsoap2.cpp PROPERTIES COMPILE_FLAGS "/W0")
# 	ENDIF(MSVC)
# 	ADD_LIBRARY(${GSNIONAME} GSNIO.cc gsn/soapGSNWebServiceSoap12BindingProxy.cpp gsn/soapC.cpp gsn/stdsoap2.cpp)
# 	IF(WIN32)
# 		IF(MINGW)
# 			ADD_DEFINITIONS(-D_isnan=isnan)
# 		ENDIF(MINGW)
# 		SET(EXTRA_LINK_LIBS "ws2_32")
# 	ENDIF(WIN32)
# 	TARGET_LINK_LIBRARIES(${GSNIONAME} ${SHAREDNAME} ${EXTRA_LINK_LIBS})
# 	SET_TARGET_PROPERTIES(${GSNIONAME} PROPERTIES
# 		PREFIX "${LIBPREFIX}"
# 		LIBRARY_OUTPUT_DIRECTORY ${PLUGIN_DIR}
# 		RUNTIME_OUTPUT_DIRECTORY ${PLUGIN_DIR}
# 		VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
# 		SOVERSION ${VERSION_MAJOR})
# 	#necessary because gsoap uses long long
# 	IF(NOT MSVC)
# 		SET_TARGET_PROPERTIES(${GSNIONAME} PROPERTIES COMPILE_FLAGS "-Wno-long-long")
# 	ENDIF(NOT MSVC)
# 	SET(MYPLUGINS ${MYPLUGINS} ${GSNIONAME})
# ENDIF(PLUGIN_GSNIO)
