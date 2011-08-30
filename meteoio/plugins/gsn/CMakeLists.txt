IF(COMMAND cmake_policy)
	CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

IF(PLUGIN_GSNIO)
	FIND_PROGRAM(WSDL2H wsdl2h)
	FIND_PROGRAM(SOAPCPP2 soapcpp2)
	FIND_LIBRARY(GSOAPPP gsoap++)

	ADD_CUSTOM_TARGET(wsdl2cpp
		#wsdl2h -g A3DWebService.wsdl
		#COMMAND soapcpp2 -2 -C -L -x -w -i A3DWebService.h
		#COMMAND cmake -E remove *Soap11*
		command echo ${GSOAPPP_LIBRARIES}
	)

	#SET_SOURCE_FILES_PROPERTIES(soapC.cpp PROPERTIES GENERATED TRUE)
	#SET_SOURCE_FILES_PROPERTIES(soapA3DWebServiceSoap12BindingProxy.cpp PROPERTIES GENERATED TRUE)

	#ADD_EXECUTABLE(gsnclient gsnclient.cc soapA3DWebServiceSoap12BindingProxy.cpp soapC.cpp)
	#TARGET_LINK_LIBRARIES(gsnclient ${GSOAPPP} stdc++)
	#ADD_DEPENDENCIES(gsnclient wsdl2cpp)
ENDIF(PLUGIN_GSNIO)

SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "*~;soap*;*.nsmap;gsnclient;A3DWebService.h")
