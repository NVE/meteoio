IF(GSNIO)
  find_program(WSDL2H wsdl2h)
  find_program(SOAPCPP2 soapcpp2)
  find_library(GSOAPPP gsoap++)

  add_custom_target(wsdl2cpp 
    wsdl2h -g A3DWebService.wsdl
    COMMAND soapcpp2 -2 -C -L -x -w -i A3DWebService.h
    COMMAND cmake -E remove *Soap11*
    command echo ${GSOAPPP_LIBRARIES})

  SET_SOURCE_FILES_PROPERTIES(soapC.cpp PROPERTIES GENERATED TRUE)
  SET_SOURCE_FILES_PROPERTIES(soapA3DWebServiceSoap12BindingProxy.cpp PROPERTIES GENERATED TRUE)

  add_executable(gsnclient test.cc soapA3DWebServiceSoap12BindingProxy.cpp soapC.cpp)

  target_link_libraries(gsnclient ${GSOAPPP})
  add_dependencies(gsnclient wsdl2cpp)

ENDIF(GSNIO)

SET_DIRECTORY_PROPERTIES(PROPERTIES 
  ADDITIONAL_MAKE_CLEAN_FILES "*~;soap*;*.nsmap;gsnclient;A3DWebService.h")
