INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/src" "${PROJECT_SOURCE_DIR}/filter")

SET(filter_sources 
  FilterBase1Stn.cc  
  FilterBase.cc  
  FilterFacade.cc  
  FilterValue.cc  
  MaxChangeRate.cc  
  MaxValue.cc  
  MinMaxValue.cc  
  MinValue.cc  
  NoObservedChange.cc
)

ADD_LIBRARY(filter STATIC ${filter_sources})
SET_TARGET_PROPERTIES(filter PROPERTIES 
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib
  )

#install targets:
INSTALL(TARGETS filter
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

FILE(GLOB files "${PROJECT_SOURCE_DIR}/src/filter/*.h")
INSTALL(FILES ${files} DESTINATION include/filter)
