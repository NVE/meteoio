#SPDX-License-Identifier: LGPL-3.0-or-later
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")
INCLUDE_DIRECTORIES(SYSTEM "${PROJECT_SOURCE_DIR}/meteoio/thirdParty/Eigen/")
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/meteoio/thirdParty/ctsa/src)

SET(thirdParty_sources
	thirdParty/tinyexpr.c
)

add_subdirectory(${PROJECT_SOURCE_DIR}/meteoio/thirdParty/ctsa/src)
install(DIRECTORY "${PROJECT_SOURCE_DIR}/meteoio/thirdParty/ctsa/src" 
		DESTINATION include/meteoio/thirdParty/ctsa 
		FILES_MATCHING PATTERN "*.h")
SET(thirdParty_libs ctsalib)

# Add the Eigen library as part of the install target
INSTALL(
	DIRECTORY "${PROJECT_SOURCE_DIR}/meteoio/thirdParty/Eigen/"
	DESTINATION include/meteoio/thirdParty/Eigen
        COMPONENT headers
)
