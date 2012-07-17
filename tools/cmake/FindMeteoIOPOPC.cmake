INCLUDE(LibFindMacros)

# Finally the library itself
GET_FILENAME_COMPONENT(SRC_DIR ${CMAKE_SOURCE_DIR} PATH) #ie goes up one level
STRING(REPLACE " " "\\ " SRC_DIR ${SRC_DIR})

IF(WIN32)
	MESSAGE("POPC is not supported under Windows!")
ELSE(WIN32)
	IF(APPLE)
		FIND_LIBRARY(METEOIO_LIBRARY
		NAMES meteoiopopc
		PATHS
			"/Applications/MeteoIO/lib"
			ENV LD_LIBRARY_PATH
			ENV DYLD_FALLBACK_LIBRARY_PATH
			"~/usr/lib"
			"/usr/local/lib"
			"/usr/lib"
			"/opt/lib"
			${SRC_DIR}/meteoio/lib
		DOC "Location of the libmeteoio, like /usr/lib/libmeteoiopopc.dylib"
		)
	ELSE(APPLE)
		FIND_LIBRARY(METEOIO_LIBRARY
		NAMES meteoiopopc
		PATHS
			ENV LD_LIBRARY_PATH
			"~/usr/lib"
			"/usr/local/lib"
			"/usr/lib"
			"/opt/lib"
			${SRC_DIR}/meteoio/lib
		DOC "Location of the libmeteoio, like /usr/lib/libmeteoiopopc.so"
		)
	ENDIF(APPLE)
ENDIF(WIN32)

#build METEOIO_ROOT so we can provide a hint for searching for the header file
GET_FILENAME_COMPONENT(meteoio_libs_root ${METEOIO_LIBRARY} PATH)
SET(METEOIO_ROOT "${meteoio_libs_root}/../")
STRING(REPLACE  " " "\\ " METEOIO_ROOT ${METEOIO_ROOT})

# locate main header file
FIND_PATH(METEOIO_INCLUDE_DIR
  NAMES meteoio/MeteoIO.h
  #HINTS ${METEOIO_ROOT}/include
  PATHS
	"${METEOIO_ROOT}/include"
	"${METEOIO_ROOT}"
	"~/usr/include"
	"/usr/local/include"
	"/usr/include"
	"/opt/include"
  DOC "Location of the meteoio headers, like /usr/include"
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
SET(METEOIO_PROCESS_INCLUDES METEOIO_INCLUDE_DIR)
SET(METEOIO_PROCESS_LIBS METEOIO_LIBRARY)
libfind_process(METEOIO)
