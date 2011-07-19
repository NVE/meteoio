#OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")

SET(meteostats_sources
	meteostats/libfit1D.cc
	meteostats/libfit1DCore.cc
	meteostats/libinterpol1D.cc
	meteostats/libinterpol2D.cc
)
