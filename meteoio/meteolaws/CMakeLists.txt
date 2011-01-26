#OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")

SET(meteolaws_sources
	meteolaws/Atmosphere.cc
	meteolaws/Suntrajectory.cc
	meteolaws/Sun.cc
)
