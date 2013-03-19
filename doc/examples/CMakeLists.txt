CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../../tools/cmake/")

PROJECT(meteoio_examples)

FIND_PACKAGE(MeteoIO)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${METEOIO_INCLUDE_DIR})

#####################################################################
# add all the executables
ADD_EXECUTABLE(2D_interpolations 2D_interpolations.cc)
TARGET_LINK_LIBRARIES(2D_interpolations ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(coordinates coordinates.cc)
TARGET_LINK_LIBRARIES(coordinates ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(data_converter data_converter.cc)
TARGET_LINK_LIBRARIES(data_converter ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(dem_reading dem_reading.cc)
TARGET_LINK_LIBRARIES(dem_reading ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(grid2d_reading grid2d_reading.cc)
TARGET_LINK_LIBRARIES(grid2d_reading ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(matrix matrix.cc)
TARGET_LINK_LIBRARIES(matrix ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(meteo_reading meteo_reading.cc)
TARGET_LINK_LIBRARIES(meteo_reading ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(sun sun.cc)
TARGET_LINK_LIBRARIES(sun ${METEOIO_LIBRARIES})

ADD_EXECUTABLE(time time.cc)
TARGET_LINK_LIBRARIES(time ${METEOIO_LIBRARIES})

#####################################################################
# add a distclean target
ADD_CUSTOM_TARGET(distclean make clean
	COMMAND cmake -E remove *~ *.o
	COMMAND cmake -E remove CMakeCache.txt cmake_install.cmake
	COMMAND cmake -E remove_directory CMakeFiles
	)