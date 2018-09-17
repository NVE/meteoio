#############################
## Global Configuration    ##
#############################

# set Memory check command
set(CTEST_MEMORYCHECK_COMMAND valgrind)

#get everytyhing needed to make executables
FIND_PACKAGE(MeteoIO)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${METEOIO_INCLUDE_DIR})
SET(LIBRARIES ${METEOIO_LIBRARIES})

# go back to the source dir to have all .h files as would be installed
INCLUDE_DIRECTORIES(..)

###################
## Add Tests     ##
###################
ADD_SUBDIRECTORY(meteo_reading_interpol)
ADD_SUBDIRECTORY(meteo_reading_no_interpol)
ADD_SUBDIRECTORY(sun)
ADD_SUBDIRECTORY(dem_reading)
ADD_SUBDIRECTORY(2D_interpolations)
ADD_SUBDIRECTORY(arrays)
ADD_SUBDIRECTORY(coords)
ADD_SUBDIRECTORY(stats)
