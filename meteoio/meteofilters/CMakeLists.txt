INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")

SET(meteofilters_sources
	meteofilters/FilterMin.cc
	meteofilters/FilterMax.cc
	meteofilters/FilterMinMax.cc
	meteofilters/FilterMedianAvg.cc
	meteofilters/FilterMeanAvg.cc
	meteofilters/FilterWindAvg.cc
	meteofilters/FilterStdDev.cc
	meteofilters/FilterTukey.cc
	meteofilters/RateFilter.cc
	meteofilters/FilterBlock.cc
	meteofilters/WindowedFilter.cc
	meteofilters/ProcessingBlock.cc
	meteofilters/ProcessingStack.cc
)
