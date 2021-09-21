#include <meteoio/Grid1DInterpolator.h>

#include <iostream>
#include <utility>

using namespace std;

namespace mio {

const std::string Grid1DInterpolator::section_name("GridInterpolations1D");

Grid1DInterpolator::Grid1DInterpolator(const Config& in_cfg) : algorithm_map(), cfg(in_cfg)
{
	cfg.getValue("GRID_WINDOW_SIZE", section_name, grid_window_size, IOUtils::nothrow);
	if (grid_window_size <= 1.)
		throw IOException("GRID_WINDOW_SIZE not valid, it should be a duration in seconds at least greater than 1", AT);
	grid_window_size /= 86400.; //user uses seconds, internally julian day is used

	cfg.getValue("ENABLE_GRID_RESAMPLING", section_name, enable_grid_resampling, IOUtils::nothrow);

	//create the grid resampling algorithms for each MeteoData::Parameters entry:
	for (size_t ii = 0; ii < MeteoData::nrOfParameters; ++ii) { //loop over all MeteoData member variables
		const std::string parname( MeteoData::getParameterName(ii) ); //current semantic parameter name
		const std::string algo_name( IOUtils::strToUpper(getGridAlgorithmForParameter(parname)) );

		const std::vector< std::pair<std::string, std::string> > vecArgs( cfg.getArgumentsForAlgorithm(parname, algo_name, section_name) );
		algorithm_map[parname] = GridResamplingAlgorithmsFactory::getAlgorithm(algo_name, parname, grid_window_size, vecArgs);
	}
}

Grid1DInterpolator::~Grid1DInterpolator()
{
	for (auto it = algorithm_map.begin(); it != algorithm_map.end(); ++it)
		delete it->second;
}

bool Grid1DInterpolator::resampleData(const Date& date, const MeteoGrids::Parameters& parameter,
	const std::map<Date, Grid2DObject>& available_grids, Grid2DObject& resampled_grid)
{
	const std::string parname( MeteoGrids::getParameterName(parameter) );
	const size_t idx = MeteoData().getParameterIndex(parname);
	const MeteoData::Parameters mpar( static_cast<MeteoData::Parameters>(idx) );
	const std::string mparname( MeteoData::getParameterName(mpar) );
	if (algorithm_map.find(parname) == algorithm_map.end())
		return false;
	algorithm_map[mparname]->resample(date, available_grids, resampled_grid);
	return true; //successfull resampling
}

std::string Grid1DInterpolator::getGridAlgorithmForParameter(const std::string& parname) const
{
	std::string algorithm( "linear" ); //default value
	cfg.getValue(parname + "::resample", section_name, algorithm, IOUtils::nothrow);
	return algorithm;
}

} //namespace
