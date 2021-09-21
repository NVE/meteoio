#include <meteoio/gridResampling/GridResamplingAlgorithms.h>
#include <meteoio/gridResampling/GridLinearResampling.h>
#include <meteoio/gridResampling/GridTimeseriesResampling.h>

namespace mio {

GridResamplingAlgorithm::GridResamplingAlgorithm(const std::string& i_algorithm, const std::string& i_parname, const double& dflt_window_size,
	const std::vector< std::pair<std::string, std::string> >& /*vecArgs*/)
	: algo(i_algorithm), parname(i_parname), grid_window_size(dflt_window_size)
{
	//do nothing
}

GridResamplingAlgorithm* GridResamplingAlgorithmsFactory::getAlgorithm(const std::string& i_algorithm, const std::string& parname,
	const double& grid_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string algorithm( IOUtils::strToUpper(i_algorithm) );
	if (algorithm == "LINEAR")
		return new GridLinearResampling(algorithm, parname, grid_window_size, vecArgs);
	else if (algorithm == "TIMESERIES")
		return new GridTimeseriesResampling(algorithm, parname, grid_window_size, vecArgs);
	else
		throw IOException("The grid resampling algorithm '" + algorithm + "' is not implemented", AT);

	return nullptr;
}

} //namespace

