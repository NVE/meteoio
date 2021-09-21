#ifndef GRIDRESAMPLINGALGORITHM_H
#define GRIDRESAMPLINGALGORITHM_H

#include <meteoio/dataClasses/Grid2DObject.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/IOUtils.h>

#include <string>
#include <vector>

namespace mio {

class GridResamplingAlgorithm {

	public:
		GridResamplingAlgorithm(const std::string& i_algorithm, const std::string& i_parname, const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& /*vecArgs*/);
		virtual void resample(const Date& date, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid) = 0;

 	protected:
		const std::string algo, parname;
		double grid_window_size;
};

class GridResamplingAlgorithmsFactory {
	public:
		static GridResamplingAlgorithm* getAlgorithm(const std::string& i_algorithm, const std::string& parname, const double& grid_window_size,
			const std::vector< std::pair<std::string, std::string> >& vecArgs);
};

} //end namespace

#endif
