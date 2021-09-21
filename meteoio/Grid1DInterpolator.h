#ifndef GRID1DINTERPOLATOR_H
#define GRID1DINTERPOLATOR_H

#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/dataClasses/Grid2DObject.h>
#include <meteoio/Config.h>
#include <meteoio/gridResampling/GridResamplingAlgorithms.h>

#include <string>
#include <vector>
#include <map>

namespace mio {

class Grid1DInterpolator {
	public:
		Grid1DInterpolator(const Config& in_cfg);
		~Grid1DInterpolator();
		Grid1DInterpolator& operator=(const Grid1DInterpolator&); ///<Assignement operator
		bool resampleData(const Date& date, const MeteoGrids::Parameters& parameter, const std::map<Date, Grid2DObject>& available_grids, Grid2DObject& resampled_grid);
		double getWindowSize() const { return grid_window_size; };

 	private:
		std::string getGridAlgorithmForParameter(const std::string& parname) const;

		static const std::string section_name;
		std::map<std::string, GridResamplingAlgorithm*> algorithm_map; //per parameter interpolation algorithms
		const Config& cfg;
		double grid_window_size = 86400.; ///< in seconds, default is 2 Julian days
		bool enable_grid_resampling = true; ///< easy way to turn grid resampling on/off
};

} //end namespace

#endif
