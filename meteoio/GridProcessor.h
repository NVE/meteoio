#ifndef GRIDPROCESSOR_H
#define GRIDPROCESSOR_H

#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/Config.h>
#include <meteoio/Grid1DInterpolator.h>

#include <map>
#include <string>
#include <vector>

namespace mio {

class GridProcessor {
	public:
		GridProcessor(const Config& cfg);
		bool resample(const Date& date, const MeteoGrids::Parameters& parameter, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid);
		static std::map<Date, Grid2DObject>::const_iterator seek(const Date& date, const std::map<Date, Grid2DObject>& grids, const bool& exact_match = false);
		static std::map<Date, Grid2DObject>::const_iterator seek_before(const Date& date, const std::map<Date, Grid2DObject>& grids);
		static std::map<Date, Grid2DObject>::const_iterator seek_after(const Date& date, const std::map<Date, Grid2DObject>& grids);
		double getWindowSize() const { return gi1d.getWindowSize(); };

 	private:
		static std::set<std::string> getParameters(const Config& cfg);

		Grid1DInterpolator gi1d;
		bool enable_grid_filtering = false;
};

} //end namespace

#endif
