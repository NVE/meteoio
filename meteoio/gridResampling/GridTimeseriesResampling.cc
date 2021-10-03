#include <meteoio/gridResampling/GridTimeseriesResampling.h>
#include <meteoio/meteoResampling/ResamplingAlgorithms.h>
#include <meteoio/GridProcessor.h>
#include <meteoio/IOUtils.h>

#include <sstream>

namespace mio {

GridTimeseriesResampling::GridTimeseriesResampling(const std::string& i_algoname, const std::string& i_parname,
	const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs)
	: GridResamplingAlgorithm(i_algoname, i_parname, dflt_window_size, vecArgs), vecArgs_(vecArgs), base_algorithm_("LINEAR")
{
	for (size_t ii = 0; ii < vecArgs.size(); ++ii) {
		if (vecArgs[ii].first == "ALGORITHM")
			base_algorithm_ = vecArgs[ii].second;
	}
}

std::string GridTimeseriesResampling::toString() const
{
	std::ostringstream ss;
	ss << std::right << std::setw(10) << parname << "::"  << std::left << std::setw(15) << algo << "[ ]";
	return ss.str();
}

void GridTimeseriesResampling::resample(const Date& date, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid)
{
	resampled_grid.set(all_grids.begin()->second, IOUtils::nodata);
	ResamplingAlgorithms* ts_interpolator = ResamplingAlgorithmsFactory::getAlgorithm(base_algorithm_, parname, 86400., vecArgs_);

	StationData point_meta; //fill this object with what's available of grid point metadata

	for (int xx = 0; xx < resampled_grid.getNx(); ++xx) {
		for (int yy = 0; yy < resampled_grid.getNy(); ++yy) {
			Coords point_coords;
			point_coords.setGridIndex(xx, yy, IOUtils::inodata); //TODO: get altitude (e. g. for solar resampling)
			resampled_grid.gridify(point_coords);
			point_meta.setStationData(point_coords, "_GRS_ID_", "_GRID_RES_STAT_"); //some unique dummy ID
			std::vector<MeteoData> vecM;
			size_t index = -1;
			size_t counter = 0;

			MeteoData resampled_pt; //point at which to resample
			for (auto it = all_grids.begin(); it != all_grids.end(); ++it) { //
				MeteoData md( it->first, point_meta );
				md(parname) = it->second(xx, yy);
				if (it->first > date && index == -1) { //put a nodata point at the date to be resampled
					resampled_pt = md;
					resampled_pt.setDate(date);
					resampled_pt.reset();
					vecM.push_back(resampled_pt);
					index = counter; //remember index of nodata point
				}
				vecM.push_back(md);
				counter++;
			}
			ts_interpolator->resample(point_meta.getHash(), index, ResamplingAlgorithms::exact_match,
				MeteoData().getParameterIndex(parname),	vecM, resampled_pt);
			resampled_grid(xx, yy) = resampled_pt(parname);
		} //endfor yy
	} //endfor xx

	delete ts_interpolator;
}

} //namespace
