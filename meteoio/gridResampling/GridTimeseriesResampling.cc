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

	StationData mockup_meta;
	mockup_meta.setStationData(Coords(), "GRID_RES_STAT", "GRID_RES_STAT");

	for (int xx = 0; xx < resampled_grid.getNx(); ++xx) {
		for (int yy = 0; yy < resampled_grid.getNy(); ++yy) {
			std::vector<MeteoData> vecM;
			size_t index = 0;
			size_t counter = 0;
			for (auto it = all_grids.begin(); it != all_grids.end(); ++it) { //
				counter++;
				MeteoData md(it->first, mockup_meta);
				md(parname) = it->second(xx, yy);
				if (it->first > date) {
					MeteoData nodata_point( md );
					md.reset();
					vecM.push_back(nodata_point);
					index = counter;
				}
				vecM.push_back(md);
			}
			MeteoData resampled;
			ts_interpolator->resample(mockup_meta.getHash(), index, ResamplingAlgorithms::exact_match,
				MeteoData().getParameterIndex(parname),	vecM, resampled);
			resampled_grid(xx, yy) = resampled(parname);
		} //endfor yy
	} //endfor xx

	delete ts_interpolator;
}

} //namespace
