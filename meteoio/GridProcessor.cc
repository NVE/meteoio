#include <meteoio/GridProcessor.h>

#include <algorithm>

namespace mio {

GridProcessor::GridProcessor(const Config& cfg) : gi1d(cfg)
{
	cfg.getValue("ENABLE_GRID_FILTERING", "GridFilters", enable_grid_filtering, IOUtils::nothrow);
	if (enable_grid_filtering)
		std::cout << "[W] Grid filtering is not implemented yet." << std::endl;
}

bool GridProcessor::resample(const Date& date, const MeteoGrids::Parameters& parameter, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid)
{
	return gi1d.resampleData(date, parameter, all_grids, resampled_grid);
}

std::map<Date, Grid2DObject>::const_iterator GridProcessor::seek(const Date& date, const std::map<Date, Grid2DObject>& grids, const bool& exact_match)
{
	if (grids.empty())
		return grids.end();
	std::map<Date, Grid2DObject>::const_iterator it = grids.find(date);
	if (exact_match)
		return it;
	for (it = ++grids.begin(); it != grids.end(); ++it) {
		if (it->first > date)
			return --it;
	}
	return grids.end();
}

std::map<Date, Grid2DObject>::const_iterator GridProcessor::seek_before(const Date& date, const std::map<Date, Grid2DObject>& grids)
{
	if (grids.empty())
		return grids.end();
	if (date < grids.begin()->first) //there is no element before date
		return grids.end();
	for (auto it = ++grids.begin(); it != grids.end(); ++it) {
		if (it->first > date)
			return --it;
	}
	return grids.end();
}

std::map<Date, Grid2DObject>::const_iterator GridProcessor::seek_after(const Date& date, const std::map<Date, Grid2DObject>& grids)
{
	if (grids.empty())
		return grids.end();
	for (auto it = grids.begin(); it != grids.end(); ++it) {
		if (it->first > date)
			return it;
	}
	return grids.end();
}

} //namespace
