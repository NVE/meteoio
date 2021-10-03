#include <meteoio/gridResampling/GridLinearResampling.h>
#include <meteoio/GridProcessor.h>
#include <meteoio/IOUtils.h>

#include <sstream>

namespace mio {

GridLinearResampling::GridLinearResampling(const std::string& i_algoname, const std::string& i_parname,
	const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs)
	: GridResamplingAlgorithm(i_algoname, i_parname, dflt_window_size, vecArgs)
{
	//do nothing
}

std::string GridLinearResampling::toString() const
{
	std::ostringstream ss;
	ss << std::right << std::setw(10) << parname << "::"  << std::left << std::setw(15) <<
		algo << "[ window_size=" << grid_window_size << " ]";
	return ss.str();
}

void GridLinearResampling::resample(const Date& date, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid)
{
	auto it_before( GridProcessor::seek_before(date, all_grids) );
	auto it_after( GridProcessor::seek_after(date, all_grids) );
	if (it_before == all_grids.end() || it_after == all_grids.end() || all_grids.size() == 1)
		throw IOException("Grids not loaded to cover linear interpolation date (is your buffer size big enough?)", AT);
	Grid2DObject grid_before = it_before->second;
	Grid2DObject grid_after = it_after->second;
	resampled_grid.set(grid_before, IOUtils::nodata);

	//solve:
	//(y - y1)/(x - x1) = (y2 - y1)/(x2 - x1)
	//==> y = y1 + (y2 - y1)/(x2 - x1) * (x - x1)
	const double x1 = it_before->first.getJulian();
	const double x2 = it_after->first.getJulian();
	const double xx = date.getJulian();
	if (x1 == x2)
		throw IOException("Equal start and end date for grid linear interpolation", AT);

	for (size_t jj = 0; jj < grid_before.size(); ++jj) {
		const double y1 = grid_before(jj);
		const double y2 = grid_after(jj);
		if ((y1 == IOUtils::nodata) || (y2 == IOUtils::nodata))
			continue; //already at nodata
		const double aa = (y2 - y1) / (x2 - x1);
		resampled_grid(jj) = y1 + aa * (xx - x1);
	} //endfor jj
}

} //namespace
