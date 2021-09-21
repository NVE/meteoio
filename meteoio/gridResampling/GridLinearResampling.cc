#include <meteoio/gridResampling/GridLinearResampling.h>
#include <meteoio/GridProcessor.h>
#include <meteoio/IOUtils.h>

#include <sstream>

namespace mio {

GridLinearResampling::GridLinearResampling(const std::string& i_algoname, const std::string& i_parname,
	const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs)
	: GridResamplingAlgorithm(i_algoname, i_parname, dflt_window_size, vecArgs)
{
	const std::string where( "GridInterpolations1D::"+i_parname+"::"+i_algoname );
	if (!vecArgs.empty()) //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments for \""+where+"\"", AT);
}

std::string GridLinearResampling::toString() const
{
	std::ostringstream ss;
	ss << std::right << std::setw(10) << parname << "::"  << std::left << std::setw(15) << algo << "[ ]";
	return ss.str();
}

void GridLinearResampling::resample(const Date& date, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid)
{
	auto it_before( GridProcessor::seek_before(date, all_grids) );
	auto it_after( GridProcessor::seek_after(date, all_grids) );
	Grid2DObject grid_before = it_before->second;
	Grid2DObject grid_after = it_after->second;
	resampled_grid.set(grid_before, IOUtils::nodata);

	const double ts = it_before->first.getJulian();
	const double te = it_after->first.getJulian();
	const double tt = date.getJulian();
	if (ts == te)
		throw IOException("Equal start and end date for grid linear interpolation", AT);

	for (int xx = 0; xx < grid_before.getNx(); ++xx) {
		for (int yy = 0; yy < grid_before.getNy(); ++yy) {
			if (!grid_before(xx, yy) == IOUtils::nodata && !grid_after(xx, yy) == IOUtils::nodata) {
				const double aa = (grid_after(xx, yy) - grid_before(xx, yy) / (te - ts));
				const double bb = grid_after(xx, yy) - aa * te;
				resampled_grid(xx, yy) = aa * tt + bb;
			}
		} //endfor yy
	} //endfor xx
}

} //namespace
