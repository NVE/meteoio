/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <meteoio/meteofilters/ProcessingBlock.h>
#include <meteoio/meteofilters/FilterMin.h>
#include <meteoio/meteofilters/FilterMax.h>
#include <meteoio/meteofilters/FilterMinMax.h>
#include <meteoio/meteofilters/FilterMeanAvg.h>
#include <meteoio/meteofilters/FilterMedianAvg.h>
#include <meteoio/meteofilters/FilterWindAvg.h>
#include <meteoio/meteofilters/FilterStdDev.h>
#include <meteoio/meteofilters/FilterRate.h>
#include <meteoio/meteofilters/FilterUnheatedHNW.h>
#include <meteoio/meteofilters/FilterTukey.h>
#include <meteoio/meteofilters/FilterMAD.h>
#include <meteoio/meteofilters/ProcButterworth.h>
#include <meteoio/meteofilters/ProcUndercatch_WMO.h>
#include <meteoio/meteofilters/ProcUndercatch_Hamon.h>
#include <meteoio/meteofilters/ProcUnventilatedT.h>
#include <meteoio/meteofilters/ProcUnshade.h>
#include <meteoio/meteofilters/ProcAdd.h>
#include <meteoio/meteofilters/ProcMult.h>
#include <meteoio/meteofilters/ProcExpSmoothing.h>
#include <meteoio/meteofilters/ProcWMASmoothing.h>

namespace mio {
/**
 * @page processing Processing overview
 * The pre-processing infrastructure is described in ProcessingBlock (for its API). The goal of this page is to give an overview of the available filters and processing elements and their usage.
 *
 * @section processing_modes Modes of operation
 * It should be noted that filters often have two modes of operations: soft or hard. In soft mode, all value that is rejected is replaced by the filter parameter's value. This means that for a soft min filter set at 0.0, all values less than 0.0 will be replaced by 0.0. In hard mode, all rejected values are replaced by nodata.
 *
 * Several filter take arguments describing a processing window (for example, FilterStdDev). In such a case, two values are given:
 * - the minimum time span of the window
 * - the minimum number of points in the window
 *
 * The ProcessingBlock will walk through the data, starting at the current point and adding points to the processing window.  When both of these criterias are met,
 * the window is accepted. This means that a window defined as "6 21600" is a window that contains 6 points minimum AND spans at least 21600 seconds.
 *
 * @section processing_section Filtering section
 * The filters are specified for each parameter in the [Filters] section. This section contains
 * a list of the various meteo parameters (see MeteoData) with their associated choice of filtering algorithms and
 * optional parameters.The filters are applied serialy, in the order they are given in. An example of such section is given below:
 * @code
 * [Filters]
 * TA::filter1	= min_max
 * TA::arg1	= 230 330
 *
 * RH::filter1	= min_max
 * RH::arg1	= -0.2 1.2
 * RH::filter2	= min_max
 * RH::arg2	= soft 0.0 1.0
 *
 * HNW::filter1	= min
 * HNW::arg1	= -0.1
 * HNW::filter2	= min
 * HNW::arg2	= soft 0.
 * @endcode
 *
 * @section processing_available Available processing elements
 * New filters can easily be developed. The filters that are currently available are the following:
 * - MIN: minimum check filter, see FilterMin
 * - MAX: maximum check filter, see FilterMax
 * - MIN_MAX: range check filter, see FilterMinMax
 * - RATE: rate of change filter, see FilterRate
 * - STD_DEV: reject data outside mean +/- k*stddev, see FilterStdDev
 * - MAD: median absolute deviation, see FilterMAD
 * - TUKEY: Tukey53H spike detection, based on median, see FilterTukey
 * - UNHEATED_RAINGAUGE: detection of snow melting in a rain gauge, see FilterUnheatedHNW
 *
 * A few data transformations are also supported besides filtering:
 * - EXP_SMOOTHING: exponential smoothing of data, see ProcExpSmoothing
 * - WMA_SMOOTHING: weighted moving average smoothing of data, see ProcWMASmoothing
 * - BUTTERWORTH: low pass butterworth filter, see ProcButterworth
 * - MEDIAN_AVG: running median average over a given window, see FilterMedianAvg
 * - MEAN_AVG: running mean average over a given window, see FilterMeanAvg
 * - WIND_AVG: vector average over a given window, see FilterWindAvg (currently, getting both vw AND dw is broken)
 * - ADD: adds a given offset to the data, see ProcAdd
 * - MULT: multiply the data by a given factor, see ProcMult
 * - UNDERCATCH_WMO: WMO rain gauge correction for undercatch, using various correction models, see ProcUndercatch_WMO
 * - UNDERCATCH_HAMON: Hamon1973 rain gauge correction for undercatch, see ProcUndercatch_Hamon
 * - UNVENTILATED_T: unventilated temperature sensor correction, see ProcUnventilatedT
 *
 */

ProcessingBlock* BlockFactory::getBlock(const std::string& blockname, const std::vector<std::string>& vec_args)
{
	if (blockname == "MIN"){
		return new FilterMin(vec_args, blockname);
	} else if (blockname == "MAX"){
		return new FilterMax(vec_args, blockname);
	} else if (blockname == "MIN_MAX"){
		return new FilterMinMax(vec_args, blockname);
	} else if (blockname == "MEAN_AVG"){
		return new FilterMeanAvg(vec_args, blockname);
	} else if (blockname == "MEDIAN_AVG"){
		return new FilterMedianAvg(vec_args, blockname);
	} else if (blockname == "WIND_AVG"){
		return new FilterWindAvg(vec_args, blockname);
	} else if (blockname == "STD_DEV"){
		return new FilterStdDev(vec_args, blockname);
	} else if (blockname == "RATE"){
		return new FilterRate(vec_args, blockname);
	} else if (blockname == "TUKEY"){
		return new FilterTukey(vec_args, blockname);
	} else if (blockname == "MAD"){
		return new FilterMAD(vec_args, blockname);
	} else if (blockname == "BUTTERWORTH"){
		return new ProcButterworth(vec_args, blockname);
	} else if (blockname == "UNHEATED_RAINGAUGE"){
		return new FilterUnheatedHNW(vec_args, blockname);
	} else if (blockname == "UNDERCATCH_WMO"){
		return new ProcUndercatch_WMO(vec_args, blockname);
	} else if (blockname == "UNDERCATCH_HAMON"){
		return new ProcUndercatch_Hamon(vec_args, blockname);
	} else if (blockname == "UNVENTILATED_T"){
		return new ProcUnventilatedT(vec_args, blockname);
	} else if (blockname == "UNSHADE"){
		return new ProcUnshade(vec_args, blockname);
	} else if (blockname == "MULT"){
		return new ProcMult(vec_args, blockname);
	} else if (blockname == "ADD"){
		return new ProcAdd(vec_args, blockname);
	} else if (blockname == "EXP_SMOOTHING"){
		return new ProcExpSmoothing(vec_args, blockname);
	} else if (blockname == "WMA_SMOOTHING"){
		return new ProcWMASmoothing(vec_args, blockname);
	} else {
		throw IOException("The processing block '"+blockname+"' does not exist! " , AT);
	}

}

ProcessingBlock::ProcessingBlock(const std::string& name) : properties(), block_name(name)
{}

void ProcessingBlock::convert_args(const size_t& min_nargs, const size_t& max_nargs,
                               const std::vector<std::string>& vec_args, std::vector<double>& dbl_args) const
{
	const size_t nr_of_args = vec_args.size();
	if ((nr_of_args < min_nargs) || (nr_of_args > max_nargs))
		throw InvalidArgumentException("Wrong number of arguments for filter/processing element \"" + getName() + "\"", AT);

	for (size_t ii=0; ii<nr_of_args; ii++){
		double tmp = IOUtils::nodata;
		IOUtils::convertString(tmp, vec_args[ii]);
		dbl_args.push_back(tmp);
	}
}

bool ProcessingBlock::is_soft(std::vector<std::string>& vec_args) {
	if (!vec_args.empty()){
		if (vec_args.front() == "soft"){
			vec_args.erase(vec_args.begin());
			return true;
		}
	}

	return false;
}

std::string ProcessingBlock::getName() const {
	return block_name;
}

ProcessingBlock::~ProcessingBlock() {}

const std::string ProcessingBlock::toString() const {
	std::ostringstream os;
	os << "[" << block_name << " ";
	os << properties.toString();
	os << "]";
	return os.str();
}

const ProcessingProperties& ProcessingBlock::getProperties() const {
	return properties;
}

const std::string ProcessingProperties::toString() const
{
	std::ostringstream os;
	const double h_before = time_before.getJulian()*24.;
	const double h_after = time_after.getJulian()*24.;
	const size_t p_before = points_before;
	const size_t p_after = points_after;

	os << "{";
	if(h_before>0. || h_after>0.) os << "-" << h_before << " +" << h_after << " h; ";
	if(p_before>0 || p_after>0) os << "-" << p_before << " +" << p_after << " pts; ";
	if(stage==ProcessingProperties::first)
		os << "p¹";
	if(stage==ProcessingProperties::second)
		os << "p²";
	if(stage==ProcessingProperties::both)
		os << "p½";
	os << "}";
	return os.str();
}

} //end namespace
