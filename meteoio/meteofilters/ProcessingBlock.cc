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
#include <meteoio/meteofilters/FilterHNWMelt.h>
#include <meteoio/meteofilters/FilterTukey.h>
#include <meteoio/meteofilters/FilterMAD.h>
#include <meteoio/meteofilters/ProcUndercatch.h>
#include <meteoio/meteofilters/ProcPassiveT.h>
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
 * - RATE: rate of change filter, see FilterRate
 * - MIN_MAX: range check filter, see FilterMinMax
 * - MIN: minimum check filter, see FilterMin
 * - MAX: maximum check filter, see FilterMax
 * - STD_DEV: reject data outside mean +/- k*stddev, see FilterStdDev
 * - MAD: median absolute deviation, see FilterMAD
 * - TUKEY: Tukey53H spike detection, based on median, see FilterTukey
 * - HNW_MELT: detection of snow melting in a rain gauge, see FilterHNWMelt
 *
 * A few data transformations are also supported besides filtering:
 * - EXP_SMOOTHING: exponential smoothing of data, see ProcExpSmoothing
 * - WMA_SMOOTHING window moving average smoothing of data, see ProcWMASmoothing
 * - MEDIAN_AVG: running median average over a given window, see FilterMedianAvg
 * - MEAN_AVG: running mean average over a given window, see FilterMeanAvg
 * - WIND_AVG: vector average over a given window, see FilterWindAvg (currently, getting both vw AND dw is broken)
 * - ADD: adds a given offset to the data, see ProcAdd
 * - MULT: multiply the data by a given factor, see ProcMult
 * - UNDERCATCH: rain gauge correction for undercatch, using various correction models, see ProcUndercatch
 * - PASSIVE_T: unventilated temperature sensor correction, see ProcPassiveT
 *
 */

std::set<std::string> BlockFactory::availableBlocks;
const bool BlockFactory::__init = BlockFactory::initStaticData();

bool BlockFactory::initStaticData()
{
	availableBlocks.insert("MIN");
	availableBlocks.insert("MAX");
	availableBlocks.insert("MIN_MAX");
	availableBlocks.insert("MEAN_AVG");
	availableBlocks.insert("MEDIAN_AVG");
	availableBlocks.insert("WIND_AVG");
	availableBlocks.insert("STD_DEV");
	availableBlocks.insert("RATE");
	availableBlocks.insert("TUKEY");
	availableBlocks.insert("MAD");
	availableBlocks.insert("HNW_MELT");
	availableBlocks.insert("UNDERCATCH");
	availableBlocks.insert("PASSIVE_T");
	availableBlocks.insert("ADD");
	availableBlocks.insert("MULT");
	availableBlocks.insert("EXP_SMOOTHING");
	availableBlocks.insert("WMA_SMOOTHING");
	return true;
}

ProcessingBlock* BlockFactory::getBlock(const std::string& blockname, const std::vector<std::string>& vec_args)
{
	//Check whether algorithm theoretically exists
	if (availableBlocks.find(blockname) == availableBlocks.end())
		throw UnknownValueException("The processing block '"+blockname+"' does not exist" , AT);

	if (blockname == "MIN"){
		return new FilterMin(vec_args);
	} else if (blockname == "MAX"){
		return new FilterMax(vec_args);
	} else if (blockname == "MIN_MAX"){
		return new FilterMinMax(vec_args);
	} else if (blockname == "MEAN_AVG"){
		return new FilterMeanAvg(vec_args);
	} else if (blockname == "MEDIAN_AVG"){
		return new FilterMedianAvg(vec_args);
	} else if (blockname == "WIND_AVG"){
		return new FilterWindAvg(vec_args);
	} else if (blockname == "STD_DEV"){
		return new FilterStdDev(vec_args);
	} else if (blockname == "RATE"){
		return new FilterRate(vec_args);
	} else if (blockname == "TUKEY"){
		return new FilterTukey(vec_args);
	} else if (blockname == "MAD"){
		return new FilterMAD(vec_args);
	} else if (blockname == "HNW_MELT"){
		return new FilterHNWMelt(vec_args);
	} else if (blockname == "UNDERCATCH"){
		return new ProcUndercatch(vec_args);
	} else if (blockname == "PASSIVE_T"){
		return new ProcPassiveT(vec_args);
	} else if (blockname == "MULT"){
		return new ProcMult(vec_args);
	} else if (blockname == "ADD"){
		return new ProcAdd(vec_args);
	} else if (blockname == "EXP_SMOOTHING"){
		return new ProcExpSmoothing(vec_args);
	} else if (blockname == "WMA_SMOOTHING"){
		return new ProcWMASmoothing(vec_args);
	} else {
		throw IOException("The processing block '"+blockname+"' has not been declared! " , AT);
	}

}

ProcessingBlock::ProcessingBlock(const std::string& name) : block_name(name)
{}

void ProcessingBlock::convert_args(const unsigned int& min_nargs, const unsigned int& max_nargs,
                               const std::vector<std::string>& vec_args, std::vector<double>& dbl_args)
{
	if ((vec_args.size() < min_nargs) || (vec_args.size() > max_nargs))
		throw InvalidArgumentException("Wrong number of arguments for filter/processing element \"" + getName() + "\"", AT);

	for (size_t ii=0; ii<vec_args.size(); ii++){
		double tmp = IOUtils::nodata;
		IOUtils::convertString(tmp, vec_args[ii]);
		dbl_args.push_back(tmp);
	}
}

std::string ProcessingBlock::getName() const {
	return block_name;
}

ProcessingBlock::~ProcessingBlock() {}

std::ostream& operator<<(std::ostream& os, const ProcessingBlock& data) {
	os << "[" << data.block_name << " ";
	os << data.properties;
	os << "]";
	return os;
}

const ProcessingProperties& ProcessingBlock::getProperties() const {
	return properties;
}

std::ostream& operator<<(std::ostream& os, const ProcessingProperties& data) {
	const double h_before = data.time_before.getJulianDate()*24.;
	const double h_after = data.time_after.getJulianDate()*24.;
	const unsigned int p_before = data.points_before;
	const unsigned int p_after = data.points_after;

	os << "{";
	if(h_before>0. || h_after>0.) os << "-" << h_before << " +" << h_after << " h; ";
	if(p_before>0 || p_after>0) os << "-" << p_before << " +" << p_after << " pts; ";
	if(data.stage==ProcessingProperties::first)
		os << "p¹";
	if(data.stage==ProcessingProperties::second)
		os << "p²";
	if(data.stage==ProcessingProperties::both)
		os << "p½";
	os << "}";
	return os;
}

} //end namespace
