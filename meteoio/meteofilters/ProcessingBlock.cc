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
#include <meteoio/meteofilters/RateFilter.h>

namespace mio {
/**
 * @page processing Processing overview
 * The pre-processing infrastructure is described in ProcessingBlock (for its API). The goal of this page is to give an overview of the available filters and processing elements and their usage.
 *
 * @section processing_modes Modes of operation
 * It should be noted that filters often have two modes of operations: soft or hard. In soft mode, all value that is rejected is replaced by the filter parameter's value. This means that for a soft min filter set at 0.0, all values less than 0.0 will be replaced by 0.0. In hard mode, all rejected values are replaced by nodata.
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
 * The filters are being ported to the new filtering infrastructure. Only the filters whose key is capitalized have been
 * ported and are ready to use in the current version.
 * The filters that are currently available are the following:
 * - RATE: rate of change filter, see RateFilter
 * - MIN_MAX: range check filter, see FilterMinMax
 * - MIN: minimum check filter, see FilterMin
 * - MAX: maximum check filter, see FilterMax
 * - STD_DEV: reject data outside mean +/- k*stddev, see FilterAlgorithms::StandardDeviationFilter
 * - mad: median absolute deviation, see FilterAlgorithms::MedianAbsoluteDeviationFilter
 * - Tukey53H: Tukey53H spike detection, based on median, see FilterAlgorithms::Tukey53HFilter
 *
 * A few data transformations are also supported besides filtering:
 * - accumulate: data accumulates over a given period, see FilterAlgorithms::AccumulateProcess
 * - exp_smoothing: exponential smoothing of data, see FilterAlgorithms::ExpSmoothingProcess
 * - wma_smoothing window moving average smoothing of data, see FilterAlgorithms::WMASmoothingProcess
 * - MEDIAN_AVG: running median average over a given window, see FilterMedianAvg
 * - MEAN_AVG: running mean average over a given window, see FilterMeanAvg
 * - WIND_AVG: vector average over a given window, see FilterWindAvg (currently, getting both vw AND dw is broken)
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
		return new RateFilter(vec_args);
	} else {
		throw IOException("The processing block '"+blockname+"' has not been declared! " , AT);		
	}

	//return NULL; //unreachable code
}

ProcessingBlock::ProcessingBlock(const std::string& name) : block_name(name)
{}

std::string ProcessingBlock::getName() const {
	return block_name;
}

ProcessingBlock::~ProcessingBlock() {}

std::ostream& operator<<(std::ostream& os, const ProcessingBlock& data) {
	os << "[" << data.block_name << " ";
	//os << data.properties;
	os << "]";
	return os;
}


const ProcessingProperties& ProcessingBlock::getProperties() const {
	return properties;
}

std::ostream& operator<<(std::ostream& os, const ProcessingProperties& data) {
	os << "{-" << data.time_before.getJulianDate()*3600. << " +";
	os << data.time_after.getJulianDate()*3600. << " h ; ";
	os << "-" << data.points_before << " +" << data.points_after << " pts";
	if(data.for_second_pass==true)
		os << " ; 2nd pass";
	os << "}";
	return os;
}

} //end namespace
