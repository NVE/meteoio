// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef PROCSHIFT_H
#define PROCSHIFT_H

//#include <meteoio/meteoFilters/WindowedFilter.h> //use this one for filters relying on a data window, for example std_dev
#include <meteoio/meteoFilters/ProcessingBlock.h> //use this one for all others

#include <vector>
#include <string>

namespace mio {

/**
 * @class ProcShift
 * @brief This (empty) class is to be used as a template for developing new filters
 * @details
 * Here, write a description of how the filter operates, references to papers, etc
 * @ingroup processing
 * @author Mathias Bavay
 * @date   2014-02-19
 * @code
 * ILWR::filter1	= TEMPLATE
 * @endcode
 */

class ProcShift : public ProcessingBlock { //use this one for simple filter that only look at one data point at a time, for example min_max
	public:
		ProcShift(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const Config& cfg);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		typedef enum INTERPOL_TYPE {
			cst,
			stepwise,
			linear
		} interpol_type;
		
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		void writeOffsets(const unsigned int& param, const std::vector<MeteoData>& ivec);
		void correctOffsets(const unsigned int& param, std::vector<MeteoData>& ovec);
		
		static bool isAllNodata(const std::vector< std::pair<Date, double> >& vecX, const size_t& startIdx, const size_t& endIdx);
		static double getMedianSampling(const size_t& param, const std::vector<MeteoData>& ivec);
		
		std::vector< std::pair<Date, double> > resampleVector(const std::vector<MeteoData>& ivec, const size_t& param) const;
		double getPearson(const std::vector< std::pair<Date, double> >& vecX, const std::vector< std::pair<Date, double> >& vecY, const size_t& curr_idx, const int& offset) const;
		int getOffsetFullScan(const std::vector< std::pair<Date, double> >& vecX, const std::vector< std::pair<Date, double> >& vecY, const size_t& curr_idx, const int& range_min, const int& range_max) const;
		double getOffset(const std::vector< std::pair<Date, double> >& vecX, const std::vector< std::pair<Date, double> >& vecY, const size_t& curr_idx) const;
		
		std::string ref_param; ///< The reference parameter to compare to
		std::string root_path;
		std::string offsets_file; ///< File name that contains the extracted offsets or the correction offsets
		double cst_offset; ///< Constant correction offset, to be provided by the user
		double sampling_rate; ///< dataset sampling rate to use, either automatically extracted or provided by the user
		double offset_range; ///< range of time offsets to consider, in days
		double width_d; ///< size of the data window in days
		size_t width_idx; ///< size of the data window over which to compute the correlation
		interpol_type offsets_interp; ///< type of interpolation to use to interpolate the provided offsets
		bool extract_offsets; ///< do not apply any correcxtion but extract the time-varying offsets from the data
};

} //end namespace

#endif
