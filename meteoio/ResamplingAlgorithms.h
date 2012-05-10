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
#ifndef __RESAMPLINGALGORITHMS_H__
#define __RESAMPLINGALGORITHMS_H__

#include <meteoio/MeteoData.h>
#include <meteoio/StationData.h>
#include <meteoio/meteostats/libinterpol1D.h>

#include <iostream>
#include <string>
#include <vector>
#include <map>

namespace mio {

/**
 * @class ResamplingAlgorithms
 * @brief Temporal resampling algorithms
 *
 * @ingroup stats
 * @author Thomas Egger
 * @date   2010-06-25
 */
class ResamplingAlgorithms {
	public:
		enum ResamplingPosition {
			exact_match,
			before,
			end
		};

		typedef void(*resamplingptr)(const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
							    const std::vector<std::string>& taskargs, const double& window_size, const std::vector<MeteoData>& vecM, MeteoData& md);

		static const resamplingptr& getAlgorithm(const std::string& algorithmname);

		//Available algorithms
		static void NoResampling(const size_t& index, const ResamplingPosition& position, const size_t& paramindex, 
							const std::vector<std::string>& taskargs, const double& window_size, const std::vector<MeteoData>& vecM, MeteoData& md);
		static void LinearResampling(const size_t& index, const ResamplingPosition& position, const size_t& paramindex, 
							    const std::vector<std::string>& taskargs, const double& window_size, const std::vector<MeteoData>& vecM, MeteoData& md);
		static void NearestNeighbour(const size_t& index, const ResamplingPosition& position, const size_t& paramindex, 
							    const std::vector<std::string>& taskargs, const double& window_size, const std::vector<MeteoData>& vecM, MeteoData& md);

		static void Accumulate(const size_t& index, const ResamplingPosition& position, const size_t& paramindex, 
						   const std::vector<std::string>& taskargs, const double& window_size, const std::vector<MeteoData>& vecM, MeteoData& md);

 	private:
		static double funcval(const size_t& position, const size_t& paramindex, const std::vector<MeteoData>& vecM,
		                      const Date& date);
		static void getNearestValidPts(const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date,
		                               const double& window_size, size_t& indexP1, size_t& indexP2);
		static double linearInterpolation(const double& x1, const double& y1,
		                                  const double& x2, const double& y2, const double& x3);

		static std::map<std::string, resamplingptr> algorithmMap;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map algorithmMap
};
} //end namespace

#endif
