/***********************************************************************************/
/*  Copyright 2017 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef TIMEFILTERS_H
#define TIMEFILTERS_H

#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  TimeSuppr
 * @ingroup processing
 * @brief Timesteps suppression filter.
 * @details
 * This filter deletes some timesteps based on the provided arguments:
 *  - SUPPR: provide a file that contains a list of station ID's and timesteps that should be suppressed;
 *  - FRAC: suppress a given fraction of the data at random. For example, <i>0.5</i> would ensure that at least <i>50%</i> of the
 * data set's points are deleted.
 *
 * @code
 * TIME::filter1     = suppr
 * TIME::arg1::suppr = ./input/meteo/suppr.dat
 * @endcode
 * 
 * The file <i>suppr.dat</i> would look like this (the time is given in the timezone declared in Input::TIME_ZONE):
 * @code
 * *WFJ 2015-10-01T12:00
 * *DAV 2015-10-02T15:00
 * *WFJ 2015-11-10T06:00
 * *WFJ 2015-12-25T01:00 2015-12-27T13:30
 * *WFJ 2015-09-01T07:15 - 2015-09-10T20:30
 * STB2 2015-10-01T21:30
 * @endcode
 * Time ranges are declared by providing two dates on the same line. For more visibility, the said two dates can be separated by " - " (which a white
 * space character on both sides, as shown in the example above).
 */
	
class TimeSuppr : public ProcessingBlock {
	public:
		TimeSuppr(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const std::string& root_path, const double& TZ);

		void process(const unsigned int& param, const std::vector<MeteoData>& ivec, std::vector<MeteoData>& ovec);

	private:
		void supprByDates(std::vector<MeteoData>& ovec) const;
		void supprFrac(std::vector<MeteoData>& ovec) const;
		
		std::map< std::string, std::vector<dates_range> > suppr_dates;
		double range;
};

} //end namespace

#endif
