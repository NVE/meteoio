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
#ifndef __BUFFER_H__
#define __BUFFER_H__

#include <meteoio/dataClasses/Date.h>
#include <meteoio/dataClasses/MeteoData.h>

namespace mio {

//we start easy: this is NOT a true ring buffer, just the same as what we previously had (for now)
class MeteoBuffer {
	public:
		MeteoBuffer() : ts_buffer(), ts_start(), ts_end() {};

		bool get(const Date& date, METEO_SET &vecMeteo) const;
		bool get(const Date& date_start, const Date& date_end, std::vector< METEO_SET > &vecMeteo) const;
		double getAvgSamplingRate() const;
		Date getBufferStart() const;
		Date getBufferEnd() const;

		bool empty() const;
		void clear();
		void push(const Date& date_start, const Date& date_end, const std::vector< METEO_SET >& vecMeteo);
		//void push(const Date& date, const METEO_SET& vecMeteo);

		const std::string toString() const;

		//HACK: these should be removed in order to hide the internals!
		//but this requires a re-write of MeteoProcessor
		std::vector< METEO_SET >& getBuffer();
		void setBufferStart(const Date& date);
		void setBufferEnd(const Date& date);
	private:
		std::vector< METEO_SET > ts_buffer; ///< stores raw data
		Date ts_start, ts_end; ///< store the beginning and the end date of the ts_buffer
};

}
#endif