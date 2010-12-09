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
#ifndef __STATIONDATA_H__
#define __STATIONDATA_H__

#include <meteoio/Coords.h>
#include <meteoio/IOUtils.h>

#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

namespace mio {

/**
 * @class StationData
 * @brief A class to represent meteo stations with attributes like longitude, latitude, etc.
 *
 * @author Thomas Egger
 * @date   2008-11-29
 */

#ifdef _POPC_
class StationData :POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class StationData {
#endif
	public:
		//Constructors
		/**
		* @brief The default constructor initializing every double attribute to nodata and strings to  ""
		*/
		StationData(void);

		/**
		* @brief A constructor that takes three to six arguments
		* @param _position Position of the station
		* @param _id Station's id (short identification)
		* @param _name Full name of the station (default "")
		*/
		StationData(const Coords& _position, const std::string& _id="", const std::string& _name="");

		//Specific getter functions
		std::string getStationID() const;
		std::string getStationName() const;
		Coords getPosition() const;
		
		/**
		* @brief General setter function, requires three to six arguments
		* @param _position Position of the station
		* @param _id Station's id (short identification)
		* @param _name Name of the station (default "")
		*/
		void setStationData(const Coords& _position, const std::string& _id="", const std::string& _name="");

		friend std::ostream& operator<<(std::ostream& os, const StationData& station);

		//Comparison operators
		/**
		* @brief Equality %operator
		* check all parameters but the station name
		* @return true or false
		*/
		bool operator==(const StationData&) const;
		bool operator!=(const StationData&) const; ///<Operator that tests for inequality

	public:
		Coords position;
		std::string stationID; ///<ID of the Station, typically a short string
		std::string stationName; ///<Name of the Station, detailled description
		
		//for Snowpack and other (1D) applications: add slope, aspect, horizon, etc
};

typedef std::vector<StationData> STATION_DATASET;

} //end namespace

#endif
