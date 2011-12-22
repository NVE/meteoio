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
#ifndef __METEODATA_H__
#define __METEODATA_H__

#include <meteoio/Date.h>
#include <meteoio/StationData.h>
#include <meteoio/IOUtils.h>

#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>

namespace mio {

/**
 * @class MeteoGrids
 * @brief A class to represent the meteorological parameters that could be contained in a grid.
 * This should be very close to MeteoData with a few additions (like the wind u,v,w)
 * @ingroup data_str
 * @author Mathias Bavay
 * @date   2011-12-22
 */

class MeteoGrids {
	public:
		/// \anchor meteogrids this enum provides names for possible meteogrids (from an ARPS file, etc)
		enum Parameters {firstparam=0,
		                 TA=firstparam, ///< Air temperature
		                 RH, ///< Relative humidity
		                 VW, ///< Wind velocity
		                 DW, ///< Wind direction
		                 VW_MAX, ///< Maximum wind velocity
		                 ISWR, ///< Incoming short wave radiation
		                 RSWR, ///< Reflected short wave radiation
		                 ILWR, ///< Incoming long wave radiation
		                 HS, ///< Snow height
		                 HNW, ///< New water equivalent height
		                 TSG, ///< Temperature ground surface
		                 TSS, ///< Temperature snow surface
		                 P, ///< Air pressure
		                 U, ///< East component of wind
		                 V, ///< North component of wind
		                 W, ///< Vertical component of wind
		                 lastparam=W};

		static const size_t nrOfParameters; ///<holds the number of meteo parameters stored in MeteoData
		static const std::string& getParameterName(const size_t& parindex);

	private:
		//static methods
		static std::vector<std::string> paramname;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map meteoparamname
};

/**
 * @class MeteoData
 * @brief A class to represent a singular measurement received from one station at a certain time (represented by the Date object)
 *
 * @ingroup data_str
 * @author Thomas Egger
 * @date   2008-12-05
 */

#ifdef _POPC_
class MeteoData : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class MeteoData {
#endif
	public:
		/// \anchor meteoparam this enum provides indexed access to meteorological fields
		enum Parameters {firstparam=0,
		                 TA=firstparam, ///< Air temperature
		                 RH, ///< Relative humidity
		                 VW, ///< Wind velocity
		                 DW, ///< Wind direction
		                 VW_MAX, ///< Maximum wind velocity
		                 ISWR, ///< Incoming short wave radiation
		                 RSWR, ///< Reflected short wave radiation
		                 ILWR, ///< Incoming long wave radiation
		                 HS, ///< Snow height
		                 HNW, ///< New water equivalent height
		                 TSG, ///< Temperature ground surface
		                 TSS, ///< Temperature snow surface
		                 P, ///< Air pressure
		                 lastparam=P};

		static const size_t nrOfParameters; ///<holds the number of meteo parameters stored in MeteoData
		static const std::string& getParameterName(const size_t& parindex);

		/**
		 * @brief The default constructor initializing every double attribute to nodata and the Date to julian==0.0
		 */
		MeteoData(void);

		/**
		* @brief A constructor that sets the measurment time
		* @param in_date A Date object representing the time of the measurement
		*/
		MeteoData(const Date& in_date);

		/**
		* @brief A setter function for the measurement date
		* @param in_date A Date object representing the time of the measurement
		*/
		void setDate(const Date& in_date);

		/**
		* @brief Add another variable to the MeteoData object,
		*        a double value will be added and the nrOfParameters increased
		* @param i_paramname A parameter name, e.g. "VSWR"
		* @return A size_t denoting the index of the the parameter added
		*/
		size_t addParameter(const std::string& i_paramname);

		/**
		* @brief Check whether a certain parameter is a part of this MeteoData instance
		* @param parname A string parameter, representing a meteo parameter, e.g. "VSWR"
		* @return A boolean indicating whether the parameter is a part of the object
		*/
		bool param_exists(const std::string& parname) const;

		/**
		 * @brief Resets all the meteo parameters to IOUtils::nodata
		 *        NOTE: member vars date and resampled are not affected
		 */
		void reset();

		bool isResampled() const;
		void setResampled(const bool&);

		void standardizeNodata(const double& plugin_nodata);

		double& operator()(const size_t& parindex);
		const double& operator()(const size_t& parindex) const;
		double& operator()(const std::string& parname);
		const double& operator()(const std::string& parname) const;

		const std::string& getNameForParameter(const size_t& parindex) const;
		size_t getParameterIndex(const std::string& parname) const;
		size_t getNrOfParameters() const;

		friend std::ostream& operator<<(std::ostream& os, const MeteoData& data);

		//Comparison operators
		bool operator==(const MeteoData&) const; ///<Operator that tests for equality
		bool operator!=(const MeteoData&) const; ///<Operator that tests for inequality

		//direct access allowed
		Date date; ///<Timestamp of the measurement
		StationData meta; ///<The meta data of the measurement

	private:
		//static methods
		static std::map<size_t, std::string> static_meteoparamname; ///<Associate a name with meteo parameters in Parameters
		static std::vector<std::string> s_default_paramname;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map meteoparamname

		//private data members
		std::vector<std::string> param_name;
		std::vector<double> data;
		size_t nrOfAllParameters;
		bool resampled;              ///<set this to true if MeteoData is result of resampling
};

typedef std::vector<MeteoData> METEO_TIMESERIE;

} //end namespace

#endif
