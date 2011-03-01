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
		                 ISWR, ///< Incoming short wave radiation
		                 RSWR, ///< Reflected short wave radiation
		                 ILWR, ///< Incoming long wave radiation
		                 HS, ///< Snow height
		                 HNW, ///< New water equivalent height
		                 TSG, ///< Temperature ground surface
		                 TSS, ///< Temperature snow surface
		                 P, ///< Air pressure
		                 lastparam=P};

		static const unsigned int nrOfParameters; ///<holds the number of meteo parameters stored in MeteoData
		static const std::string& getParameterName(const unsigned int& parindex);

		/**
		 * @brief The default constructor initializing every double attribute to nodata and the Date to julian==0.0
		 */
		MeteoData(void);

		/**
		 * @brief The copy constructor is required because of the pointers that are stored int mapparam
		 */
		MeteoData(const MeteoData&);

		MeteoData& operator=(const MeteoData&);

		/**
		* @brief A constructor that sets the measurment time
		* @param in_date A Date object representing the time of the measurement
		*/
		MeteoData(const Date& in_date);

		/**
		* @brief General setter function, requires one to eight arguments
		* @param param One element out of the Parameters enum, e.g. MeteoData::TA
		* @param value A double value to be assigned to the respective variable denoted in the first argument
		*/
		void setData(const MeteoData::Parameters& param, const double& value);

		/**
		* @brief A setter function for the measurement date
		* @param in_date A Date object representing the time of the measurement
		*/
		void setDate(const Date& in_date);

		/**
		* @brief Add another variable to the MeteoData object, 
		*        a double value will be added and the nrOfParameters increased
		* @param i_paramname A parameter name, e.g. "VSWR"
		*/
		void addParameter(const std::string& i_paramname);

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

		double& param(const unsigned int& parindex);
		const double& param(const unsigned int& parindex) const;
		double& param(const std::string& parname);
		const double& param(const std::string& parname) const;
		const std::string& getNameForParameter(const unsigned int& parindex) const;
		unsigned int getParameterIndex(const std::string& parname) const;

		friend std::ostream& operator<<(std::ostream& os, const MeteoData& data);

		//Comparison operators
		bool operator==(const MeteoData&) const; ///<Operator that tests for equality
		bool operator!=(const MeteoData&) const; ///<Operator that tests for inequality

		//direct access allowed
		Date date; ///<Timestamp of the measurement
		StationData meta; ///<The meta data of the measurement
		double ta; ///<Air temperature in Kelvin
		double vw; ///<Wind velocity in m s-1
		double dw; ///<Wind direction in degrees
		double rh; ///<Relative humidity between 0 and 1
		double hnw; ///<Precipitations in mm h-1
		double iswr; ///<Incoming shortwave radiation in W m-2
		double rswr; ///<Reflected Short Wave Radiation in W m-2
		double ilwr; ///<Incoming Long wave radiation in W m-2
		double tsg; ///<Soil or snow bottom temperature in Kelvin
		double tss; ///<Soil or snow surface temperature in Kelvin
		double hs; ///<Snow height in m
		double p;  ///<Atmospheric pressure in Pa

		unsigned int getNrOfParameters() const;
 private:
		static std::map<unsigned int, std::string> static_meteoparamname; ///<Associate a name with meteo parameters in Parameters
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map meteoparamname 

		unsigned int nrOfAllParameters;

		std::map<std::string, double> extraparameters; ///<All non-standard meteo parameters will end up in this map
		std::map<std::string, double*> mapParameterByName; ///<Associate name and meteo parameter
		std::map<unsigned int, double*> meteoparam; ///<Associate an unsigned int with every meteo parameter
		std::map<unsigned int, std::string> meteoparamname; ///<Associate a name with every meteo parameter

		void initAllParameters();
		void associateMemberVariables();
		void initParameterMap();     ///<initializes the meteoparam map that allows sequential access to meteo parameters
		bool resampled;              ///<set this to true if MeteoData is result of resampling
};

typedef std::vector<MeteoData> METEO_TIMESERIE;

} //end namespace

#endif
