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

#include "Date_IO.h"
#include "IOUtils.h"
#include <string>
#include <sstream>
#include <iomanip>

using namespace IOUtils;
/**
 * @class MeteoData
 * @brief A class to represent a singular measurement received from one station at a certain time (represented by the Date_IO object)
 *
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
		///this enum provides indexed access to meteorological fields
		enum Parameters {firstparam=0, 
					  TA=firstparam, ISWR, VW, DW, RH, LWR, HNW, TSG, TSS, HS, RSWR, P, 
					  lastparam=P};

		static const unsigned int nrOfParameters; ///<holds the number of meteo parameters stored in MeteoData

		/**
		 * @brief The default constructor initializing every double attribute to nodata and the Date_IO to julian==0.0
		 */
		MeteoData(void);

		/**
		 * @brief The copy constructor is required because of the pointers that are stored int mapparam
		 */
		MeteoData(const MeteoData&);

		MeteoData& operator=(const MeteoData&);

		/**
		* @brief A constructor that takes one to eight arguments
		* @param date_in A Date_IO object representing the time of the measurement
		* @param ta Air temperature in Kelvin (default nodata)
		* @param iswr Incoming shortwave radiation in W m-2 (default nodata)
		* @param vw Wind velocity in m s-1 (default nodata)
		* @param dw Wind direction in degrees (default nodata)
		* @param rh Relative humidity between 0 and 1 (default nodata)
		* @param lwr Long wave radiation in W m-2 (default nodata)
		* @param hnw Precipitations in mm h-1 (default nodata)
		* @param tsg Soil or snow bottom temperature in Kelvin (default nodata)
		* @param tss Soil or snow surface temperature in Kelvin (default nodata)
		* @param hs Snow height in cm (default nodata)
		* @param rswr Reflected Short Wave Radiation in W m-2 (default nodata)
		* @param p Atmospheric pressure in Pa (default nodata)
		*/
		MeteoData(const Date_IO& date_in, 
		    const double& ta=nodata, 
		    const double& iswr=nodata, 
		    const double& vw=nodata, 
		    const double& dw=nodata, 
		    const double& rh=nodata, 
		    const double& lwr=nodata, 
		    const double& hnw=nodata,
		    const double& tsg=nodata,
		    const double& tss=nodata, 
		    const double& hs=nodata,
		    const double& rswr=nodata,
		    const double& p=nodata);

		/**
		* @brief General setter function, requires one to eight arguments
		* @param date_in A Date_IO object representing the time of the measurement
		* @param ta Air temperature in Kelvin (default nodata)
		* @param iswr Incoming shortwave radiation in W m-2 (default nodata)
		* @param vw Wind velocity in m s-1 (default nodata)
		* @param dw Wind direction in degrees (default nodata)
		* @param rh Relative humidity between 0 and 1 (default nodata)
		* @param lwr Long wave radiation in W m-2 (default nodata)
		* @param hnw Precipitations in mm h-1 (default nodata)
		* @param tsg Soil or snow bottom temperature in Kelvin (default nodata)
		* @param tss Soil or snow surface temperature in Kelvin (default nodata)
		* @param hs Snow height in cm (default nodata)
		* @param rswr Reflected Short Wave Radiation in W m-2 (default nodata)
		* @param p Atmospheric pressure in Pa (default nodata)
		*/
		void setMeteoData(const Date_IO& date_in, 
		    	const double& ta=nodata, 
		    	const double& iswr=nodata, 
		    	const double& vw=nodata, 
		    	const double& dw=nodata, 
		    	const double& rh=nodata, 
		    	const double& lwr=nodata, 
		    	const double& hnw=nodata,
		    	const double& tsg=nodata,
		    	const double& tss=nodata, 
			const double& hs=nodata,
			const double& rswr=nodata,
		     const double& p=nodata);

		bool isResampled();
		void setResampled(const bool&);

		double& param(const unsigned int& parindex);
		const double& param(const unsigned int& parindex) const;
		static const std::string& getParameterName(const unsigned int& parindex);

		const std::string toString(void) const;

		//Comparison operators
		bool operator==(const MeteoData&) const; ///<Operator that tests for equality
		bool operator!=(const MeteoData&) const; ///<Operator that tests for inequality

		//direct access allowed
		Date_IO date; ///<Timestamp of the measurement
		double ta; ///<Air temperature in Kelvin
		double vw; ///<Wind velocity in m s-1
		double dw; ///<Wind direction in degrees
		double rh; ///<Relative humidity between 0 and 1
		double hnw; ///<Precipitations in mm h-1
		double iswr; ///<Incoming shortwave radiation in W m-2
		double rswr; ///<Reflected Short Wave Radiation in W m-2
		double lwr; ///<Long wave radiation in W m-2
		double tsg; ///<Soil or snow bottom temperature in Kelvin
		double tss; ///<Soil or snow surface temperature in Kelvin
		double hs; ///<Snow height in cm
		double p;  ///<Atmospheric pressure in Pa


 private:
		std::map<unsigned int, double*> meteoparam; ///<Associate an unsigned int with every meteo parameter
		static std::map<unsigned int, std::string> meteoparamname; ///<Associate a name with every meteo parameter
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map meteoparamname 
		bool resampled;              ///<set this to true if MeteoData is result of resampling
		void initParameterMap();     ///<initializes the meteoparam map that allows sequential access to meteo parameters
};

typedef std::vector<MeteoData> METEO_DATASET;

#endif
