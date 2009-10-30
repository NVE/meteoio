#ifndef __METEODATA_H__
#define __METEODATA_H__

#include "Date_IO.h"
#include "Laws.h"
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
		enum Parameters {firstparam=0, 
					  TA=firstparam, ISWR, VW, DW, RH, LWR, HNW, TSG, TSS, HS, RSWR, P, 
					  lastparam=P};

		static const unsigned int nrOfParameters; ///<holds the number of meteo parameters stored in MeteoData

		/**
		* @brief The default constructor initializing every double attribute to nodata and the Date_IO to julian==0.0
		*/
		MeteoData(void);

		/**
		* @brief A constructor that takes one to eight arguments
		* @param date_in A Date_IO object representing the time of the measurement
		* @param ta AIR TEMPERATURE in CELSIUS (default nodata)
		* @param iswr Incoming SHORTWAVE radiation in W m-2 (default nodata)
		* @param vw Wind VELOCITY in m s-1 (default nodata)
		* @param dw Wind DIRECTION in m s-1 (default nodata)
		* @param rh RELATIVE HUMIDITY (default nodata)
		* @param lwr LONG WAVE radiation in W m-2 (default nodata)
		* @param hnw Height new precipitations in kg m-2 (default nodata)
		* @param tsg Soil or snow bottom TEMPERATURE in CELSIUS (default nodata)
		* @param tss Soil or snow surface TEMPERATURE in CELSIUS (default nodata)
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
		* @param ta AIR TEMPERATURE in CELSIUS (default nodata)
		* @param iswr Incoming SHORTWAVE radiation in W m-2 (default nodata)
		* @param vw Wind VELOCITY in m s-1 (default nodata)
		* @param dw Wind DIRECTION in m s-1 (default nodata)
		* @param rh RELATIVE HUMIDITY (default nodata)
		* @param lwr LONG WAVE radiation in W m-2 (default nodata)
		* @param hnw Height new precipitations in kg m-2 (default nodata)
		* @param tsg Soil or snow bottom TEMPERATURE in CELSIUS (default nodata)
		* @param tss Soil or snow surface TEMPERATURE in CELSIUS (default nodata)
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

		/**
		* @brief Check data for plausibility and set fishy data to MeteoData::nodata
		*/
		void cleanData();

		bool isResampled();
		void setResampled(const bool&);

		const std::string toString(void) const;

		void Check_min_max(double& param, const double low_hard, const double low_soft, const double high_soft, const double high_hard);

		//Comparison operators
		bool operator==(const MeteoData&) const; ///<Operator that tests for equality
		bool operator!=(const MeteoData&) const; ///<Operator that tests for inequality

		double ta, iswr, vw, dw, rh, lwr, hnw, tsg, tss, hs, rswr, p; //direct access allowed
		Date_IO date; ///<Date_IO/Time of the measurement

		std::map<unsigned int, double*> meteoparam; ///<Associate an unsigned int with every meteo parameter
		static std::map<unsigned int, std::string> meteoparamname; ///<Associate a name with every meteo parameter
 private:
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map meteoparamname 
		bool resampled;              ///<set this to true if MeteoData is result of resampling
		void initParameterMap();     ///<initializes the meteoparam map that allows sequential access to meteo parameters
};

typedef std::vector<MeteoData> METEO_DATASET;

#endif
