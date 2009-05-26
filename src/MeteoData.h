#ifndef __METEODATA_H__
#define __METEODATA_H__

#include "Date_IO.h"
#include "Laws.h"
#include <string>
#include <sstream>
#include <iomanip>

/**
 * @class MeteoData
 * @brief A class to represent a singular measurement received from one station at a certain time (represented by the Date_IO object)
 *
 * @author Thomas Egger
 * @date   2008-12-05
 */
#ifdef _PAROC_
class MeteoData : POPBase {
  public:
    void Serialize(paroc_buffer &buf, bool pack);
#else
  class MeteoData {
#endif  
 public:
  /**
   * @brief The default constructor initializing every double attribute to nodata and the Date_IO to julian==0.0
   */
  MeteoData(void);

  /**
   * @brief A constructor that takes one to eight arguments
   * @param date_in A Date_IO object representing the time of the measurement
   * @param ta AIR TEMPERATURE in CELSIUS (default MeteoData::nodata)
   * @param iswr Incoming SHORTWAVE radiation in W m-2 (default MeteoData::nodata)
   * @param vw Wind VELOCITY in m s-1 (default MeteoData::nodata)
   * @param rh RELATIVE HUMIDITY (default MeteoData::nodata)
   * @param lwr LONG WAVE radiation in W m-2 (default MeteoData::nodata)
   * @param nswc NEW SNOW WATER EQUIVALENT in kg m-2 (default MeteoData::nodata)
   * @param ts0 Soil or snow bottom TEMPERATURE in CELSIUS (default MeteoData::nodata)
   */
  MeteoData(const Date_IO& date_in, 
	    const double& ta=nodata, 
	    const double& iswr=nodata, 
	    const double& vw=nodata, 
	    const double& rh=nodata, 
	    const double& lwr=nodata, 
	    const double& nswc=nodata,
	    const double& ts0=nodata);

  /**
   * @brief General setter function, requires one to eight arguments
   * @param date_in A Date_IO object representing the time of the measurement
   * @param ta AIR TEMERATURE in KELVIN (default MeteoData::nodata)
   * @param iswr Incoming SHORTWAVE radiation in W m-2 (default MeteoData::nodata)
   * @param vw Wind VELOCITY in m s-1 (default MeteoData::nodata)
   * @param rh RELATIVE HUMIDITY, between 0 and 1 (default MeteoData::nodata)
   * @param lwr LONG WAVE radiation in W m-2 (default MeteoData::nodata)
   * @param nswc NEW SNOW WATER EQUIVALENT in kg m-2 (default MeteoData::nodata)
   * @param ts0 Soil or snow bottom TEMPERATURE in KELVIN (default MeteoData::nodata)
   */
  void setMeteoData(const Date_IO& date_in, 
		    const double& ta=nodata, 
		    const double& iswr=nodata, 
		    const double& vw=nodata, 
		    const double& rh=nodata, 
		    const double& lwr=nodata, 
		    const double& nswc=nodata,
		    const double& ts0=nodata);

  /**
   * @brief Check data for plausibility and set fishy data to MeteoData::nodata
   */
  void cleanData();

  const std::string toString(void) const;

  void Check_min_max(double& param, const double low_hard, const double low_soft, const double high_soft, const double high_hard);

  //Comparison operators
  bool operator==(const MeteoData&) const; ///<Operator that tests for equality
  bool operator!=(const MeteoData&) const; ///<Operator that tests for inequality

  static const double nodata; ///<value == -999.0


  double ta, iswr, vw, rh, lwr, nswc, ts0; //direct access allowed
  Date_IO date;///<Date_IO/Time of the measurement

};

#endif
