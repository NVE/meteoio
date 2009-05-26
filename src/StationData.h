#ifndef __STATIONDATA_H__
#define __STATIONDATA_H__

#include "slfutils.h"
#include <string>
#include <sstream>
#include <iomanip>

/**
 * @class StationData
 * @brief A class to represent meteo stations with attributes like longitude, latitude, etc.
 *
 * @author Thomas Egger
 * @date   2008-11-29
 */
#ifdef _PAROC_
class StationData :POPBase {
  public:
    void Serialize(paroc_buffer &buf, bool pack);
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
   * @param easting_in easting coordinate (Swiss system)
   * @param northing_in northing coordinate (Swiss system)
   * @param alt_in Altitude
   * @param name_in Name of the station (default "")
   * @param lat_in Latitude (default StationData::nodata)
   * @param long_in Longitude (default StationData::nodata)
   */
  StationData(const double& easting_in, 
	      const double& northing_in, 
	      const double& alt_in, 
	      const std::string& name_in="",
	      const double& lat_in=nodata, 
	      const double& long_in=nodata);

  //Specific getter functions
  double getLatitude() const;
  double getLongitude() const;
  double getEasting() const;
  double getNorthing() const;
  double getAltitude() const;
  std::string getStationName() const;


  /**
   * @brief General setter function, requires three to six arguments
   * @param easting_in easting coordinate (Swiss system)
   * @param northing_in northing coordinate (Swiss system)
   * @param alt_in Altitude
   * @param name_in Name of the station (default "")
   * @param lat_in Latitude (default StationData::nodata)
   * @param long_in Longitude (default StationData::nodata)
   */
  void setStationData(const double& easting_in, 
		      const double& northing_in, 
		      const double& alt_in, 
		      const std::string& name_in="",
		      const double& lat_in=nodata, 
		      const double& long_in=nodata);

  /**
   * @brief General getter function, requires six arguments
   * @param easting_out easting coordinate (Swiss system)
   * @param northing_out northing coordinate (Swiss system)
   * @param alt_out Altitude
   * @param name_out Name of the station
   * @param lat_out Latitude 
   * @param long_out Longitude
   */
  void getStationData(double& easting_out, double& northing_out, 
		      double& alt_out, std::string& name_out,
		      double& lat_out, double& long_out) const;

  const std::string toString(void) const;

  //Comparison operators
  bool operator==(const StationData&) const; ///<Operator that tests for equality
  bool operator!=(const StationData&) const; ///<Operator that tests for inequality

 public:
  static const double nodata; ///<Takes on the value of -999.0

  double longitude, latitude, altitude, eastCoordinate, northCoordinate;
  std::string stationName; ///<Name of the Station
};

#endif
