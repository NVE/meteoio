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
#include "StationData.h"

using namespace std;

//Default constructor initializing every double attribute to nodata and strings to  ""
StationData::StationData()
{
	setStationData(nodata, nodata, nodata, "", nodata, nodata);
}

StationData::StationData(const double& x_in, const double& y_in, 
			 const double& alt_in, const std::string& name_in,
			 const double& lat_in, const double& long_in)
{
	setStationData(x_in, y_in, alt_in, name_in, lat_in, long_in);
}

void StationData::setStationData(const double& easting_in, const double& northing_in, 
				 const double& alt_in, const std::string& name_in,
				 const double& lat_in, const double& long_in)
{
	altitude = alt_in;
	stationName = name_in;
	eastCoordinate = easting_in;
	northCoordinate = northing_in;
	longitude = long_in;
	latitude = lat_in;
}

void StationData::getStationData(double& easting_out, double& northing_out, 
				 double& alt_out, std::string& name_out,
				 double& lat_out, double& long_out) const
{
	easting_out = eastCoordinate;
	northing_out = northCoordinate;
	alt_out = altitude;
	name_out = stationName;
	lat_out = latitude;
	long_out = longitude;
}


//Comparison operator
bool StationData::operator==(const StationData& in) const
{
	//latitude, longitude, eastCoordinate and northCoordinate are checked for equality in an epsilon environment
	
	const double earth_radius = 6371e3;				//in meters
	const double grid_epsilon = 5.;				//in meters
	const double long_epsilon = grid_epsilon / earth_radius;	//in degrees. small angle, so sin(x)=x
	const double lat_epsilon = long_epsilon/2.;			//in degrees. Since long is for 360deg and lat only 180, then epsilonj is 1/2

	return (IOUtils::checkEpsilonEquality(longitude, in.longitude, long_epsilon)
		&& IOUtils::checkEpsilonEquality(latitude, in.latitude, lat_epsilon) 
		&& IOUtils::checkEpsilonEquality(eastCoordinate, in.eastCoordinate, grid_epsilon) 
		&& IOUtils::checkEpsilonEquality(northCoordinate, in.northCoordinate, grid_epsilon)
		&& (altitude == in.altitude));
}

bool StationData::operator!=(const StationData& in) const{
  return !(*this==in);
}

//Specific Getter Functions for StationData
double StationData::getLatitude() const{
	return latitude;
}

double StationData::getLongitude() const
{
	return longitude;
}

double StationData::getEasting() const
{ 
	return eastCoordinate;
}

double StationData::getNorthing() const
{ 
	return northCoordinate;
}

double StationData::getAltitude() const
{
	return altitude;
}

string StationData::getStationName() const
{ 
	return stationName;
}


const string StationData::toString() const
{
	stringstream tmpstr;

	tmpstr << setprecision(10)
	 	<< "Longitude: " << setw(15) << longitude << setw(10) << "  Latitude: " << setw(15) << latitude << "  Altitude: " << altitude << endl
	 	<< "Easting:   " << setw(15) << eastCoordinate << setw(10) << "  Northing: " << setw(15) << northCoordinate << "  Name:     " << 			stationName;

	return tmpstr.str();
}
#ifdef _POPC_
void StationData::Serialize(POPBuffer &buf, bool pack)
{
	if (pack){
		buf.Pack(&altitude,1);
		buf.Pack(&stationName, 1);
		buf.Pack(&eastCoordinate,1);
		buf.Pack(&northCoordinate,1);
		buf.Pack(&longitude,1);
		buf.Pack(&latitude,1);
	}else{
		buf.UnPack(&altitude,1);
		buf.UnPack(&stationName, 1);
		buf.UnPack(&eastCoordinate,1);
		buf.UnPack(&northCoordinate,1);
		buf.UnPack(&longitude,1);
		buf.UnPack(&latitude,1);
	}
}
#endif

