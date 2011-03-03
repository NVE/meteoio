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
#include <meteoio/StationData.h>

using namespace std;

namespace mio {

//Default constructor initializing every double attribute to nodata and strings to  ""
StationData::StationData() : position("NULL", "NULL"), stationID(""), stationName(""),
                             slope(IOUtils::nodata), azi(IOUtils::nodata) {}

StationData::StationData(const Coords& _position, const std::string& _id, const std::string& _name)
{
	setStationData(_position, _id, _name);
	setStationSlope(IOUtils::nodata, IOUtils::nodata);
}

void StationData::setStationData(const Coords& _position, const std::string& _id, const std::string& _name)
{
	position    = _position;
	stationID   = _id;
	stationName = _name;
}

void StationData::setStationSlope(const double& in_slope_angle, const double& in_azimuth)
{
	if(in_slope_angle!=IOUtils::nodata) {
		slope = fmod(in_slope_angle, 360.);
	} else
		slope = IOUtils::nodata;

	if(in_azimuth!=IOUtils::nodata)
		azi = fmod(in_azimuth, 360.);
	else
		azi =  IOUtils::nodata;
}

//Comparison operator
bool StationData::operator==(const StationData& in) const {
	return ( (position == in.position) && (stationID == in.stationID) &&
	         (slope==in.slope) && (azi==in.azi) );// && (stationName == in.stationName));
}

bool StationData::operator!=(const StationData& in) const {
	return !(*this==in);
}

//Specific Getter Functions for stationName, stationID and position
Coords StationData::getPosition() const {
	return position;
}

std::string StationData::getStationID() const {
	return stationID;
}

std::string StationData::getStationName() const {
	return stationName;
}

double StationData::getSlopeAngle() const {
	return slope;
}

double StationData::getAzimuth() const {
	return azi;
}

std::ostream& operator<<(std::ostream& os, const StationData& station) {

	os << "<station>" << endl
	   << std::setprecision(10) << station.position 
	   << "ID:    " << station.getStationID() << endl 
	   << "Name:  " << station.getStationName() << endl
	   << "Slope: " << station.getSlopeAngle() << " bearing: " << station.getAzimuth() << endl
	   << "</station>" << endl;

	return os;
}

} //end namespace

#ifdef _POPC_
#include "marshal_meteoio.h"
using namespace mio; //HACK for POPC
void StationData::Serialize(POPBuffer &buf, bool pack)
{
	if (pack){
		marshal_Coords(buf, position, 0, FLAG_MARSHAL, NULL);
		buf.Pack(&stationID, 1);
		buf.Pack(&stationName, 1);
		buf.Pack(&slope, 1);
		buf.Pack(&azi, 1);
	}else{
		marshal_Coords(buf, position, 0, !FLAG_MARSHAL, NULL);
		buf.UnPack(&stationID, 1);
		buf.UnPack(&stationName, 1);
		buf.UnPack(&slope, 1);
		buf.UnPack(&azi, 1);
	}
}
#endif

//} //namespace //HACK for POPC
