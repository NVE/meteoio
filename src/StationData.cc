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

namespace mio {

//Default constructor initializing every double attribute to nodata and strings to  ""
StationData::StationData() : position("NULL", "NULL"), stationName(""){}

StationData::StationData(const Coords& _position, const std::string& name_in)
{
	setStationData(_position, name_in);
}

void StationData::setStationData(const Coords& _position, const std::string& name_in)
{
	position = _position;
	stationName = name_in;
}

void StationData::getStationData(Coords& position_out, std::string& name_out) {
	position_out = position;
	name_out = stationName;
}


//Comparison operator
bool StationData::operator==(const StationData& in) const {
	return ( (position == in.position)
		&& (stationName == in.stationName));
}

bool StationData::operator!=(const StationData& in) const {
  return !(*this==in);
}

//Specific Getter Functions for StationData
std::string StationData::getStationName() const {
	return stationName;
}

std::ostream& operator<<(std::ostream& os, const StationData& station) {

	os << "<station>\n";
	os << std::setprecision(10) << station.position << "  Name:     " << station.stationName << "\n";
	os << "</station>\n";

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
		buf.Pack(&stationName, 1);
	}else{
		marshal_Coords(buf, position, 0, !FLAG_MARSHAL, NULL);
		buf.UnPack(&stationName, 1);
	}
}
#endif

//} //namespace //HACK for POPC
