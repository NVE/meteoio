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
#include <meteoio/BufferedIOHandler.h>

using namespace std;

namespace mio {

#ifdef _POPC_
BufferedIOHandler::BufferedIOHandler(IOHandler& _iohandler, const Config& _cfg) 
	: iohandler(_iohandler), cfg(_cfg), meteoprocessor(_cfg), meteoBuffer(), startDateBuffer(), endDateBuffer(), mapBufferedGrids()
#else
BufferedIOHandler::BufferedIOHandler(IOHandler& _iohandler, const Config& _cfg) 
	  : IOInterface(NULL), iohandler(_iohandler), cfg(_cfg), meteoprocessor(_cfg), meteoBuffer(), startDateBuffer(), endDateBuffer(), mapBufferedGrids()
#endif
{
	setDfltBufferProperties();
}

#ifdef _POPC_
BufferedIOHandler::~BufferedIOHandler()
#else
BufferedIOHandler::~BufferedIOHandler() throw()
#endif
{
	setDfltBufferProperties();
}

void BufferedIOHandler::read2DGrid(Grid2DObject& _grid2Dobj, const std::string& _filename)
{
	std::map<std::string, Grid2DObject>::iterator it = mapBufferedGrids.find(_filename);
	if (it != mapBufferedGrids.end()) { //already in map
		_grid2Dobj = (*it).second; 
		return;
	}
	
	Grid2DObject tmpgrid2D;
	iohandler.read2DGrid(tmpgrid2D, _filename);
	mapBufferedGrids[_filename] = tmpgrid2D;
	_grid2Dobj = tmpgrid2D;
}

void BufferedIOHandler::readDEM(DEMObject& _grid2Dobj)
{
	std::map<std::string, Grid2DObject>::iterator it = mapBufferedGrids.find("/:DEM");
	if (it != mapBufferedGrids.end()) {
		//already in map. If the update properties have changed,
		//we copy the ones given in input and force the update of the object
		const DEMObject::update_type in_ppt = (DEMObject::update_type)_grid2Dobj.getUpdatePpt();
		_grid2Dobj = (*it).second;
		const DEMObject::update_type buff_ppt = (DEMObject::update_type)_grid2Dobj.getUpdatePpt();
		if(in_ppt!=buff_ppt) {
			_grid2Dobj.setUpdatePpt(in_ppt);
			_grid2Dobj.update();
		}
		return;
	}
	
	DEMObject tmpgrid2D;
	 //copy the updating policy of the destination
	tmpgrid2D.setUpdatePpt((DEMObject::update_type)_grid2Dobj.getUpdatePpt());
	iohandler.readDEM(tmpgrid2D);
	mapBufferedGrids["/:DEM"] = tmpgrid2D;
	_grid2Dobj = tmpgrid2D;
}

void BufferedIOHandler::readLanduse(Grid2DObject& _grid2Dobj)
{
	std::map<std::string, Grid2DObject>::iterator it = mapBufferedGrids.find("/:LANDUSE");
	if (it != mapBufferedGrids.end()) { //already in map
		_grid2Dobj = (*it).second; 
		return;
	}
	
	Grid2DObject tmpgrid2D;
	iohandler.readLanduse(tmpgrid2D);
	mapBufferedGrids["/:LANDUSE"] = tmpgrid2D;
	_grid2Dobj = tmpgrid2D;
}

void BufferedIOHandler::readAssimilationData(const Date& _date, Grid2DObject& _grid2Dobj)
{
	std::map<std::string, Grid2DObject>::iterator it = mapBufferedGrids.find("/:ASSIMILATIONDATA" + _date.toString(Date::FULL));
	if (it != mapBufferedGrids.end()) { //already in map
		_grid2Dobj = (*it).second; 
		return;
	}
	
	Grid2DObject tmpgrid2D;
	iohandler.readAssimilationData(_date, tmpgrid2D);
	mapBufferedGrids["/:ASSIMILATIONDATA" + _date.toString(Date::FULL)] = tmpgrid2D;
	_grid2Dobj = tmpgrid2D;
}

void BufferedIOHandler::readStationData(const Date& date, std::vector<StationData>& vecStation)
{
	iohandler.readStationData(date, vecStation);
}

#ifdef _POPC_
void BufferedIOHandler::writeMeteoData(std::vector< std::vector<MeteoData> >& vecMeteo,
							    const std::string& name)
#else 
void BufferedIOHandler::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo,
							    const std::string& name)
#endif
{
	iohandler.writeMeteoData(vecMeteo, name);
}

void BufferedIOHandler::setDfltBufferProperties()
{
	always_rebuffer = true;
	bufferbefore = Date(2.0);  //minus 2 days
	bufferafter = Date(20.0);  //plus 20 days
}

void BufferedIOHandler::readMeteoData(const Date& i_date, std::vector<MeteoData>& vecMeteo){
	/* For every station: 
	 * 1) See whether data is already buffered
	 * 2) Filter - includes resampling
	 * 3) Return the values
	 */

	vecMeteo.clear();

	if (meteoBuffer.size() == 0){ //init
		bufferAllData(i_date);
	}

	//loop through all meteo buffers, there is one for each station
	for (unsigned int ii=0; ii<meteoBuffer.size(); ii++) {
		unsigned int index = IOUtils::npos;

		if (meteoBuffer[ii].size() > 0) {//check whether meteo data for the date exists in buffer
			index = IOUtils::seek(i_date, meteoBuffer[ii], false);
		}

		if (index == IOUtils::npos) { //not in buffer
			//Check buffering strategy
			bool rebuffer = false;
			if ((startDateBuffer.at(ii) > i_date) || (endDateBuffer.at(ii) <= i_date)){
				rebuffer = true;
			} else { 
				if (always_rebuffer) rebuffer = true;
			}

			if (rebuffer){
				//cout << "[I] Station " << ii << "(" << stationID 
				//	<< ") data for date " << i_date.toString(Date::FULL) << " not in buffer ..." << endl;
				
				bool dataexists = bufferData(i_date, ii);
				if (dataexists) {//i_date is contained in buffer
					index = IOUtils::seek(i_date, meteoBuffer[ii], false);
				}
			}
		}

		//APPLY FILTERS
		MeteoData md; 
		if (index != IOUtils::npos)
			meteoprocessor.processData(i_date, meteoBuffer[ii], md);

		//Check whether StationData is meaningful, try to get meaningful meta data
		if (index != IOUtils::npos) {
			vecMeteo.push_back(md);
		} else {
			cout << "[I] No data found for station " << ii+1 << " at date " << i_date.toString(Date::FULL) 
				<< endl;
			vecMeteo.push_back(MeteoData());
			vecMeteo[ii].date = i_date; //set correct date
		}
	}

	if (vecMeteo.size() == 0) {//No data found - return one object set to i_date and nodata in all other fields
		vecMeteo.push_back(MeteoData());
		vecMeteo[0].date = i_date; //set correct date
	}	
}

void BufferedIOHandler::getNextMeteoData(const Date& _date, std::vector<MeteoData>& vecMeteo)
{
	//TODO: check whether there is something in the buffer!
	//Try to rebuffer!

	vecMeteo.clear();
	
	std::vector< std::vector<MeteoData> > meteoTmpBuffer;
	readMeteoData(_date, (_date-Date(1900,1,2,0,0)), meteoTmpBuffer);	

	unsigned int emptycounter = 0;
	for (unsigned int ii=0; ii<meteoTmpBuffer.size(); ii++){//stations
		if (meteoTmpBuffer[ii].size() > 0){
			vecMeteo.push_back(meteoTmpBuffer[ii][0]);
		} else {
			emptycounter++;
		}
	}

	if (emptycounter == meteoTmpBuffer.size())
		vecMeteo.clear();
}

void BufferedIOHandler::bufferAllData(const Date& _date){
	Date fromDate = _date - bufferbefore;
	Date toDate   = _date + bufferafter;

	readMeteoData(fromDate, toDate, meteoBuffer);

	for (unsigned int ii=0; ii<meteoBuffer.size(); ii++){
		//set the start and the end date of the interval requested for each station
		startDateBuffer.push_back(fromDate);
		endDateBuffer.push_back(toDate);
	}
}

bool BufferedIOHandler::bufferData(const Date& _date, const unsigned int& stationindex)
{
	Date fromDate = _date - bufferbefore;
	Date toDate   = _date + bufferafter;

	readMeteoData(fromDate, toDate, meteoBuffer, stationindex);
	startDateBuffer.at(stationindex) = fromDate;
	endDateBuffer.at(stationindex) = toDate;

	if (meteoBuffer.size() == 0) {
		return false;
	}
	
	if (meteoBuffer[stationindex].size() == 0) {
		return false;
	}

	if ((!((_date >= meteoBuffer[stationindex][0].date) 
		  && (meteoBuffer[stationindex][meteoBuffer[stationindex].size()-1].date >= _date)))) {
		meteoBuffer[stationindex].clear();
		return false;
	}

	//If we reach this point: Date is definitely covered
	return true;
}

void BufferedIOHandler::setBufferPolicy(const buffer_policy& policy)
{
	if (policy==RECHECK_NODATA){
		always_rebuffer=true;
	} else {
		always_rebuffer=false;
	}
}

void BufferedIOHandler::setBufferDuration(const Date& _beforeDate, const Date& _afterDate)
{
	bufferbefore = _beforeDate; 
	bufferafter  = _afterDate;
}

void BufferedIOHandler::readMeteoData(const Date& dateStart, const Date& dateEnd, 
							   std::vector< std::vector<MeteoData> >& vecMeteo,
							   const unsigned int& stationindex)
	
{
	iohandler.readMeteoData(dateStart, dateEnd, vecMeteo, stationindex);

	if (&meteoBuffer != &vecMeteo)
		meteoBuffer = vecMeteo;      //copy by value
}

void BufferedIOHandler::readSpecialPoints(std::vector<Coords>& _cpa)
{
	iohandler.readSpecialPoints(_cpa);
}

void BufferedIOHandler::write2DGrid(const Grid2DObject& _grid2Dobj, const std::string& _name)
{
	iohandler.write2DGrid(_grid2Dobj, _name);
}

void BufferedIOHandler::clearBuffer(){
	meteoBuffer.clear();
	startDateBuffer.clear();
	endDateBuffer.clear();
	mapBufferedGrids.clear();
}

std::ostream& operator<<(std::ostream& os, const BufferedIOHandler& data)
{
	os << "<BufferedIOHandler>\n";
	os << "Config cfg; (not expanded)\n";
	#ifndef _POPC_
	os << data.iohandler;
	#else
	os << data.iohandler.toString();
	#endif
	os << data.meteoprocessor;

	os << "Rebuffer if not found: " << data.always_rebuffer << "\n";
	os << "Buffer span: -" << data.bufferbefore.getJulianDate() << " days, +" << data.bufferafter.getJulianDate() << " days\n";

	
	os << "Current buffer content (" << data.meteoBuffer.size() << " stations, " << data.mapBufferedGrids.size() << " grids):\n";
	for(unsigned int i=0; i<data.meteoBuffer.size(); i++) {
		os << std::setw(10) << data.meteoBuffer[i][0].meta.stationID << " = ";
		os << data.startDateBuffer[i].toString(Date::ISO) << " - ";
		os << data.endDateBuffer[i].toString(Date::ISO) << ", ";
		os << data.meteoBuffer[i].size() << " timesteps\n";
	}

	std::map<std::string, Grid2DObject>::const_iterator it1;
	for (it1=data.mapBufferedGrids.begin(); it1 != data.mapBufferedGrids.end(); it1++){
		os << setw(10) << "Grid" << " = " << it1->first << ", ";
		os << (it1->second).ncols << " x " << (it1->second).nrows << " @ " << (it1->second).cellsize << "m\n";
	}

	os << "</BufferedIOHandler>\n";

	return os;
}

} //namespace
