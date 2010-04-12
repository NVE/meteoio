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
#include "BufferedIOHandler.h"

using namespace std;

#ifdef _POPC_
BufferedIOHandler::BufferedIOHandler(IOHandler& _iohandler, const ConfigReader& _cfg) 
	: iohandler(_iohandler), cfg(_cfg), meteoFilter(_cfg), meteoBuffer(), stationBuffer(), startDateBuffer(), endDateBuffer(), mapBufferedGrids()
#else
BufferedIOHandler::BufferedIOHandler(IOHandler& _iohandler, const ConfigReader& _cfg) 
	  : IOInterface(NULL), iohandler(_iohandler), cfg(_cfg), meteoFilter(_cfg), meteoBuffer(), stationBuffer(), startDateBuffer(), endDateBuffer(), mapBufferedGrids()
#endif
{
	setBufferProperties();
}

#ifdef _POPC_
BufferedIOHandler::~BufferedIOHandler()
#else
BufferedIOHandler::~BufferedIOHandler() throw()
#endif
{
	setBufferProperties();
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
	if (it != mapBufferedGrids.end()) { //already in map
		_grid2Dobj = (*it).second; 
		return;
	}
	
	DEMObject tmpgrid2D;
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

void BufferedIOHandler::readAssimilationData(const Date_IO& _date, Grid2DObject& _grid2Dobj)
{
	std::map<std::string, Grid2DObject>::iterator it = mapBufferedGrids.find("/:ASSIMILATIONDATA" + _date.toString(Date_IO::FULL));
	if (it != mapBufferedGrids.end()) { //already in map
		_grid2Dobj = (*it).second; 
		return;
	}
	
	Grid2DObject tmpgrid2D;
	iohandler.readAssimilationData(_date, tmpgrid2D);
	mapBufferedGrids["/:ASSIMILATIONDATA" + _date.toString(Date_IO::FULL)] = tmpgrid2D;
	_grid2Dobj = tmpgrid2D;
}

void BufferedIOHandler::readStationData(const Date_IO& date, std::vector<StationData>& vecStation)
{
	iohandler.readStationData(date, vecStation);
}

#ifdef _POPC_
void BufferedIOHandler::writeMeteoData(std::vector< std::vector<MeteoData> >& vecMeteo,
							    std::vector< std::vector<StationData> >& vecStation,
							    const std::string& name)
#else 
void BufferedIOHandler::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo,
							    const std::vector< std::vector<StationData> >& vecStation,
							    const std::string& name)
#endif
{
	iohandler.writeMeteoData(vecMeteo, vecStation, name);
}

void BufferedIOHandler::readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo)
{
	std::vector< std::vector<StationData> > vecStation;
	readMeteoData(dateStart, dateEnd, vecMeteo, vecStation);
}

void BufferedIOHandler::setBufferProperties()
{
	always_rebuffer = true;
	bufferbefore = Date_IO(2.0);  //minus 2 days
	bufferafter = Date_IO(20.0);  //plus 20 days

	/*
	string bufferstrategy;
	vector<string> bufferperiod;
	cfg.getValue("BUFFERSTRATEGY", bufferstrategy, ConfigReader::nothrow);
	cfg.getValue("BUFFERPERIOD", bufferperiod, ConfigReader::nothrow);

	if (bufferperiod.size() != 2){
		bufferbefore = Date_IO(2.0);  //minus 2 days
		bufferafter = Date_IO(20.0);  //plus 20 days
	} else {
		double before=0.0, after=0.0;
		if (!IOUtils::convertString(before, bufferperiod[0], std::dec) || 
		    !IOUtils::convertString(after, bufferperiod[1], std::dec))
			throw ConversionFailedException("Key BUFFERPERIOD has invalid values", AT);
		
		if ((before<0) || (after<0))
			throw ConversionFailedException("Key BUFFERPERIOD has invalid values", AT);

		bufferbefore = Date_IO(before); 
		bufferafter  = Date_IO(after);
	}

	//buffer strategy
	always_rebuffer = true;
	if (bufferstrategy == "newinterval")
		always_rebuffer = false;
	*/
}

void BufferedIOHandler::readMeteoData(const Date_IO& date_in, std::vector<MeteoData>& vecMeteo, std::vector<StationData>& vecStation){
	/* For every station: 
	 * 1) See whether data is already buffered
	 * 2) Filter - includes resampling
	 * 3) Return the values
	 */

	vecMeteo.clear();
	vecStation.clear();

	if (meteoBuffer.size() == 0){ //init
		bufferAllData(date_in);
	}

	//loop through all meteo buffers, there is one for each station
	for (unsigned int ii=0; ii<meteoBuffer.size(); ii++) {
		unsigned int index = BufferedIOHandler::npos;

		StationData tmpsd;
		std::string stationName = "";
		if (stationBuffer[ii].size() > 0){
			stationName = stationBuffer[ii][0].stationName;
			tmpsd = stationBuffer[ii][0];
		}

		if (meteoBuffer[ii].size() > 0) {//check whether meteo data for the date exists in buffer
			index = seek(date_in, meteoBuffer[ii]);
		}

		if (index == BufferedIOHandler::npos) { //not in buffer
			//Check buffering strategy
			bool rebuffer = false;
			if ((startDateBuffer.at(ii) > date_in) || (endDateBuffer.at(ii) < date_in)){
				rebuffer = true;
			} else { 
				if (always_rebuffer) rebuffer = true;
			}

			if (rebuffer){
				cout << "[I] Station " << ii << "(" << stationName 
					<< ") data for date " << date_in.toString(Date_IO::FULL) << " not in buffer ..." << endl;
				
				bool dataexists = bufferData(date_in, ii);
				if (dataexists) {//date_in is contained in buffer
					index = seek(date_in, meteoBuffer[ii]);
				}
			}
		}

		//APPLY FILTERS
		MeteoData md; StationData sd;
		if (index != BufferedIOHandler::npos) {
			vector<MeteoData> mBuffer;
			std::vector<StationData> sBuffer;
			meteoFilter.filterData(meteoBuffer[ii], stationBuffer[ii], index, date_in, md, sd);
		}

		if (index != BufferedIOHandler::npos) {
			//vecMeteo.push_back(meteoBuffer[ii][index]);
			//vecStation.push_back(stationBuffer[ii][index]);
			vecMeteo.push_back(md);
			vecStation.push_back(sd);
		} else {
			cout << "[I] Buffering data for Station " << stationName << " at date " 
				<< date_in.toString(Date_IO::FULL) << " failed" << endl;
			vecMeteo.push_back(MeteoData());
			vecMeteo[ii].date = date_in; //set correct date

			vecStation.push_back(sd);
		}
	}

	if (vecMeteo.size() == 0) {//No data found - return one object set to date_in and nodata in all other fields
		vecMeteo.push_back(MeteoData());
		vecMeteo[0].date = date_in; //set correct date

		vecStation.push_back(StationData());
		//throw IOException("[E] No data for any station for date " + date_in.toString(Date_IO::FULL) + " found", AT);
	}	
}

void BufferedIOHandler::getNextMeteoData(const Date_IO& _date, std::vector<MeteoData>& vecMeteo, std::vector<StationData>& vecStation)
{
	//TODO: check whether there is something in the buffer!
	//Try to rebuffer!

	vecMeteo.clear();
	vecStation.clear();
	
	std::vector< std::vector<MeteoData> > meteoTmpBuffer;
	std::vector< std::vector<StationData> > stationTmpBuffer;
	readMeteoData(_date, (_date-Date_IO(1900,1,2)), meteoTmpBuffer, stationTmpBuffer);	

	unsigned int emptycounter = 0;
	for (unsigned int ii=0; ii<meteoTmpBuffer.size(); ii++){//stations
		if ((meteoTmpBuffer[ii].size() > 0) && (stationTmpBuffer[ii].size() > 0)){
			vecMeteo.push_back(meteoTmpBuffer[ii][0]);
			vecStation.push_back(stationTmpBuffer[ii][0]);
		} else {
			emptycounter++;
		}
	}

	if (emptycounter == meteoTmpBuffer.size()){
		vecMeteo.clear();
		vecStation.clear();
	}
}

void BufferedIOHandler::bufferAllData(const Date_IO& _date){
	Date_IO fromDate = _date - bufferbefore;
	Date_IO toDate   = _date + bufferafter;

	readMeteoData(fromDate, toDate, meteoBuffer, stationBuffer);

	for (unsigned int ii=0; ii<meteoBuffer.size(); ii++){
		//set the start and the end date of the interval requested for each station
		startDateBuffer.push_back(fromDate);
		endDateBuffer.push_back(toDate);
	}
}

bool BufferedIOHandler::bufferData(const Date_IO& _date, const unsigned int& stationindex)
{
	Date_IO fromDate = _date - bufferbefore;
	Date_IO toDate   = _date + bufferafter;

	readMeteoData(fromDate, toDate, meteoBuffer, stationBuffer, stationindex);
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

	//If we reach this point: Date_IO is definitely covered
	return true;
}

void BufferedIOHandler::bufferAlways(const bool& bufferalways)
{
	always_rebuffer = bufferalways;
}

void BufferedIOHandler::setBufferDuration(const Date_IO& _beforeDate, const Date_IO& _afterDate)
{
	bufferbefore = _beforeDate; 
	bufferafter  = _afterDate;
}

void BufferedIOHandler::readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							   std::vector< std::vector<MeteoData> >& vecMeteo,
							   std::vector< std::vector<StationData> >& vecStation,
							   const unsigned int& stationindex)
	
{
	iohandler.readMeteoData(dateStart, dateEnd, vecMeteo, vecStation, stationindex);
	if ((&meteoBuffer != &vecMeteo) && (&stationBuffer != &vecStation)){
		meteoBuffer = vecMeteo;      //copy by value
		stationBuffer = vecStation;  //copy by value
	}
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
	stationBuffer.clear();
	startDateBuffer.clear();
	endDateBuffer.clear();
	mapBufferedGrids.clear();
}

unsigned int BufferedIOHandler::seek(const Date_IO& date_in, std::vector<MeteoData>& vecM){ //TODO: binary search
	//returns index of element, if element does not exist it returns closest index
	//the element needs to be an exact hit or embedded between two measurments

	unsigned int ii = 1;

	if (vecM.size() <= 0) {//no elements in buffer
		return BufferedIOHandler::npos;
	}

	//if we reach this point: at least one element in buffer
	if (vecM[0].date > date_in) {
		return BufferedIOHandler::npos;
	}

	if (vecM[vecM.size()-1].date < date_in) {//last element is earlier, return npos
		return BufferedIOHandler::npos;
	}

	if (vecM[0].date == date_in) {//closest element
		return 0;
	}

	//if we reach this point: the date is spanned by the buffer and there are at least two elements
	while ((ii < vecM.size())) {
		if ((vecM[ii].date >= date_in) && (vecM[ii-1].date < date_in)) {
			return ii;
		}

		ii++;
	}

	return BufferedIOHandler::npos;
}
