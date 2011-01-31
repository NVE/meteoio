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
	: iohandler(_iohandler), cfg(_cfg), mapBufferedGrids()
#else
BufferedIOHandler::BufferedIOHandler(IOHandler& _iohandler, const Config& _cfg) 
	  : IOInterface(NULL), iohandler(_iohandler), cfg(_cfg), mapBufferedGrids()
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
	default_chunk_size = Date(15.0); //15 days
}

void BufferedIOHandler::setBufferPolicy(const buffer_policy& policy)
{
	if (policy==RECHECK_NODATA){
		always_rebuffer=true;
	} else {
		always_rebuffer=false;
	}
}

double BufferedIOHandler::getAvgSamplingRate()
{
	if (vec_buffer_meteo.size() > 0){
		unsigned int sum = 0;
		for (unsigned int ii=0; ii<vec_buffer_meteo.size(); ii++){
			//count all data
			sum += vec_buffer_meteo[ii].size();
		}
		if (sum > 0){
			double days = buffer_end.getJulianDate() - buffer_start.getJulianDate();
			return ((double)sum / days);
		}
	}

	return IOUtils::nodata;
}

const std::vector<METEO_DATASET>& BufferedIOHandler::get_complete_buffer(Date& start, Date& end)
{
	start = buffer_start;
	end   = buffer_end;

	return vec_buffer_meteo; //return reference
}

void BufferedIOHandler::readMeteoData(const Date& date_start, const Date& date_end, 
							   std::vector< std::vector<MeteoData> >& vecMeteo,
							   const unsigned int& /*stationindex*/)
	
{
	vecMeteo.clear();

	Date current_buffer_end(date_start + default_chunk_size);
	vector< vector<MeteoData> > tmp_meteo_buffer;

	//Read MeteoData for requested interval in chunks, furthermore buffer it
	//Try to buffer after the requested chunk for subsequent calls

	//0. initialize if not already initialized
	if (vec_buffer_meteo.size() == 0) //init
		bufferAllData(date_start, current_buffer_end);

	unsigned int buffer_size = vec_buffer_meteo.size();

	//1. Check whether data is in buffer already, and buffer it if not
	if ((date_start >= buffer_start) && (date_end <= buffer_end)){
		//copy data and we're done
	} else {
		//rebuffer data
		if (current_buffer_end == buffer_end){
			//only append
		} else {
			//rebuffer for real
			bufferAllData(date_start, current_buffer_end);
			buffer_size = vec_buffer_meteo.size();
		}
		
		while (date_end > current_buffer_end){
			iohandler.readMeteoData(current_buffer_end, current_buffer_end+default_chunk_size, tmp_meteo_buffer);			

			if (tmp_meteo_buffer.size() != buffer_size)
				throw IOException("God damn it!", AT);
			
			//Loop through stations and append data
			for (unsigned int ii=0; ii<buffer_size; ii++){
				unsigned int station_size = vec_buffer_meteo[ii].size();

				if ((station_size > 0) && (tmp_meteo_buffer[ii].size() > 0)){
					//check if the last element equals the first one
					if (vec_buffer_meteo[ii][station_size-1].date >= tmp_meteo_buffer[ii][0].date)
						vec_buffer_meteo[ii].pop_back(); //delete the element with the same date
				}

				vec_buffer_meteo[ii].insert(vec_buffer_meteo[ii].end(), tmp_meteo_buffer[ii].begin(), tmp_meteo_buffer[ii].end());
			}
			current_buffer_end += default_chunk_size;
			buffer_end = current_buffer_end;				
		}
	}

	//2. Copy appropriate data into vecMeteo
	for (unsigned int ii=0; ii<buffer_size; ii++){
		vecMeteo.push_back(vector<MeteoData>()); //insert one empty vector of MeteoData

		if (vec_buffer_meteo[ii].size() == 0) continue; //no data in buffer for this station

		unsigned int pos_start = IOUtils::seek(date_start, vec_buffer_meteo[ii], false);
		if (pos_start == IOUtils::npos) pos_start = 0;

		unsigned int pos_end = IOUtils::seek(date_end, vec_buffer_meteo[ii], false);//HACK:: edit IOUtils::seek to accept an offset
		if (pos_end == IOUtils::npos)	pos_end = vec_buffer_meteo[ii].size() - 1; //just copy until the end of the buffer
		//cout << "Station " << ii << ": pos_start=" << pos_start << "  pos_end=" << pos_end << endl; 
		if (vec_buffer_meteo[ii][pos_end].date > date_end){
			if (pos_end > pos_start)	pos_end--;
		} else {
			pos_end++;
		}
		//cout << "Station " << ii << ": pos_start=" << pos_start << "  pos_end=" << pos_end << endl; 
		vecMeteo[ii].insert(vecMeteo[ii].begin(), vec_buffer_meteo[ii].begin()+pos_start, vec_buffer_meteo[ii].begin()+pos_end);
	}
}

void BufferedIOHandler::bufferAllData(const Date& date_start, const Date& date_end){
	vec_buffer_meteo.clear();
	iohandler.readMeteoData(date_start, date_end, vec_buffer_meteo);
	buffer_start = date_start;
	buffer_end   = date_end;
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
	mapBufferedGrids.clear();
}

std::ostream& operator<<(std::ostream& os, const BufferedIOHandler& data)
{
	os << "<BufferedIOHandler>\n";
	os << "Config cfg = " << hex << &data.cfg << "\n";
	os << "IOHandler &iohandler = " << hex << &data.iohandler << "\n";

	os << "Rebuffer if not found: " << data.always_rebuffer << "\n";
	
	os << "Current buffer content (" << data.vec_buffer_meteo.size() << " stations, " 
	   << data.mapBufferedGrids.size() << " grids):\n";

	for(unsigned int ii=0; ii<data.vec_buffer_meteo.size(); ii++) {
		if (data.vec_buffer_meteo[ii].size() > 0){
			os << std::setw(10) << data.vec_buffer_meteo[ii][0].meta.stationID << " = "
			   << data.vec_buffer_meteo[ii][0].date.toString(Date::ISO) << " - "
			   << data.vec_buffer_meteo[ii][data.vec_buffer_meteo[ii].size()-1].date.toString(Date::ISO) << ", "
			   << data.vec_buffer_meteo[ii].size() << " timesteps" << endl;
		}
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
