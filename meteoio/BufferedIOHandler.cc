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
BufferedIOHandler::BufferedIOHandler(IOHandler& in_iohandler, const Config& in_cfg) 
	: iohandler(in_iohandler), cfg(in_cfg), mapBufferedGrids()
#else
BufferedIOHandler::BufferedIOHandler(IOHandler& in_iohandler, const Config& in_cfg) 
	  : IOInterface(NULL), iohandler(in_iohandler), cfg(in_cfg), mapBufferedGrids()
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

void BufferedIOHandler::read2DGrid(Grid2DObject& in_grid2Dobj, const std::string& in_filename)
{
	std::map<std::string, Grid2DObject>::iterator it = mapBufferedGrids.find(in_filename);
	if (it != mapBufferedGrids.end()) { //already in map
		in_grid2Dobj = (*it).second; 
		return;
	}
	
	Grid2DObject tmpgrid2D;
	iohandler.read2DGrid(tmpgrid2D, in_filename);
	mapBufferedGrids[in_filename] = tmpgrid2D;
	in_grid2Dobj = tmpgrid2D;
}

void BufferedIOHandler::readDEM(DEMObject& in_grid2Dobj)
{
	std::map<std::string, Grid2DObject>::iterator it = mapBufferedGrids.find("/:DEM");
	if (it != mapBufferedGrids.end()) {
		//already in map. If the update properties have changed,
		//we copy the ones given in input and force the update of the object
		const DEMObject::update_type in_ppt = (DEMObject::update_type)in_grid2Dobj.getUpdatePpt();
		in_grid2Dobj = (*it).second;
		const DEMObject::update_type buff_ppt = (DEMObject::update_type)in_grid2Dobj.getUpdatePpt();
		if(in_ppt!=buff_ppt) {
			in_grid2Dobj.setUpdatePpt(in_ppt);
			in_grid2Dobj.update();
		}
		return;
	}
	
	DEMObject tmpgrid2D;
	 //copy the updating policy of the destination
	tmpgrid2D.setUpdatePpt((DEMObject::update_type)in_grid2Dobj.getUpdatePpt());
	iohandler.readDEM(tmpgrid2D);
	mapBufferedGrids["/:DEM"] = tmpgrid2D;
	in_grid2Dobj = tmpgrid2D;
}

void BufferedIOHandler::readLanduse(Grid2DObject& in_grid2Dobj)
{
	std::map<std::string, Grid2DObject>::iterator it = mapBufferedGrids.find("/:LANDUSE");
	if (it != mapBufferedGrids.end()) { //already in map
		in_grid2Dobj = (*it).second; 
		return;
	}
	
	Grid2DObject tmpgrid2D;
	iohandler.readLanduse(tmpgrid2D);
	mapBufferedGrids["/:LANDUSE"] = tmpgrid2D;
	in_grid2Dobj = tmpgrid2D;
}

void BufferedIOHandler::readAssimilationData(const Date& in_date, Grid2DObject& in_grid2Dobj)
{
	std::map<std::string, Grid2DObject>::iterator it = mapBufferedGrids.find("/:ASSIMILATIONDATA" + in_date.toString(Date::FULL));
	if (it != mapBufferedGrids.end()) { //already in map
		in_grid2Dobj = (*it).second; 
		return;
	}
	
	Grid2DObject tmpgrid2D;
	iohandler.readAssimilationData(in_date, tmpgrid2D);
	mapBufferedGrids["/:ASSIMILATIONDATA" + in_date.toString(Date::FULL)] = tmpgrid2D;
	in_grid2Dobj = tmpgrid2D;
}

void BufferedIOHandler::readStationData(const Date& date, STATION_TIMESERIE& vecStation)
{
	iohandler.readStationData(date, vecStation);
}

#ifdef _POPC_
void BufferedIOHandler::writeMeteoData(std::vector< METEO_TIMESERIE >& vecMeteo,
                                       const std::string& name)
#else 
void BufferedIOHandler::writeMeteoData(const std::vector< METEO_TIMESERIE >& vecMeteo,
                                       const std::string& name)
#endif
{
	iohandler.writeMeteoData(vecMeteo, name);
}

void BufferedIOHandler::setDfltBufferProperties()
{
	always_rebuffer = false;

	double chunk_size_days = 15.; //default chunk size value
	chunks=1;
	cfg.getValue("BUFF_CHUNK_SIZE", "General", chunk_size_days,Config::nothrow); //in days
	cfg.getValue("BUFF_CHUNKS", "General", chunks,Config::nothrow);
	chunk_size = Duration(chunk_size_days, 0);

	//get buffer centering options
	double buff_centering = -1.;
	double buff_start = -1.;
	cfg.getValue("BUFF_CENTERING", "General", buff_centering, Config::nothrow);
	cfg.getValue("BUFF_BEFORE", "General", buff_start, Config::nothrow);
	if(buff_centering!=-1. && buff_start!=-1.) {
		throw InvalidArgumentException("Please do NOT provide both BUFF_CENTERING and BUFF_BEFORE!!", AT);
	}

	if(buff_start!=-1.) {
		buff_before = Duration(buff_start, 0);
	} else {
		if(buff_centering!=-1.) {
			if(buff_centering<0. || buff_centering>1.) {
				throw InvalidArgumentException("BUFF_CENTERING must be between 0 and 1", AT);
			}
			buff_before = chunk_size * buff_centering;
		} else {
			buff_before = chunk_size * 0.1; //10% centering by default
		}
	}
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

const std::vector< METEO_TIMESERIE >& BufferedIOHandler::get_complete_buffer(Date& start, Date& end)
{
	start = buffer_start;
	end   = buffer_end;

	return vec_buffer_meteo; //return reference
}

void BufferedIOHandler::readMeteoData(const Date& date_start, const Date& date_end, 
                                      std::vector< METEO_TIMESERIE >& vecMeteo,
                                      const unsigned int& /*stationindex*/)
{
	vecMeteo.clear();
	const Date new_buffer_start(date_start-buff_before); //taking centering into account
	Date new_buffer_end(new_buffer_start + chunk_size*chunks);
	vector< vector<MeteoData> > tmp_meteo_buffer; //it must be here -> adresses copied in 2. are still valid

	//Read MeteoData for requested interval in chunks, furthermore buffer it
	//Try to buffer after the requested chunk for subsequent calls

	//0. initialize if not already initialized
	if (vec_buffer_meteo.size() == 0) //init
		bufferData(new_buffer_start, new_buffer_end, vec_buffer_meteo);

	unsigned int buffer_size = vec_buffer_meteo.size();

	//1. Check whether data is in buffer already, and buffer it if not
	if ((date_start < buffer_start) || (date_end > buffer_end)){
		//rebuffer data
		if ((new_buffer_end != buffer_end) || (new_buffer_start != buffer_start)){
			//rebuffer for real
			bufferData(new_buffer_start, new_buffer_end, vec_buffer_meteo);
			buffer_size = vec_buffer_meteo.size();
		}

		while (date_end > new_buffer_end){
			//if the requested interval is bigger than a normal buffer, we have to increase the buffer anyway...
			iohandler.readMeteoData(new_buffer_end, new_buffer_end+chunk_size*chunks, tmp_meteo_buffer);

			if (tmp_meteo_buffer.size() != buffer_size)
				throw IOException("The number of stations changed over time, this is not handled yet!", AT);
			
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
			new_buffer_end += chunk_size*chunks;
			buffer_end = new_buffer_end;
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

void BufferedIOHandler::bufferData(const Date& date_start, const Date& date_end, std::vector< METEO_TIMESERIE >& vecvecMeteo){
	//TODO: implement reading data by chunks. It has to be done the same way as rebuffering
	vecvecMeteo.clear(); //the plugins do it internally anyway, but this is cheap and safe...
	iohandler.readMeteoData(date_start, date_end, vecvecMeteo);
	buffer_start = date_start;
	buffer_end   = date_end;
}

void BufferedIOHandler::readSpecialPoints(std::vector<Coords>& in_cpa)
{
	iohandler.readSpecialPoints(in_cpa);
}

void BufferedIOHandler::write2DGrid(const Grid2DObject& in_grid2Dobj, const std::string& in_name)
{
	iohandler.write2DGrid(in_grid2Dobj, in_name);
}

void BufferedIOHandler::clearBuffer(){
	mapBufferedGrids.clear();
}

std::ostream& operator<<(std::ostream& os, const BufferedIOHandler& data)
{
	os << "<BufferedIOHandler>\n";
	os << "Config& cfg = " << hex << &data.cfg << dec << "\n";
	os << "IOHandler &iohandler = " << hex << &data.iohandler << dec << "\n";

	os << "Rebuffer if not found: " << data.always_rebuffer << "\n";
	os << "Buffering " << data.chunks << " chunk(s) of " <<data.chunk_size.getJulianDate() << " days\n";
	
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
