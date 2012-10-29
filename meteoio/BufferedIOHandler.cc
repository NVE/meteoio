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
        : iohandler(in_iohandler), cfg(in_cfg), vec_buffer_meteo(), mapBufferedGrids(), IndexBufferedGrids(),
          buffer_start(), buffer_end(), chunk_size(), buff_before(), chunks(1), max_grids(10)
#else
BufferedIOHandler::BufferedIOHandler(IOHandler& in_iohandler, const Config& in_cfg)
        : IOInterface(NULL), iohandler(in_iohandler), cfg(in_cfg), vec_buffer_meteo(), mapBufferedGrids(), IndexBufferedGrids(),
          buffer_start(), buffer_end(), chunk_size(), buff_before(), chunks(1), max_grids(10)

#endif
{
	setDfltBufferProperties();
}

#ifdef _POPC_
BufferedIOHandler::~BufferedIOHandler()
#else
BufferedIOHandler::~BufferedIOHandler() throw()
#endif
{}

void BufferedIOHandler::bufferGrid(const Grid2DObject& in_grid2Dobj, const std::string& in_filename)
{
	if(IndexBufferedGrids.size() >= max_grids) { //we need to remove the oldest grid
		mapBufferedGrids.erase( mapBufferedGrids.find( IndexBufferedGrids[0] ) );
		IndexBufferedGrids.erase( IndexBufferedGrids.begin() );
	}
	mapBufferedGrids[in_filename] = in_grid2Dobj;
	IndexBufferedGrids.push_back( in_filename  );
}

void BufferedIOHandler::read2DGrid(Grid2DObject& in_grid2Dobj, const std::string& in_filename)
{
	if(max_grids>0) {
		std::map<std::string, Grid2DObject>::iterator it = mapBufferedGrids.find(in_filename);
		if (it != mapBufferedGrids.end()) { //already in map
			in_grid2Dobj = (*it).second;
			return;
		}

		Grid2DObject tmpgrid2D; //HACK: why don't we read directly into in_grid2Dobj?
		iohandler.read2DGrid(tmpgrid2D, in_filename);
		bufferGrid(tmpgrid2D, in_filename);
		in_grid2Dobj = tmpgrid2D;
	} else {
		iohandler.read2DGrid(in_grid2Dobj, in_filename);
	}
}

void BufferedIOHandler::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if(max_grids>0) {
		const string buffer_name = date.toString(Date::ISO)+"::"+MeteoGrids::getParameterName(parameter);
		std::map<std::string, Grid2DObject>::iterator it = mapBufferedGrids.find(buffer_name);
		if (it != mapBufferedGrids.end()) { //already in map
			grid_out = (*it).second;
			return;
		}

		Grid2DObject tmpgrid2D;
		iohandler.read2DGrid(tmpgrid2D, parameter, date);
		bufferGrid(tmpgrid2D, buffer_name);
		grid_out = tmpgrid2D;
	} else {
		iohandler.read2DGrid(grid_out, parameter, date);
	}
}

void BufferedIOHandler::readDEM(DEMObject& in_grid2Dobj)
{
	if(max_grids>0) {
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
	} else {
		iohandler.readDEM(in_grid2Dobj);
	}
}

void BufferedIOHandler::readLanduse(Grid2DObject& in_grid2Dobj)
{
	if(max_grids>0) {
		std::map<std::string, Grid2DObject>::iterator it = mapBufferedGrids.find("/:LANDUSE");
		if (it != mapBufferedGrids.end()) { //already in map
			in_grid2Dobj = (*it).second;
			return;
		}

		Grid2DObject tmpgrid2D;
		iohandler.readLanduse(tmpgrid2D);
		mapBufferedGrids["/:LANDUSE"] = tmpgrid2D;
		in_grid2Dobj = tmpgrid2D;
	} else {
		iohandler.readLanduse(in_grid2Dobj);
	}
}

//HACK: manage buffering of assimilation grids! Why not considering them normal grids?
void BufferedIOHandler::readAssimilationData(const Date& in_date, Grid2DObject& in_grid2Dobj)
{
	if(max_grids>0) {
		std::map<std::string, Grid2DObject>::iterator it = mapBufferedGrids.find("/:ASSIMILATIONDATA" + in_date.toString(Date::FULL));
		if (it != mapBufferedGrids.end()) { //already in map
			in_grid2Dobj = (*it).second;
			return;
		}

		Grid2DObject tmpgrid2D;
		iohandler.readAssimilationData(in_date, tmpgrid2D);
		mapBufferedGrids["/:ASSIMILATIONDATA" + in_date.toString(Date::FULL)] = tmpgrid2D;
		in_grid2Dobj = tmpgrid2D;
	} else {
		iohandler.readAssimilationData(in_date, in_grid2Dobj);
	}
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
	double chunk_size_days = 15.; //default chunk size value
	chunks = 1;
	cfg.getValue("BUFF_CHUNK_SIZE", "General", chunk_size_days, Config::nothrow); //in days
	cfg.getValue("BUFF_CHUNKS", "General", chunks, Config::nothrow);
	chunk_size = Duration(chunk_size_days, 0);

	//get buffer centering options
	double buff_centering = -1.;
	double buff_start = -1.;
	cfg.getValue("BUFF_CENTERING", "General", buff_centering, Config::nothrow);
	cfg.getValue("BUFF_BEFORE", "General", buff_start, Config::nothrow);
	if ((buff_centering != -1.) && (buff_start != -1.))
		throw InvalidArgumentException("Please do NOT provide both BUFF_CENTERING and BUFF_BEFORE!!", AT);

	if (buff_start != -1.){
		buff_before = Duration(buff_start, 0);
	} else {
		if (buff_centering != -1.){
			if ((buff_centering < 0.) || (buff_centering > 1.))
				throw InvalidArgumentException("BUFF_CENTERING must be between 0 and 1", AT);

			buff_before = chunk_size * buff_centering;
		} else {
			buff_before = chunk_size * 0.1; //10% centering by default
		}
	}

	//if buff_before>chunk_size, we will have a problem (ie: we won't ever read the whole data we need)
	if(buff_before>chunk_size) chunk_size = buff_before;
	//BUG: if we do this, we still have the meteo1d window in the way
	//-> we end up not reading enough data and rebuffering...

	max_grids = 10; //default number of grids to keep in buffer
	cfg.getValue("BUFF_GRIDS", "General", max_grids, Config::nothrow);
}

void BufferedIOHandler::setMinBufferRequirements(const double& i_chunk_size, const double& i_buff_before)
{
	if(i_buff_before!=IOUtils::nodata) {
		const Duration app_buff_before(i_buff_before, 0);
		if(app_buff_before>buff_before) buff_before = app_buff_before;
	}
	if(i_chunk_size!=IOUtils::nodata) {
		const Duration app_chunk_size(i_chunk_size, 0);
		if(app_chunk_size>chunk_size) chunk_size = app_chunk_size;
	}

	//if buff_before>chunk_size, we will have a problem (ie: we won't ever read the whole data we need)
	if(buff_before>chunk_size) chunk_size = buff_before;
}

double BufferedIOHandler::getAvgSamplingRate()
{
	if (vec_buffer_meteo.size() > 0){
		unsigned long sum = 0;
		for (size_t ii=0; ii<vec_buffer_meteo.size(); ii++){
			//count all data
			sum += (unsigned long)vec_buffer_meteo[ii].size();
		}
		if (sum > 0){
			double days = buffer_end.getJulian() - buffer_start.getJulian();
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
                                      const size_t& /*stationindex*/)
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

	size_t buffer_size = vec_buffer_meteo.size();

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
			tmp_meteo_buffer.reserve(buffer_size);
			iohandler.readMeteoData(new_buffer_end, new_buffer_end+chunk_size*chunks, tmp_meteo_buffer);

			if (tmp_meteo_buffer.size() != buffer_size) {
				stringstream ss;
				ss << "The number of stations changed over time from " << buffer_size << " to " << tmp_meteo_buffer.size() << ", ";
				ss << "this is not handled yet!";
				throw IOException(ss.str(), AT);
			}

			//Loop through stations and append data
			for (size_t ii=0; ii<buffer_size; ii++){ //loop through stations
				const size_t station_size = vec_buffer_meteo[ii].size();

				if ((station_size > 0) && (tmp_meteo_buffer[ii].size() > 0)){
					//check if the last element equals the first one
					if (vec_buffer_meteo[ii][station_size-1].date >= tmp_meteo_buffer[ii][0].date)
						vec_buffer_meteo[ii].pop_back(); //delete the element with the same date
				}

				vec_buffer_meteo[ii].reserve(vec_buffer_meteo[ii].size()+tmp_meteo_buffer[ii].size());
				vec_buffer_meteo[ii].insert(vec_buffer_meteo[ii].end(), tmp_meteo_buffer[ii].begin(), tmp_meteo_buffer[ii].end());
			}
			new_buffer_end += chunk_size*chunks;
			buffer_end = new_buffer_end;
		}
	}

	//2. Copy appropriate data into vecMeteo
	vecMeteo.reserve(buffer_size);
	for (size_t ii=0; ii<buffer_size; ii++){ //loop through stations
		vecMeteo.push_back(vector<MeteoData>()); //insert one empty vector of MeteoData

		if (vec_buffer_meteo[ii].size() == 0) continue; //no data in buffer for this station

		size_t pos_start = IOUtils::seek(date_start, vec_buffer_meteo[ii], false);
		if (pos_start == IOUtils::npos) pos_start = 0;

		size_t pos_end = IOUtils::seek(date_end, vec_buffer_meteo[ii], false);//HACK:: edit IOUtils::seek to accept an offset
		if (pos_end == IOUtils::npos) pos_end = vec_buffer_meteo[ii].size() - 1; //just copy until the end of the buffer

		if (vec_buffer_meteo[ii][pos_end].date > date_end){
			if (pos_end > pos_start) pos_end--;
		} else {
			pos_end++;
		}
		vecMeteo[ii].reserve(pos_end-pos_start+1); //weird that the "insert" does not handle it internally...
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

/**
 * @brief Push a vector of time series of MeteoData objects into the local buffer.
 *        This overwrites the local buffer. This method is a way to bypass the internal reading
 *        of MeteoData from a certain source.
 * @param date_start Representing the beginning of the data
 * @param date_end Representing the end of the data
 * @param vecMeteo The actual data being pushed into vec_buffer_meteo
 */
void BufferedIOHandler::push_meteo_data(const Date& date_start, const Date& date_end,
                                        const std::vector< METEO_TIMESERIE >& vecMeteo)
{
	//perform check on date_start and date_end
	if (date_end < date_start)
		throw InvalidArgumentException("date_start cannot be greater than date_end", AT);

	buffer_start     = date_start;
	buffer_end       = date_end;
	vec_buffer_meteo = vecMeteo;
}

void BufferedIOHandler::readSpecialPoints(std::vector<Coords>& in_cpa)
{
	iohandler.readSpecialPoints(in_cpa);
}

void BufferedIOHandler::write2DGrid(const Grid2DObject& grid_in, const std::string& in_name)
{
	iohandler.write2DGrid(grid_in, in_name);
}

void BufferedIOHandler::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	iohandler.write2DGrid(grid_in, parameter, date);
}

void BufferedIOHandler::clearBuffer(){
	mapBufferedGrids.clear();
}

std::ostream& operator<<(std::ostream& os, const BufferedIOHandler& data)
{
	os << "<BufferedIOHandler>\n";
	os << "Config& cfg = " << hex << &data.cfg << dec << "\n";
	os << "IOHandler &iohandler = " << hex << &data.iohandler << dec << "\n";

	os << "Buffering " << data.chunks << " chunk(s) of " <<data.chunk_size.getJulian() << " day(s) with "
	   << data.buff_before.getJulian() << " day(s) pre-buffering\n";

	os << "Current buffer content (" << data.vec_buffer_meteo.size() << " stations, "
	   << data.mapBufferedGrids.size() << " grids):\n";

	for(size_t ii=0; ii<data.vec_buffer_meteo.size(); ii++) {
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
