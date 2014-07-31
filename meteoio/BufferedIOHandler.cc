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

BufferedIOHandler::BufferedIOHandler(IOHandler& in_iohandler, const Config& in_cfg)
	: iohandler(in_iohandler), cfg(in_cfg), meteo_buffer(), mapBufferedGrids(), dem_buffer(),
       IndexBufferedGrids(), buffer_start(), buffer_end(), chunk_size(), buff_before(), max_grids(10)
{
	setDfltBufferProperties();
}

BufferedIOHandler& BufferedIOHandler::operator=(const BufferedIOHandler& source) {
	if(this != &source) {
		iohandler = source.iohandler;
		meteo_buffer = source.meteo_buffer;
		mapBufferedGrids = source.mapBufferedGrids;
		dem_buffer = source.dem_buffer;
		IndexBufferedGrids = source.IndexBufferedGrids;
		buffer_start = source.buffer_start;
		buffer_end = source.buffer_end;
		chunk_size = source.chunk_size;
		buff_before = source.buff_before;
		max_grids = source.max_grids;
	}
	return *this;
}

void BufferedIOHandler::addToBuffer(const Grid2DObject& in_grid2Dobj, const std::string& grid_hash)
{
	if (max_grids==0) return;

	if(IndexBufferedGrids.size() >= max_grids) { //we need to remove the oldest grid
		mapBufferedGrids.erase( mapBufferedGrids.find( IndexBufferedGrids.front() ) );
		IndexBufferedGrids.erase( IndexBufferedGrids.begin() );
	}
	mapBufferedGrids[ grid_hash ] = in_grid2Dobj;
	IndexBufferedGrids.push_back( grid_hash  );
}

bool BufferedIOHandler::getFromBuffer(const std::string& grid_hash, Grid2DObject& grid) const
{
	if (IndexBufferedGrids.empty())
		return false;

	const std::map<std::string, Grid2DObject>::const_iterator it = mapBufferedGrids.find( grid_hash );
	if (it != mapBufferedGrids.end()) { //already in map
		grid = (*it).second;
		return true;
	}

	return false;
}

void BufferedIOHandler::read2DGrid(Grid2DObject& in_grid2Dobj, const std::string& in_filename)
{
	if (getFromBuffer(in_filename, in_grid2Dobj))
		return;

	iohandler.read2DGrid(in_grid2Dobj, in_filename);
	addToBuffer(in_grid2Dobj, in_filename);
}

void BufferedIOHandler::read2DGrid(Grid2DObject& in_grid2Dobj, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if (max_grids>0) {
		const string grid_hash = date.toString(Date::ISO)+"::"+MeteoGrids::getParameterName(parameter);
		if (getFromBuffer(grid_hash, in_grid2Dobj))
			return;

		iohandler.read2DGrid(in_grid2Dobj, parameter, date);
		addToBuffer(in_grid2Dobj, grid_hash); //the STL containers make a copy
	} else {
		iohandler.read2DGrid(in_grid2Dobj, parameter, date);
	}
}

void BufferedIOHandler::readDEM(DEMObject& demobj)
{
	if (max_grids>0) {
		if (dem_buffer.size() == 1) {
			//already in buffer. If the update properties have changed,
			//we copy the ones given in input and force the update of the object
			const DEMObject::update_type in_ppt = (DEMObject::update_type)demobj.getUpdatePpt();
			const DEMObject::slope_type in_slope_alg = (DEMObject::slope_type)demobj.getDefaultAlgorithm();

			demobj = dem_buffer[0];
			const DEMObject::update_type buff_ppt = (DEMObject::update_type)demobj.getUpdatePpt();
			const DEMObject::slope_type buff_slope_alg = (DEMObject::slope_type)demobj.getDefaultAlgorithm();

			if (in_ppt!=buff_ppt || in_slope_alg!=buff_slope_alg) {
				demobj.setDefaultAlgorithm(in_slope_alg);
				demobj.setUpdatePpt(in_ppt);
				demobj.update();
			}

			return;
		}

		iohandler.readDEM(demobj);
		dem_buffer.push_back(demobj); //the STL containers make a copy
	} else {
		iohandler.readDEM(demobj);
	}
}

void BufferedIOHandler::readLanduse(Grid2DObject& in_grid2Dobj)
{
	if (getFromBuffer("/:LANDUSE", in_grid2Dobj))
		return;

	iohandler.readLanduse(in_grid2Dobj);
	addToBuffer(in_grid2Dobj, "/:LANDUSE");
}

//HACK: manage buffering of assimilation grids! Why not considering them normal grids?
void BufferedIOHandler::readAssimilationData(const Date& date, Grid2DObject& in_grid2Dobj)
{
	if(max_grids>0) {
		const string grid_hash = "/:ASSIMILATIONDATA"+date.toString(Date::ISO);
		if (getFromBuffer(grid_hash, in_grid2Dobj))
			return;

		iohandler.readAssimilationData(date, in_grid2Dobj);
		addToBuffer(in_grid2Dobj, grid_hash); //the STL containers make a copy
	} else {
		iohandler.readAssimilationData(date, in_grid2Dobj);
	}
}

void BufferedIOHandler::readStationData(const Date& date, STATIONS_SET& vecStation)
{
	iohandler.readStationData(date, vecStation);
}

void BufferedIOHandler::writeMeteoData(const std::vector< METEO_SET >& vecMeteo,
                                       const std::string& name)
{
	iohandler.writeMeteoData(vecMeteo, name);
}

void BufferedIOHandler::setDfltBufferProperties()
{
	double chunk_size_days = 15.; //default chunk size value
	cfg.getValue("BUFF_CHUNK_SIZE", "General", chunk_size_days, IOUtils::nothrow); //in days
	chunk_size = Duration(chunk_size_days, 0);

	//get buffer centering options
	double buff_centering = -1.;
	double buff_start = -1.;
	cfg.getValue("BUFF_CENTERING", "General", buff_centering, IOUtils::nothrow);
	cfg.getValue("BUFF_BEFORE", "General", buff_start, IOUtils::nothrow);
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
	cfg.getValue("BUFF_GRIDS", "General", max_grids, IOUtils::nothrow);
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

double BufferedIOHandler::getAvgSamplingRate() const
{
	if (!meteo_buffer.empty()){
		const size_t nr_stations = meteo_buffer.size();
		double sum = 0;
		for (size_t ii=0; ii<nr_stations; ii++){ //loop over all stations
			if(!meteo_buffer[ii].empty()) {
				const std::vector<MeteoData>& curr_station = meteo_buffer[ii];
				const double days = curr_station.back().date.getJulian() - curr_station.front().date.getJulian();

				//add the average sampling rate for this station
				const size_t nr_data_pts = meteo_buffer[ii].size();
				if(days>0.) sum += (double)(nr_data_pts-1) / days; //the interval story: 2 points define 1 interval!
			}
		}
		if (sum > 0.){
			return ((double)sum / (double)(nr_stations*24*3600)); //in points per seconds, ie Hz
		}
	}

	return IOUtils::nodata;
}

/**
 * @brief return all the buffered data as well as the start and end dates
 * @param start start date of the buffer
 * @param end end date of the buffer
 * @return complete buffer
 */
const std::vector< METEO_SET >& BufferedIOHandler::getFullBuffer(Date& start, Date& end)
{
	start = buffer_start;
	end   = buffer_end;

	return meteo_buffer; //return reference
}

/**
 * @brief return all the buffered data between the given dates
 * @param date_start requested start date of the buffer
 * @param date_end requested end date of the buffer
 * @param data vector to fill with the buffered data
 */
void BufferedIOHandler::getFromBuffer(const Date& date_start, const Date& date_end, std::vector< METEO_SET > &vecMeteo)
{
	//1. Prepare the output vector
	const size_t buffer_size = meteo_buffer.size();
	vecMeteo.clear();
	vecMeteo.reserve(buffer_size);

	//2. Copy appropriate data into vecMeteo for each station
	for (size_t ii=0; ii<buffer_size; ii++){ //loop through stations
		vecMeteo.push_back(vector<MeteoData>()); //insert one empty vector of MeteoData

		if (meteo_buffer[ii].empty()) continue; //no data in buffer for this station

		size_t pos_start = IOUtils::seek(date_start, meteo_buffer[ii], false);
		if (pos_start == IOUtils::npos) pos_start = 0;

		size_t pos_end = IOUtils::seek(date_end, meteo_buffer[ii], false);//HACK:: edit IOUtils::seek to accept an offset
		if (pos_end == IOUtils::npos) pos_end = meteo_buffer[ii].size() - 1; //just copy until the end of the buffer

		if (meteo_buffer[ii][pos_end].date > date_end){
			if (pos_end > pos_start) pos_end--;
		} else {
			pos_end++;
		}
		vecMeteo[ii].reserve(pos_end-pos_start+1); //weird that the "insert" does not handle it internally...
		vecMeteo[ii].insert(vecMeteo[ii].begin(), meteo_buffer[ii].begin()+pos_start, meteo_buffer[ii].begin()+pos_end);
	}
}

void BufferedIOHandler::fillBuffer(const Date& date_start, const Date& date_end,
                                   const size_t& /*stationindex*/)
{
	const Date new_buffer_start(date_start-buff_before); //taking centering into account
	Date new_buffer_end(new_buffer_start + chunk_size);

	//Read MeteoData for requested interval in chunks, furthermore buffer it
	//Try to buffer after the requested chunk for subsequent calls

	//0. initialize if not already initialized
	if (meteo_buffer.empty()) {
		iohandler.readMeteoData(new_buffer_start, new_buffer_end, meteo_buffer);
		buffer_start = new_buffer_start;
		buffer_end   = new_buffer_end;
	}

	//1. Check whether data is in buffer already, and buffer it if not
	if ((date_start < buffer_start) || (date_end > buffer_end)) {
		//rebuffer data
		if ((new_buffer_end != buffer_end) || (new_buffer_start != buffer_start)) { //rebuffer for real
			meteo_buffer.clear(); //the plugins do it internally anyway, but this is cheap and safe...
			iohandler.readMeteoData(new_buffer_start, new_buffer_end, meteo_buffer);
			buffer_start = new_buffer_start;
			buffer_end   = new_buffer_end;
		}

		const size_t buffer_size = meteo_buffer.size();
		vector< vector<MeteoData> > tmp_meteo_buffer;
		while (date_end > new_buffer_end){
			//if the requested interval is bigger than a normal buffer, we have to increase the buffer anyway...
			tmp_meteo_buffer.reserve(buffer_size);
			iohandler.readMeteoData(new_buffer_end, new_buffer_end+chunk_size, tmp_meteo_buffer);

			if (tmp_meteo_buffer.size() != buffer_size) {
				ostringstream ss;
				ss << "The number of stations changed over time from " << buffer_size << " to " << tmp_meteo_buffer.size() << ", ";
				ss << "this is not handled yet!";
				throw IOException(ss.str(), AT);
			}

			//Loop through stations and append data
			for (size_t ii=0; ii<buffer_size; ii++){ //loop through stations
				if ((!meteo_buffer[ii].empty()) && (!tmp_meteo_buffer[ii].empty())){
					//check if the last element equals the first one
					if (meteo_buffer[ii].back().date >= tmp_meteo_buffer[ii].front().date)
						meteo_buffer[ii].pop_back(); //delete the element with the same date
				}

				meteo_buffer[ii].reserve(meteo_buffer[ii].size()+tmp_meteo_buffer[ii].size());
				meteo_buffer[ii].insert(meteo_buffer[ii].end(), tmp_meteo_buffer[ii].begin(), tmp_meteo_buffer[ii].end());
			}
			new_buffer_end += chunk_size;
			buffer_end = new_buffer_end;
		}
	}
}

void BufferedIOHandler::readMeteoData(const Date& date_start, const Date& date_end,
                                      std::vector< METEO_SET >& vecMeteo,
                                      const size_t& stationindex)
{
	fillBuffer(date_start, date_end, stationindex);
	getFromBuffer(date_start, date_end, vecMeteo);
}

/**
 * @brief Push a vector of time series of MeteoData objects into the local buffer.
 *        This overwrites the local buffer. This method is a way to bypass the internal reading
 *        of MeteoData from a certain source.
 * @param date_start Representing the beginning of the data
 * @param date_end Representing the end of the data
 * @param vecMeteo The actual data being pushed into meteo_buffer
 */
void BufferedIOHandler::push_meteo_data(const Date& date_start, const Date& date_end,
                                        const std::vector< METEO_SET >& vecMeteo)
{
	//perform check on date_start and date_end
	if (date_end < date_start)
		throw InvalidArgumentException("date_start cannot be greater than date_end", AT);

	buffer_start     = date_start;
	buffer_end       = date_end;
	meteo_buffer = vecMeteo;
}

void BufferedIOHandler::readPOI(std::vector<Coords>& in_cpa)
{
	iohandler.readPOI(in_cpa);
}

void BufferedIOHandler::write2DGrid(const Grid2DObject& grid_in, const std::string& in_name)
{
	iohandler.write2DGrid(grid_in, in_name);
}

void BufferedIOHandler::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	iohandler.write2DGrid(grid_in, parameter, date);
}

void BufferedIOHandler::clearBuffer() {
	meteo_buffer.clear();
	buffer_start = Date(0., 0.);
	buffer_end = Date(0., 0.);
	mapBufferedGrids.clear();
	IndexBufferedGrids.clear();
}

const std::string BufferedIOHandler::toString() const
{
	std::ostringstream os;
	os << "<BufferedIOHandler>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";
	os << "IOHandler &iohandler = " << hex << &iohandler << dec << "\n";

	os << "Buffering " <<chunk_size.getJulian() << " day(s) with "
	   << buff_before.getJulian() << " day(s) pre-buffering\n";

	os << "Current buffer content (" << meteo_buffer.size() << " stations, "
	   << mapBufferedGrids.size() << " grids):\n";

	for(size_t ii=0; ii<meteo_buffer.size(); ii++) {
		if (!meteo_buffer[ii].empty()){
			os << std::setw(10) << meteo_buffer[ii].front().meta.stationID << " = "
			   << meteo_buffer[ii].front().date.toString(Date::ISO) << " - "
			   << meteo_buffer[ii].back().date.toString(Date::ISO) << ", "
			   << meteo_buffer[ii].size() << " timesteps" << endl;
		}
	}

	std::map<std::string, Grid2DObject>::const_iterator it1;
	for (it1=mapBufferedGrids.begin(); it1 != mapBufferedGrids.end(); ++it1){
		os << setw(10) << "Grid" << " = " << it1->first << ", ";
		os << (it1->second).ncols << " x " << (it1->second).nrows << " @ " << (it1->second).cellsize << "m\n";
	}

	os << "</BufferedIOHandler>\n";

	return os.str();
}

} //namespace
