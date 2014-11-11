/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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

#include <meteoio/TimeSeriesManager.h>

using namespace std;

namespace mio {

TimeSeriesManager::TimeSeriesManager(IOHandler& in_iohandler, const Config& in_cfg) : cfg(in_cfg), iohandler(in_iohandler),
                                            meteoprocessor(in_cfg), dataGenerator(in_cfg),
                                            proc_properties(), point_cache(), filtered_cache(), raw_buffer(),
                                            fcache_start(Date(0.0, 0.)), fcache_end(Date(0.0, 0.)), //this should not matter, since 0 is still way back before any real data...
                                            raw_start(Date(0.0, 0.)), raw_end(Date(0.0, 0.)), chunk_size(), buff_before(),
                                            processing_level(IOUtils::filtered | IOUtils::resampled | IOUtils::generated)
{
	meteoprocessor.getWindowSize(proc_properties);
	setDfltBufferProperties();
}

void TimeSeriesManager::setDfltBufferProperties()
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
}

void TimeSeriesManager::setMinBufferRequirements(const double& i_chunk_size, const double& i_buff_before)
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

void TimeSeriesManager::setProcessingLevel(const unsigned int& i_level)
{
	if (i_level >= IOUtils::num_of_levels)
		throw InvalidArgumentException("The processing level is invalid", AT);

	if (((i_level & IOUtils::raw) == IOUtils::raw)
	    && ((i_level & IOUtils::filtered) == IOUtils::filtered))
		throw InvalidArgumentException("The processing level is invalid (raw and filtered at the same time)", AT);

	processing_level = i_level;
}

void TimeSeriesManager::push_meteo_data(const IOUtils::ProcessingLevel& level, const Date& date_start, const Date& date_end,
                                const std::vector< METEO_SET >& vecMeteo)
{
	//perform check on date_start and date_end
	if (date_end < date_start) {
		std::ostringstream ss;
		ss << "Trying to push data set from " << date_start.toString(Date::ISO) << " to " << date_end.toString(Date::ISO) << ". ";
		ss << " Obviously, date_start should be less than date_end!";
		throw InvalidArgumentException(ss.str(), AT);
	}

	if (level == IOUtils::filtered) {
		fcache_start   = date_start;
		fcache_end     = date_end;
		filtered_cache = vecMeteo;
	} else if (level == IOUtils::raw) {
		fcache_start = fcache_end = Date(0.0, 0.);
		filtered_cache.clear();
		raw_start     = date_start;
		raw_end       = date_end;
		raw_buffer = vecMeteo;
	} else {
		throw InvalidArgumentException("The processing level is invalid (should be raw OR filtered)", AT);
	}

	point_cache.clear(); //clear point cache, so that we don't return resampled values of deprecated data
}

size_t TimeSeriesManager::getStationData(const Date& date, STATIONS_SET& vecStation)
{
	vecStation.clear();

	if (processing_level == IOUtils::raw){
		iohandler.readStationData(date, vecStation);
	} else {
		iohandler.readStationData(date, vecStation);
	}

	return vecStation.size();
}

//for an interval of data: decide whether data should be filtered or raw
size_t TimeSeriesManager::getMeteoData(const Date& dateStart, const Date& dateEnd, std::vector< METEO_SET >& vecVecMeteo)
{
	vecVecMeteo.clear();

	if (processing_level == IOUtils::raw){
		iohandler.readMeteoData(dateStart, dateEnd, vecVecMeteo);
	} else {
		const bool success = read_filtered_cache(dateStart, dateEnd, vecVecMeteo);

		if (!success){
			vector< vector<MeteoData> > tmp_meteo;
			fillRawBuffer(dateStart, dateEnd);
			getFromRawBuffer(dateStart, dateEnd, tmp_meteo);

			//now it needs to be secured that the data is actually filtered, if configured
			if ((IOUtils::filtered & processing_level) == IOUtils::filtered){
				//we don't use tmp_meteo, but calling fillRawBuffer has filled the buffer for us
				//and fill_filtered_cache will directly use raw_buffer
				//HACK: if raw_buffer can not hold all data between start and end
				//then this would not work
				fill_filtered_cache();
				read_filtered_cache(dateStart, dateEnd, vecVecMeteo);
			} else {
				vecVecMeteo = tmp_meteo;
			}
		}

		if ((IOUtils::generated & processing_level) == IOUtils::generated){
			dataGenerator.createParameters(vecVecMeteo);
			dataGenerator.fillMissing(vecVecMeteo);
		}
	}

	return vecVecMeteo.size(); //equivalent with the number of stations that have data
}

size_t TimeSeriesManager::getMeteoData(const Date& i_date, METEO_SET& vecMeteo)
{
	vecMeteo.clear();
	vector< vector<MeteoData> > vec_cache;

	//1. Check whether user wants raw data or processed data
	//The first case: we are looking at raw data directly, only unresampled values are considered, exact date match
	if (processing_level == IOUtils::raw) {
		iohandler.readMeteoData(i_date-Duration(1./(24.*3600.), 0.), i_date+Duration(1./(24.*3600.), 0.), vec_cache);
		for (size_t ii=0; ii<vec_cache.size(); ii++){ //for every station
			const size_t index = IOUtils::seek(i_date, vec_cache[ii], true);
			if (index != IOUtils::npos)
				vecMeteo.push_back(vec_cache[ii][index]); //Insert station into vecMeteo
		}
		return vecMeteo.size();
	}

	//2.  Check which data point is available, buffered locally
	const map<Date, vector<MeteoData> >::const_iterator it = point_cache.find(i_date);
	if (it != point_cache.end()){
		vecMeteo = it->second;
		return vecMeteo.size();
	}

	//Let's make sure we have the data we need, in the filtered_cache or in vec_cache
	const Date buffer_start( i_date-proc_properties.time_before ), buffer_end( i_date+proc_properties.time_after );
	vector< vector<MeteoData> >* data = NULL; //reference to either filtered_cache or vec_cache
	if ((IOUtils::filtered & processing_level) == IOUtils::filtered){
		const bool cached = (fcache_start <= buffer_start) && (fcache_end >= buffer_end);
		if (!cached) {
			//explicit caching, rebuffer if necessary
			fillRawBuffer(buffer_start, buffer_end);
			fill_filtered_cache();
		}
		data = &filtered_cache;
	} else { //data to be resampled should be IOUtils::raw
		fillRawBuffer(buffer_start, buffer_end);
		getFromRawBuffer(buffer_start, buffer_end, vec_cache);
		data = &vec_cache;
	}

	if ((IOUtils::resampled & processing_level) != IOUtils::resampled) { //no resampling required
		for (size_t ii=0; ii<(*data).size(); ii++) { //for every station
			const size_t index = IOUtils::seek(i_date, (*data)[ii], true); //needs to be an exact match
			if (index != IOUtils::npos)
				vecMeteo.push_back((*data)[ii][index]); //Insert station into vecMeteo
		}
	} else { //resampling required
		MeteoData md;
		for (size_t ii=0; ii<(*data).size(); ii++) { //for every station
			const bool success = meteoprocessor.resample(i_date, (*data)[ii], md);
			if (success) vecMeteo.push_back(md);
		}
	}

	if ((IOUtils::generated & processing_level) == IOUtils::generated) {
		dataGenerator.createParameters(vecMeteo);
		dataGenerator.fillMissing(vecMeteo);
	}

	add_to_points_cache(i_date, vecMeteo); //Store result in the local cache

	return vecMeteo.size();
}

void TimeSeriesManager::writeMeteoData(const std::vector< METEO_SET >& vecMeteo, const std::string& name)
{
	if (processing_level == IOUtils::raw){
		iohandler.writeMeteoData(vecMeteo, name);
	} else {
		iohandler.writeMeteoData(vecMeteo, name);
	}
}

double TimeSeriesManager::getAvgSamplingRate() const
{
	if (processing_level == IOUtils::raw || raw_buffer.empty())
		return IOUtils::nodata;

	const size_t nr_stations = raw_buffer.size();
	double sum = 0;
	for (size_t ii=0; ii<nr_stations; ii++){ //loop over all stations
		if(!raw_buffer[ii].empty()) {
			const std::vector<MeteoData>& curr_station = raw_buffer[ii];
			const double days = curr_station.back().date.getJulian() - curr_station.front().date.getJulian();

			//add the average sampling rate for this station
			const size_t nr_data_pts = raw_buffer[ii].size();
			if(days>0.) sum += (double)(nr_data_pts-1) / days; //the interval story: 2 points define 1 interval!
		}
	}
	if (sum > 0.){
		return ((double)sum / (double)(nr_stations*24*3600)); //in points per seconds, ie Hz
	}

	return IOUtils::nodata;
}

/**
 * @brief Filter the whole raw meteo data buffer
 */
void TimeSeriesManager::fill_filtered_cache()
{
	if ((IOUtils::filtered & processing_level) == IOUtils::filtered){
		//ask the bufferediohandler for the whole buffer
		fcache_start = raw_start;
		fcache_end = raw_end;
		const vector< METEO_SET >& buffer( raw_buffer );
		meteoprocessor.process(buffer, filtered_cache);
	}
}

/**
 * @brief Try to cut out a chunk of the time series stored in filtered_cache
 * @param start_date The start date of the chunk to be cut out (inclusive)
 * @param start_date The end date of the chunk to be cut out (inclusive)
 * @param vec_meteo  A vector to store the chunk cut out
 * @return true if the requested chunk was contained by filtered_cache, false otherwise
 */
bool TimeSeriesManager::read_filtered_cache(const Date& start_date, const Date& end_date, std::vector< METEO_SET >& vec_meteo)
{
	if ((start_date >= fcache_start) && (end_date <= fcache_end)){
		//it's already in the filtered_cache, so just copy the requested slice
		for (size_t ii=0; ii<filtered_cache.size(); ii++){ //loop over stations
			size_t startpos = IOUtils::seek(start_date, filtered_cache[ii], false);
			if (startpos == IOUtils::npos){
				if (!filtered_cache[ii].empty()){
					if (filtered_cache[ii][0].date <= end_date){
						startpos = 0;
					}
				}
			}

			if (startpos != IOUtils::npos){
				vec_meteo.push_back(vector<MeteoData>());
				for (size_t jj=startpos; jj<filtered_cache[ii].size(); jj++){
					const MeteoData& md = filtered_cache[ii][jj];
					if (md.date <= end_date){
						vec_meteo.back().push_back(md);
					} else {
						break;
					}
				}
			}
		}

		return true;
	}

	return false;
}

void TimeSeriesManager::add_to_points_cache(const Date& i_date, const METEO_SET& vecMeteo)
{
	//Check cache size, delete oldest elements if necessary
	if (point_cache.size() > 2000) {
		point_cache.clear(); //HACK: implement a true ring buffer!
	}

	point_cache[i_date] = vecMeteo;
}

void TimeSeriesManager::clear_cache()
{
	raw_buffer.clear();
	raw_start = Date(0., 0.);
	raw_end = Date(0., 0.);

	filtered_cache.clear();
	fcache_start = Date(0., 0.);
	fcache_end = Date(0., 0.);

	point_cache.clear();
}

/**
 * @brief return all the buffered data between the given dates
 * @param date_start requested start date of the buffer
 * @param date_end requested end date of the buffer
 * @param data vector to fill with the buffered data
 */
void TimeSeriesManager::getFromRawBuffer(const Date& date_start, const Date& date_end, std::vector< METEO_SET > &vecMeteo)
{
	//1. Prepare the output vector
	const size_t buffer_size = raw_buffer.size();
	vecMeteo.clear();
	vecMeteo.reserve(buffer_size);

	//2. Copy appropriate data into vecMeteo for each station
	for (size_t ii=0; ii<buffer_size; ii++){ //loop through stations
		vecMeteo.push_back(vector<MeteoData>()); //insert one empty vector of MeteoData

		if (raw_buffer[ii].empty()) continue; //no data in buffer for this station

		size_t pos_start = IOUtils::seek(date_start, raw_buffer[ii], false);
		if (pos_start == IOUtils::npos) pos_start = 0;

		size_t pos_end = IOUtils::seek(date_end, raw_buffer[ii], false);//HACK:: edit IOUtils::seek to accept an offset
		if (pos_end == IOUtils::npos) pos_end = raw_buffer[ii].size() - 1; //just copy until the end of the buffer

		if (raw_buffer[ii][pos_end].date > date_end){
			if (pos_end > pos_start) pos_end--;
		} else {
			pos_end++;
		}
		vecMeteo[ii].reserve(pos_end-pos_start+1); //weird that the "insert" does not handle it internally...
		vecMeteo[ii].insert(vecMeteo[ii].begin(), raw_buffer[ii].begin()+pos_start, raw_buffer[ii].begin()+pos_end);
	}
}

void TimeSeriesManager::fillRawBuffer(const Date& date_start, const Date& date_end)
{
	const Date new_buffer_start(date_start-buff_before); //taking centering into account
	Date new_buffer_end(new_buffer_start + chunk_size);

	//Read MeteoData for requested interval in chunks, furthermore buffer it
	//Try to buffer after the requested chunk for subsequent calls

	//0. initialize if not already initialized
	if (raw_buffer.empty()) {
		iohandler.readMeteoData(new_buffer_start, new_buffer_end, raw_buffer);
		raw_start = new_buffer_start;
		raw_end   = new_buffer_end;
	}

	//1. Check whether data is in buffer already, and buffer it if not
	if ((date_start < raw_start) || (date_end > raw_end)) {
		//rebuffer data
		if ((new_buffer_end != raw_end) || (new_buffer_start != raw_start)) { //rebuffer for real
			raw_buffer.clear(); //the plugins do it internally anyway, but this is cheap and safe...
			iohandler.readMeteoData(new_buffer_start, new_buffer_end, raw_buffer);
			raw_start = new_buffer_start;
			raw_end   = new_buffer_end;
		}

		const size_t buffer_size = raw_buffer.size();
		vector< vector<MeteoData> > tmp_raw_buffer;
		while (date_end > new_buffer_end){
			//if the requested interval is bigger than a normal buffer, we have to increase the buffer anyway...
			tmp_raw_buffer.reserve(buffer_size);
			iohandler.readMeteoData(new_buffer_end, new_buffer_end+chunk_size, tmp_raw_buffer);

			if (tmp_raw_buffer.size() != buffer_size) {
				ostringstream ss;
				ss << "The number of stations changed over time from " << buffer_size << " to " << tmp_raw_buffer.size() << ", ";
				ss << "this is not handled yet!";
				throw IOException(ss.str(), AT);
			}

			//Loop through stations and append data
			for (size_t ii=0; ii<buffer_size; ii++){ //loop through stations
				if ((!raw_buffer[ii].empty()) && (!tmp_raw_buffer[ii].empty())){
					//check if the last element equals the first one
					if (raw_buffer[ii].back().date >= tmp_raw_buffer[ii].front().date)
						raw_buffer[ii].pop_back(); //delete the element with the same date
				}

				raw_buffer[ii].reserve(raw_buffer[ii].size()+tmp_raw_buffer[ii].size());
				raw_buffer[ii].insert(raw_buffer[ii].end(), tmp_raw_buffer[ii].begin(), tmp_raw_buffer[ii].end());
			}
			new_buffer_end += chunk_size;
			raw_end = new_buffer_end;
		}
	}
}

const std::string TimeSeriesManager::toString() const {
	ostringstream os;
	os << "<TimeSeriesManager>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";
	os << "IOHandler& iohandler = " << hex << &iohandler << dec << "\n";
	os << meteoprocessor.toString();
	os << "Processing level = " << processing_level << "\n";
	os << dataGenerator.toString();

	//display raw_buffer
	os << "RawBuffer content (" << raw_buffer.size() << " stations)\n";
	for(size_t ii=0; ii<raw_buffer.size(); ii++) {
		if (!raw_buffer[ii].empty()){
			os << std::setw(10) << raw_buffer[ii].front().meta.stationID << " = "
			   << raw_buffer[ii].front().date.toString(Date::ISO) << " - "
			   << raw_buffer[ii].back().date.toString(Date::ISO) << ", "
			   << raw_buffer[ii].size() << " timesteps" << endl;
		}
	}

	//display filtered_cache
	os << "Filteredcache content (" << filtered_cache.size() << " stations)\n";
	for(size_t ii=0; ii<filtered_cache.size(); ii++) {
		if (!filtered_cache[ii].empty()){
			os << std::setw(10) << filtered_cache[ii].front().meta.stationID << " = "
			   << filtered_cache[ii].front().date.toString(Date::ISO) << " - "
			   << filtered_cache[ii].back().date.toString(Date::ISO) << ", "
			   << filtered_cache[ii].size() << " timesteps" << endl;
		}
	}

	//display meteocache
	size_t count=0;
	size_t min_stations=std::numeric_limits<size_t>::max();
	size_t max_stations=0;
	std::map<Date, std::vector<MeteoData> >::const_iterator iter = point_cache.begin();
	for (; iter != point_cache.end(); ++iter) {
		const size_t nb_stations = iter->second.size();
		if(nb_stations>max_stations) max_stations=nb_stations;
		if(nb_stations<min_stations) min_stations=nb_stations;
		count++;
	}

	if(count==0) {
		os << "Resampled cache is empty\n";
	}
	if(count==1) {
		os << "Resampled cache content (";
		if(max_stations==min_stations)
			os << min_stations;
		else
			os << min_stations << " to " << max_stations;
		os << " station(s))\n";
		os << std::setw(22) << point_cache.begin()->first.toString(Date::ISO) << " - 1 timestep\n";
	}
	if(count>1) {
		const double avg_sampling = ( (point_cache.rbegin()->first.getJulian()) - (point_cache.begin()->first.getJulian()) ) / (double)(count-1);

		os << "Resampled cache content (";
		if(max_stations==min_stations)
			os << min_stations;
		else
			os << min_stations << " to " << max_stations;
		os << " station(s))\n";
		os << std::setw(22) << point_cache.begin()->first.toString(Date::ISO);
		os << " - " << point_cache.rbegin()->first.toString(Date::ISO);
		os << " - " << count << " timesteps (" << setprecision(3) << fixed << avg_sampling*24.*3600. << " s sampling rate)";
	}

	os << "</TimeSeriesManager>\n";
	return os.str();
}

} //namespace
