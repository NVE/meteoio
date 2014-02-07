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

#include <meteoio/IOManager.h>

using namespace std;

namespace mio {

IOManager::IOManager(const Config& i_cfg) : cfg(i_cfg), rawio(cfg), bufferedio(rawio, cfg),
                                            meteoprocessor(cfg), interpolator(cfg), dataGenerator(cfg),
                                            v_params(), v_coords(), v_stations(),
                                            proc_properties(), point_cache(), filtered_cache(),
                                            fcache_start(Date(0.0, 0.)), fcache_end(Date(0.0, 0.)), //this should not matter, since 0 is still way back before any real data...
                                            processing_level(IOManager::filtered | IOManager::resampled | IOManager::generated),
                                            virtual_stations(false), skip_virtual_stations(false)
{
	meteoprocessor.getWindowSize(proc_properties);
	interpolator.setIOManager(*this); //because "*this" does not necessarily exist in the initialization list...
	cfg.getValue("Virtual_stations", "Input", virtual_stations, IOUtils::nothrow);
	if(virtual_stations) {
		initVirtualStations();
	}
}

void IOManager::initVirtualStations()
{
	if(!cfg.keyExists("DEM", "Input"))
		throw NoAvailableDataException("In order to use virtual stations, please provide a DEM!", AT);
	DEMObject dem;
	bufferedio.readDEM(dem);

	//get virtual stations coordinates
	std::string coordin, coordinparam, coordout, coordoutparam;
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);

	std::vector<std::string> vecStation;
	cfg.getValues("Vstation", "INPUT", vecStation);
	for(size_t ii=0; ii<vecStation.size(); ii++) {
		Coords tmp(coordin, coordinparam, vecStation[ii]);
		if(!tmp.isNodata())
			v_coords.push_back( tmp );
	}

	//create stations' metadata
	for(size_t ii=0; ii<v_coords.size(); ii++) {
		if(!dem.gridify(v_coords[ii])) {
			ostringstream ss;
			ss << "Virtual station \"" << vecStation[ii] << "\" is not contained is provided DEM";
			throw NoAvailableDataException(ss.str(), AT);
		}

		const size_t i = v_coords[ii].getGridI(), j = v_coords[ii].getGridJ();
		v_coords[ii].setAltitude(dem(i,j), false);

		ostringstream name;
		name << "Virtual_Station_" << ii+1;
		ostringstream id;
		id << "VIR" << ii+1;
		StationData sd(v_coords[ii], id.str(), name.str());
		sd.setSlope(dem.slope(i,j), dem.azi(i,j));

		v_stations.push_back( sd );
	}
}

void IOManager::setProcessingLevel(const unsigned int& i_level)
{
	if (i_level >= IOManager::num_of_levels)
		throw InvalidArgumentException("The processing level is invalid", AT);

	if (((i_level & IOManager::raw) == IOManager::raw)
	    && ((i_level & IOManager::filtered) == IOManager::filtered))
		throw InvalidArgumentException("The processing level is invalid (raw and filtered at the same time)", AT);

	processing_level = i_level;
}

void IOManager::setMinBufferRequirements(const double& buffer_size, const double& buff_before) {
	bufferedio.setMinBufferRequirements(buffer_size, buff_before);
}

double IOManager::getAvgSamplingRate() const
{
	if (processing_level == IOManager::raw){
		return IOUtils::nodata;
	} else {
		return bufferedio.getAvgSamplingRate();
	}
}

const Config IOManager::getConfig() const
{
	return cfg;
}

void IOManager::push_meteo_data(const ProcessingLevel& level, const Date& date_start, const Date& date_end,
                                const std::vector< METEO_SET >& vecMeteo)
{
	//perform check on date_start and date_end
	if (date_end < date_start) {
		std::ostringstream ss;
		ss << "Trying to push data set from " << date_start.toString(Date::ISO) << " to " << date_end.toString(Date::ISO) << ". ";
		ss << " Obviously, date_start should be less than date_end!";
		throw InvalidArgumentException(ss.str(), AT);
	}

	if (level == IOManager::filtered){
		fcache_start   = date_start;
		fcache_end     = date_end;
		filtered_cache = vecMeteo;
	} else if (level == IOManager::raw){
		//push data into the BufferedIOHandler
		fcache_start = fcache_end = Date(0.0, 0.);
		filtered_cache.clear();
		bufferedio.push_meteo_data(date_start, date_end, vecMeteo);
	} else {
		throw InvalidArgumentException("The processing level is invalid (should be raw OR filtered)", AT);
	}

	point_cache.clear(); //clear point cache, so that we don't return resampled values of deprecated data
}

size_t IOManager::getStationData(const Date& date, STATIONS_SET& vecStation)
{
	vecStation.clear();

	if (processing_level == IOManager::raw){
		rawio.readStationData(date, vecStation);
	} else {
		bufferedio.readStationData(date, vecStation);
	}

	return vecStation.size();
}

//for an interval of data: decide whether data should be filtered or raw
size_t IOManager::getMeteoData(const Date& dateStart, const Date& dateEnd, std::vector< METEO_SET >& vecVecMeteo)
{
	vecVecMeteo.clear();

	if (processing_level == IOManager::raw){
		rawio.readMeteoData(dateStart, dateEnd, vecVecMeteo);
	} else {
		const bool success = read_filtered_cache(dateStart, dateEnd, vecVecMeteo);

		if (!success){
			vector< vector<MeteoData> > tmp_meteo;
			bufferedio.readMeteoData(dateStart, dateEnd, tmp_meteo);

			//now it needs to be secured that the data is actually filtered, if configured
			if ((IOManager::filtered & processing_level) == IOManager::filtered){
				//we don't use tmp_meteo, but calling readMeteoData has filled the buffer for us
				//and fill_filtered_cache will directly use the BufferedIO buffer
				//HACK: if BufferedIO's buffer can not hold all data between start and end
				//then this would not work
				fill_filtered_cache();
				read_filtered_cache(dateStart, dateEnd, vecVecMeteo);
			} else {
				vecVecMeteo = tmp_meteo;
			}
		}

		if ((IOManager::generated & processing_level) == IOManager::generated){
			dataGenerator.fillMissing(vecVecMeteo);
		}
	}

	return vecVecMeteo.size(); //equivalent with the number of stations that have data
}

/**
 * @brief Filter the whole meteo data buffer provided by bufferedio
 */
void IOManager::fill_filtered_cache()
{
	if ((IOManager::filtered & processing_level) == IOManager::filtered){
		//ask the bufferediohandler for the whole buffer
		const vector< METEO_SET >& buffer = bufferedio.getFullBuffer(fcache_start, fcache_end);
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
bool IOManager::read_filtered_cache(const Date& start_date, const Date& end_date, std::vector< METEO_SET >& vec_meteo)
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

void IOManager::add_to_cache(const Date& i_date, const METEO_SET& vecMeteo)
{
	//Check cache size, delete oldest elements if necessary
	if (point_cache.size() > 2000){
		point_cache.clear();
	}

	point_cache[i_date] = vecMeteo;
}

void IOManager::clear_cache()
{
	point_cache.clear();
}

size_t IOManager::getTrueMeteoData(const Date& i_date, METEO_SET& vecMeteo)
{
	vecMeteo.clear();
	vector< vector<MeteoData> > vec_cache;

	//1. Check whether user wants raw data or processed data
	//The first case: we are looking at raw data directly, only unresampled values are considered, exact date match
	if (processing_level == IOManager::raw) {
		rawio.readMeteoData(i_date-Duration(1./(24.*3600.), 0.), i_date+Duration(1./(24.*3600.), 0.), vec_cache);
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
	vector< vector<MeteoData> >* data = NULL; //reference to either filtered_cache or vec_cache
	if ((IOManager::filtered & processing_level) == IOManager::filtered){
		const bool cached = (fcache_start <= i_date-proc_properties.time_before) && (fcache_end >= i_date+proc_properties.time_after);
		if (!cached) {
			//explicit caching, this forces the bufferediohandler to rebuffer, if necessary
			bufferedio.readMeteoData(i_date-proc_properties.time_before, i_date+proc_properties.time_after);
			fill_filtered_cache();
		}
		data = &filtered_cache;
	} else { //data to be resampled should be IOManager::raw
		bufferedio.readMeteoData(i_date-proc_properties.time_before, i_date+proc_properties.time_after, vec_cache);
		data = &vec_cache;
	}

	if ((IOManager::resampled & processing_level) != IOManager::resampled) { //no resampling required
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

	if ((IOManager::generated & processing_level) == IOManager::generated){
		dataGenerator.fillMissing(vecMeteo);
	}

	return vecMeteo.size();
}

size_t IOManager::getVirtualMeteoData(const Date& i_date, METEO_SET& vecMeteo)
{
	vecMeteo.clear();
	METEO_SET vecTrueMeteo;
	getTrueMeteoData(i_date, vecTrueMeteo);
	if(vecTrueMeteo.empty()) return 0;

	if(v_params.empty()) {
		//get parameters to interpolate if not already done
		//we need valid data in order to handle extra parameters
		std::vector<std::string> vecStr;
		cfg.getValue("Virtual_parameters", "Input", vecStr);
		for(size_t ii=0; ii<vecStr.size(); ii++) {
			v_params.push_back( vecTrueMeteo[0].getParameterIndex(vecStr[ii]) );
		}
	}

	DEMObject dem;
	bufferedio.readDEM(dem); //this is not a big deal since it will be in the buffer

	//create stations without measurements
	for(size_t ii=0; ii<v_stations.size(); ii++) {
		MeteoData md(i_date, v_stations[ii]);
		vecMeteo.push_back( md );
	}

	//fill meteo parameters
	string info_string;
	for(size_t param=0; param<v_params.size(); param++) {
		std::vector<double> result;
		interpolate(i_date, dem, static_cast<MeteoData::Parameters>(v_params[param]), v_coords, result, info_string);
		for(size_t ii=0; ii<v_coords.size(); ii++)
			vecMeteo[ii](v_params[param]) = result[ii];
	}

	return vecMeteo.size();
}

//data can be raw or processed (filtered, resampled)
size_t IOManager::getMeteoData(const Date& i_date, METEO_SET& vecMeteo)
{
	vecMeteo.clear();

	if(!virtual_stations || skip_virtual_stations)
		getTrueMeteoData(i_date, vecMeteo);
	else
		getVirtualMeteoData(i_date, vecMeteo);

	//Store result in the local cache
	add_to_cache(i_date, vecMeteo);
	return vecMeteo.size();
}

void IOManager::writeMeteoData(const std::vector< METEO_SET >& vecMeteo, const std::string& name)
{
	if (processing_level == IOManager::raw){
		rawio.writeMeteoData(vecMeteo, name);
	} else {
		bufferedio.writeMeteoData(vecMeteo, name);
	}
}

bool IOManager::getMeteoData(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                  Grid2DObject& result)
{
	string info_string;
	const bool status = getMeteoData(date, dem, meteoparam, result, info_string);
	cerr << "[i] Interpolating " << MeteoData::getParameterName(meteoparam);
	cerr << " (" << info_string << ") " << endl;
	return status;
}

bool IOManager::getMeteoData(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                  Grid2DObject& result, std::string& info_string)
{
	skip_virtual_stations = true;
	interpolator.interpolate(date, dem, meteoparam, result, info_string);
	skip_virtual_stations = false;
	return (!result.isEmpty());
}

void IOManager::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                            const std::vector<Coords>& in_coords, std::vector<double>& result)
{
	string info_string;
	interpolate(date, dem, meteoparam, in_coords, result, info_string);
	cerr << "[i] Interpolating " << MeteoData::getParameterName(meteoparam);
	cerr << " (" << info_string << ") " << endl;
}

void IOManager::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                            const std::vector<Coords>& in_coords, std::vector<double>& result, std::string& info_string)
{
	result.clear();
	skip_virtual_stations = true;

	vector<Coords> vec_coords = in_coords;

	for (size_t ii=0; ii<vec_coords.size(); ii++) {
		const bool gridify_success = dem.gridify(vec_coords[ii]);
		if (!gridify_success)
			throw InvalidArgumentException("Coordinate given to interpolate is outside of dem", AT);

		//Make new DEM with just one point, namely the one specified by vec_coord[ii]
		//Copy all other properties of the big DEM into the new one
		DEMObject one_point_dem(dem, (unsigned)vec_coords[ii].getGridI(), (unsigned)vec_coords[ii].getGridJ(), 1, 1, false);

		one_point_dem.min_altitude = dem.min_altitude;
		one_point_dem.max_altitude = dem.max_altitude;
		one_point_dem.min_slope = dem.min_slope;
		one_point_dem.max_slope = dem.max_slope;
		one_point_dem.min_curvature = dem.min_curvature;
		one_point_dem.max_curvature = dem.max_curvature;

		Grid2DObject result_grid;
		interpolator.interpolate(date, one_point_dem, meteoparam, result_grid, info_string);
		result.push_back(result_grid(0,0));
	}
	skip_virtual_stations = false;
}

void IOManager::read2DGrid(Grid2DObject& grid2D, const std::string& filename)
{
	if (processing_level == IOManager::raw){
		rawio.read2DGrid(grid2D, filename);
	} else {
		bufferedio.read2DGrid(grid2D, filename);
	}
}

void IOManager::read2DGrid(Grid2DObject& grid2D, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if (processing_level == IOManager::raw){
		rawio.read2DGrid(grid2D, parameter, date);
	} else {
		bufferedio.read2DGrid(grid2D, parameter, date);
	}
}

void IOManager::readDEM(DEMObject& grid2D)
{
	if (processing_level == IOManager::raw){
		rawio.readDEM(grid2D);
	} else {
		bufferedio.readDEM(grid2D);
	}
}

void IOManager::readLanduse(Grid2DObject& grid2D)
{
	if (processing_level == IOManager::raw){
		rawio.readLanduse(grid2D);
	} else {
		bufferedio.readLanduse(grid2D);
	}
}

void IOManager::readAssimilationData(const Date& date, Grid2DObject& grid2D)
{
	if (processing_level == IOManager::raw){
		rawio.readAssimilationData(date, grid2D);
	} else {
		bufferedio.readAssimilationData(date, grid2D);
	}
}

void IOManager::readPOI(std::vector<Coords>& cpa)
{
	if (processing_level == IOManager::raw){
		rawio.readPOI(cpa);
	} else {
		bufferedio.readPOI(cpa);
	}
}

void IOManager::write2DGrid(const Grid2DObject& grid2D, const std::string& name)
{
	if (processing_level == IOManager::raw){
		rawio.write2DGrid(grid2D, name);
	} else {
		bufferedio.write2DGrid(grid2D, name);
	}
}

void IOManager::write2DGrid(const Grid2DObject& grid2D, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if (processing_level == IOManager::raw){
		rawio.write2DGrid(grid2D, parameter, date);
	} else {
		bufferedio.write2DGrid(grid2D, parameter, date);
	}
}

const std::string IOManager::toString() const {
	ostringstream os;
	os << "<IOManager>\n";
	os << "Config cfg = " << hex << &cfg << dec << "\n";
	os << rawio.toString();
	os << bufferedio.toString();
	os << meteoprocessor.toString();
	os << "Processing level = " << processing_level << "\n";
	os << interpolator.toString();

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
		os << point_cache.begin()->first.toString(Date::ISO) << " - 1 timestep\n";
	}
	if(count>1) {
		const double avg_sampling = ( (point_cache.rbegin()->first.getJulian()) - (point_cache.begin()->first.getJulian()) ) / (double)(count-1);

		os << "Resampled cache content (";
		if(max_stations==min_stations)
			os << min_stations;
		else
			os << min_stations << " to " << max_stations;
		os << " station(s))\n";
		os << point_cache.begin()->first.toString(Date::ISO);
		os << " - " << point_cache.rbegin()->first.toString(Date::ISO);
		os << " - " << count << " timesteps (" << setprecision(3) << fixed << avg_sampling*24.*3600. << " s sampling rate)";
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

	os << "</IOManager>\n";
	return os.str();
}

} //namespace
