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

IOManager::IOManager(const Config& i_cfg) : cfg(i_cfg), rawio(i_cfg), bufferedio(rawio, i_cfg), meteoprocessor(i_cfg)
{
	setProcessingLevel(IOManager::filtered | IOManager::resampled);

	fcache_start = fcache_end = Date(0.0, 0.); //this should not matter, since 0 is still way back before any real data...

	meteoprocessor.getWindowSize(proc_properties);
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

double IOManager::getAvgSamplingRate()
{
	if (processing_level == IOManager::raw){
		return IOUtils::nodata;
	} else {
		return bufferedio.getAvgSamplingRate();
	}
}

unsigned int IOManager::getStationData(const Date& date, STATION_TIMESERIE& vecStation)
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
unsigned int IOManager::getMeteoData(const Date& dateStart, const Date& dateEnd, 
                                     std::vector< METEO_TIMESERIE >& vecMeteo)
{
	vecMeteo.clear();

	if (processing_level == IOManager::raw){
		rawio.readMeteoData(dateStart, dateEnd, vecMeteo);
	} else {
		bool success = read_filtered_cache(dateStart, dateEnd, vecMeteo);

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
				read_filtered_cache(dateStart, dateEnd, vecMeteo);
			} else {
				vecMeteo = tmp_meteo;
			}
		}
	}

	return vecMeteo.size(); //equivalent with the number of stations that have data
}

/**
 * @brief Filter the whole meteo data buffer provided by bufferedio
 */
void IOManager::fill_filtered_cache()
{
	if ((IOManager::filtered & processing_level) == IOManager::filtered){
		//ask the bufferediohandler for the whole buffer
		const vector< METEO_TIMESERIE >& buffer = bufferedio.get_complete_buffer(fcache_start, fcache_end);

		//cout << "Now filtering ..." << endl;
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
bool IOManager::read_filtered_cache(const Date& start_date, const Date& end_date,
                                    std::vector< METEO_TIMESERIE >& vec_meteo)
{
	if ((start_date >= fcache_start) && (end_date <= fcache_end)){
		//it's already in the filtered_cache, so just copy the requested slice
		for (unsigned int ii=0; ii<filtered_cache.size(); ii++){
			unsigned int startpos = IOUtils::seek(start_date, filtered_cache[ii], false);
			if (startpos == IOUtils::npos){
				if (filtered_cache[ii].size() > 0){
					if (filtered_cache[ii][0].date <= end_date){
						startpos = 0;
					}
				}
			}

			if (startpos != IOUtils::npos){
				vec_meteo.push_back(vector<MeteoData>());
				unsigned int index = vec_meteo.size()-1;
				for (unsigned int jj=startpos; jj<filtered_cache[ii].size(); jj++){
					const MeteoData& md = filtered_cache[ii][jj];
					if (md.date <= end_date){
						vec_meteo[index].push_back(md);
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

void IOManager::add_to_cache(const Date& i_date, const METEO_TIMESERIE& vecMeteo)
{
	//Check cache size, delete oldest elements if necessary
	if (resampled_cache.size() > 200){
		resampled_cache.clear();
		//resampled_cache.erase(resampled_cache.begin(), resampled_cache.begin()+50);
	}

	resampled_cache[i_date] = vecMeteo;
}

//data can be raw or processed (filtered, resampled)
unsigned int IOManager::getMeteoData(const Date& i_date, METEO_TIMESERIE& vecMeteo)
{
	vecMeteo.clear();

	vector< vector<MeteoData> > vec_cache;

	//1. Check whether user wants raw data or processed data
	if (processing_level == IOManager::raw){
		rawio.readMeteoData(i_date-Duration(0.001, 0.), i_date+Duration(0.001, 0.), vec_cache);
		for (unsigned int ii=0; ii<vec_cache.size(); ii++){
			unsigned int index = IOUtils::seek(i_date, vec_cache[ii], true);
			if (index != IOUtils::npos)
				vecMeteo.push_back(vec_cache[ii][index]); //Insert station into vecMeteo
		}
		
		return vecMeteo.size();
	}


	//2.  Check which data point is available, buffered locally
	map<Date, vector<MeteoData> >::const_iterator it = resampled_cache.find(i_date);
	if (it != resampled_cache.end()){
		vecMeteo = it->second;
		return vecMeteo.size();
	}

	//request an appropriate window of filtered or unfiltered data
	getMeteoData(i_date-proc_properties.time_before, i_date+proc_properties.time_after, vec_cache);
	//vec_cache is either filtered or unfiltered, in any case it is wise to buffer it

	for (unsigned int ii=0; ii<vec_cache.size(); ii++){//resampling for every station
		if ((IOManager::resampled & processing_level) == IOManager::resampled){
			//cout << "Resampling data for station " << ii << " (" << vec_cache[ii].size() << " elements)" << endl;
			unsigned int position = meteoprocessor.resample(i_date, vec_cache[ii]);
			vecMeteo.push_back(vec_cache[ii][position]);
		} else { //only filtering activated
			unsigned int index = IOUtils::seek(i_date, vec_cache[ii], true);
			if (index != IOUtils::npos)
				vecMeteo.push_back(vec_cache[ii][index]); //Insert station into vecMeteo
		}
	}

	//Store result in the local cache
	add_to_cache(i_date, vecMeteo);

	return vecMeteo.size();
}
#ifdef _POPC_ //HACK popc
void IOManager::writeMeteoData(/*const*/ std::vector< METEO_TIMESERIE >& vecMeteo, /*const*/ std::string& name)
#else
void IOManager::writeMeteoData(const std::vector< METEO_TIMESERIE >& vecMeteo, const std::string& name)
#endif
{
	if (processing_level == IOManager::raw){
		rawio.writeMeteoData(vecMeteo, name);
	} else {
		bufferedio.writeMeteoData(vecMeteo, name);
	}
}

#ifdef _POPC_ //HACK popc
void IOManager::interpolate(/*const*/ Date& date, /*const*/ DEMObject& dem, /*const*/ MeteoData::Parameters& meteoparam,
                            Grid2DObject& result)
#else
void IOManager::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam, 
                            Grid2DObject& result)
#endif
{
	string info_string;
	interpolate(date, dem, meteoparam, result, info_string);
}

#ifdef _POPC_ //HACK popc
void IOManager::interpolate(/*const*/ Date& date, /*const*/ DEMObject& dem, /*const*/ MeteoData::Parameters& meteoparam,
                            Grid2DObject& result, std::string& info_string)
#else
void IOManager::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                            Grid2DObject& result, std::string& info_string)
#endif
{
	Meteo2DInterpolator mi(cfg, *this);
	mi.interpolate(date, dem, meteoparam, result, info_string);
}

void IOManager::read2DGrid(Grid2DObject& grid2D, const std::string& filename)
{
	if (processing_level == IOManager::raw){
		rawio.read2DGrid(grid2D, filename);
	} else {
		bufferedio.read2DGrid(grid2D, filename);
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

void IOManager::readSpecialPoints(std::vector<Coords>& cpa)
{
	if (processing_level == IOManager::raw){
		rawio.readSpecialPoints(cpa);
	} else {
		bufferedio.readSpecialPoints(cpa);
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

std::string IOManager::toString() const {
	stringstream os;

	os << "<IOManager>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";
#ifndef _POPC_ //HACK popc
	os << rawio;
#endif
	os << bufferedio;
	os << meteoprocessor;
	os << "Processing level = " << processing_level << "\n";

	//display meteocache
	unsigned int count=0;
	unsigned int min_stations=std::numeric_limits<unsigned int>::max();
	unsigned int max_stations=-std::numeric_limits<unsigned int>::max();
	std::map<Date, std::vector<MeteoData> >::const_iterator iter = resampled_cache.begin();
	for (; iter != resampled_cache.end(); iter++) {
		const unsigned int nb_stations = iter->second.size();
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
		os << resampled_cache.begin()->first.toString(Date::ISO) << " - 1 timestep\n";
	}
	if(count>1) {
		const double avg_sampling = ( (resampled_cache.rbegin()->first.getJulianDate()) - (resampled_cache.begin()->first.getJulianDate()) ) / (double)(count-1);

		os << "Resampled cache content (";
		if(max_stations==min_stations)
			os << min_stations;
		else
			os << min_stations << " to " << max_stations;
		os << " station(s))\n";
		os << resampled_cache.begin()->first.toString(Date::ISO);
		os << " - " << resampled_cache.rbegin()->first.toString(Date::ISO);
		os << " - " << count << " timesteps (" << setprecision(3) << fixed << avg_sampling*24.*3600. << " s sampling rate)";
	}

	//display filtered_cache
	os << "Filteredcache content (" << filtered_cache.size() << " stations)\n";
	for(unsigned int ii=0; ii<filtered_cache.size(); ii++) {
		if (filtered_cache[ii].size() > 0){
			os << std::setw(10) << filtered_cache[ii][0].meta.stationID << " = "
			   << filtered_cache[ii][0].date.toString(Date::ISO) << " - "
			   << filtered_cache[ii][filtered_cache[ii].size()-1].date.toString(Date::ISO) << ", "
			   << filtered_cache[ii].size() << " timesteps" << endl;
		}
	}

	os << "</IOManager>\n";
	return os.str();
}

//#ifndef _POPC_
std::ostream& operator<<(std::ostream& os, const IOManager& io)
{
	os << io.toString();
	return os;
}
//#endif

} //namespace
