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
	processing_level = IOManager::filtered | IOManager::resampled;
}

void IOManager::setProcessingLevel(const unsigned int& i_level)
{
	if (i_level >= IOManager::num_of_levels)
		throw InvalidArgumentException("The processing level is invalid", AT);

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

unsigned int IOManager::getStationData(const Date& date, std::vector<StationData>& vecStation)
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
                                     std::vector< std::vector<MeteoData> >& vecMeteo)
{
	vecMeteo.clear();

	vector< vector<MeteoData> > tmp_meteo;

	if (processing_level == IOManager::raw){
		rawio.readMeteoData(dateStart, dateEnd, vecMeteo);
	} else {
		bufferedio.readMeteoData(dateStart, dateEnd, tmp_meteo);
		//now it needs to secured that the data is actually filtered, if configured

		if ((IOManager::filtered & processing_level) == IOManager::filtered){
			//cout << "Now filtering ..." << endl;
			meteoprocessor.process(tmp_meteo, vecMeteo);
		} else {
			vecMeteo = tmp_meteo;
		}
	}

	return vecMeteo.size(); //equivalent with the number of stations that have data
}

void IOManager::add_to_cache(const Date& i_date, const std::vector<MeteoData>& vecMeteo)
{
	//Check cache size, delete oldest elements if necessary
	if (meteo_cache.size() > 200){
		meteo_cache.clear();
		//meteo_cache.erase(meteo_cache.begin(), meteo_cache.begin()+50);
	}

	meteo_cache[i_date] = vecMeteo;
}

//data can be raw or processed (filtered, resampled)
unsigned int IOManager::getMeteoData(const Date& i_date, std::vector<MeteoData>& vecMeteo)
{
	vecMeteo.clear();

	vector< vector<MeteoData> > vec_cache;

	ProcessingProperties properties;
	meteoprocessor.getWindowSize(properties);

	//1. Check whether user wants raw data or processed data
	if (processing_level == IOManager::raw){
		rawio.readMeteoData(i_date-Date(0.001), i_date+Date(0.001), vec_cache);
		for (unsigned int ii=0; ii<vec_cache.size(); ii++){
			unsigned int index = IOUtils::seek(i_date, vec_cache[ii], true);
			if (index != IOUtils::npos)
				vecMeteo.push_back(vec_cache[ii][index]); //Insert station into vecMeteo
		}
		
		return vecMeteo.size();
	}


	//2.  Check which data is available, buffered locally
	map<Date, vector<MeteoData> >::const_iterator it = meteo_cache.find(i_date);
	if (it != meteo_cache.end()){
		vecMeteo = it->second;
		return vecMeteo.size();
	}

	//    request an appropriate window of data from bufferedio
	//    Hand window of data over to meteo processor
	bufferedio.readMeteoData(i_date-properties.time_before, i_date+properties.time_after, vec_cache);
	//vec_cache is either filtered or unfiltered

	for (unsigned int ii=0; ii<vec_cache.size(); ii++){//resampling for every station
		cout << "Resampling data for station " << ii << " (" << vec_cache[ii].size() << " elements)" << endl;
		unsigned int position = meteoprocessor.resample(i_date, vec_cache[ii]);
		vecMeteo.push_back(vec_cache[ii][position]);
	}

	//Store result in the local cache
	add_to_cache(i_date, vecMeteo);

	return vecMeteo.size();
}

void IOManager::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string& name)
{
	if (processing_level == IOManager::raw){
		rawio.writeMeteoData(vecMeteo, name);
	} else {
		bufferedio.writeMeteoData(vecMeteo, name);
	}
}

void IOManager::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam, 
                            Grid2DObject& result)
{
	string info_string;
	interpolate(date, dem, meteoparam, result, info_string);
}
	
void IOManager::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
					   Grid2DObject& result, std::string& info_string)
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

std::ostream& operator<<(std::ostream& os, const IOManager& io)
{
	os << "<IOManager>\n";
	os << "Config cfg = " << hex << &io.cfg << "\n";
	os << io.rawio;
	os << io.bufferedio;
	os << io.meteoprocessor;
	os << "Processing level = " << io.processing_level << "\n";

	unsigned int count=0;
	unsigned int min_stations=std::numeric_limits<unsigned int>::max();
	unsigned int max_stations=-std::numeric_limits<unsigned int>::max();
	std::map<Date, std::vector<MeteoData> >::const_iterator iter = io.meteo_cache.begin();
	for (; iter != io.meteo_cache.end(); iter++) {
		const unsigned int nb_stations = iter->second.size();
		if(nb_stations>max_stations) max_stations=nb_stations;
		if(nb_stations<min_stations) min_stations=nb_stations;
		count++;
	}

	if(count==0) {
		os << "Meteo cache is empty\n";
	}
	if(count==1) {
		os << "Meteo cache contains 1 element at " << io.meteo_cache.begin()->first.toString(Date::ISO);
		os << " for ";
		if(max_stations==min_stations)
			os << min_stations << " station(s)\n";
		else
			os << "between " << min_stations << " and " << max_stations << " stations\n";
	}
	if(count>1) {
		const double avg_sampling = ( (io.meteo_cache.rbegin()->first.getJulianDate()) - (io.meteo_cache.begin()->first.getJulianDate()) ) / (double)(count-1);
		os << "Meteo cache goes from " << io.meteo_cache.begin()->first.toString(Date::ISO);
		os << " to " << io.meteo_cache.rbegin()->first.toString(Date::ISO);
		os << " with " << count << " timesteps (" << setprecision(3) << fixed << avg_sampling*24.*3600. << " s sampling rate)";
		os << " for ";
		if(max_stations==min_stations)
			os << min_stations << " station(s)\n";
		else
			os << "between " << min_stations << " and " << max_stations << " stations\n";
	}

	os << "</IOManager>\n";
	return os;
}

} //namespace
