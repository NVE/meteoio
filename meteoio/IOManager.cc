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

	if (processing_level == IOManager::raw){
		rawio.readMeteoData(dateStart, dateEnd, vecMeteo);
	} else {
		bufferedio.readMeteoData(dateStart, dateEnd, vecMeteo);
		//now it needs to secured that the data is actually filtered, if configured
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

	Date time_before, time_after;
	unsigned int npoints;
	meteoprocessor.getWindowSize(time_before, time_after, npoints);

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
	bufferedio.readMeteoData(i_date-time_before, i_date+time_after, vec_cache);
	for (unsigned int ii=0; ii<vec_cache.size(); ii++){
		MeteoData tmpmd;
		meteoprocessor.processData(i_date, vec_cache[ii], tmpmd);
		vecMeteo.push_back(tmpmd);
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

} //namespace
