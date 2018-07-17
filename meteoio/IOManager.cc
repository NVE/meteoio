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
/**
 * @page spatial_resampling Spatial resampling handling
 * It is possible to use spatially interpolated meteorological fields or time series of 2D grids to extract meteorological time series for a set of points.
 * This is handled as "spatial resampling" and the data **will seem to originate from these virtual stations points** where no station is present. This obviously
 * comes at the cost of much higher run times. Several strategies are available (with the *RESAMPLING_STRATEGY* keyword, after setting RESAMPLING to TRUE):
 *     + VSTATIONS: points measurements are spatially interpolated at the chosen locations;
 *     + GRID_EXTRACT: gridded values are extracted for the cells containing the given locations;
 *     + GRID_SMART: the four nodes surrounding each given locations are extracted and potential duplicates are removed (so it is ready for performing spatial interpolations);
 *     + GRID_ALL: all grid points are extracted.
 * 
 * Currently, it is necessary to provide a hint on how often the data should be extrated versus temporally interpolated between extracted point. This is described
 * by providing a refresh rate and an offset (both in seconds, with the VSTATIONS_REFRESH_RATE and VSTATIONS_REFRESH_OFFSET keywords, respectively)
 * \image html vstations_sampling.png "Resampling workflow"
 * \image latex vstations_sampling.eps "Resampling workflow" width=0.9\textwidth
 *
 * @section vstations VSTATIONS
 * The data from real input stations (as read by the plugin defined with the METEO key in the [input] section) is filtered/processed, temporally interpolated and 
 * spatially interpolated as defined in the configuration file. Then time series are reconstructed from these grids at a set of defined points (which will receive
 * station IDs such as <i>VIR#</i> for each station). This behavior is configured by the following keys (in the [Input] section):
 *    + RESAMPLING set to *true*;
 *    + RESAMPLING_STRATEGY set to VSTATIONS;
 *    + VSTATION# : provide the lat, lon and altitude or easting, northing and altitude for a virtual station (see \link Coords::Coords(const std::string& in_coordinatesystem, const std::string& in_parameters, std::string coord_spec) Coords()\endlink for the syntax);
 *    + VIRTUAL_PARAMETERS: list of MeteoData::Parameters that have to be interpolated to populate the virtual stations;
 *    + VSTATIONS_REFRESH_RATE: how often to rebuild the spatial interpolations, in seconds;
 *    + VSTATIONS_REFRESH_OFFSET: time offset to the stations' refresh rate, in seconds;
 *    + INTERPOL_USE_FULL_DEM: should the spatial interpolations be performed on the whole DEM? (this is necessary for some algorithms, for example WINSTAL);
 * 
 * Currently, a DEM also has to be provided since this will be used to retrieve the elevation, slope and azimuth of the virtual stations.
 * 
 * In the example provided below, 4 stations provide the original data that will be spatially interpolated at 2 points (or virtual stations, VIR1 and VIR2) for
 * 7 meteorological parameters. Every 6 hours, with starting offset of on hour, the original data will be spatially interpolated (so at 01:00, 07:00, 13:00 and 19:00).
 * Any data requested at other time steps will be temporally resampled from the spatially interpolated data.
 * @code
 * DEM = ARC
 * DEMFILE = ./input/surface-grids/davos.asc
 *
 * #here, the real data as measured by some stations 
 * METEO		= IMIS
 * DBNAME		= sdbo
 * DBUSER		= xxx
 * DBPASS		= xxx
 * STATION1	= *WFJ
 * STATION2	= STB2
 * STATION3	= WFJ2
 * STATION4	= *DAV
 * 
 * #here the locations where the data will be generated. The caller will only see these stations!
 * RESAMPLING = true
 * RESAMPLING_STRATEGY = VSTATIONS
 * VSTATION1 = latlon 46.793029 9.821343 1987
 * VSTATION2 = latlon 46.793031 9.831572 2300
 * Virtual_parameters = TA RH PSUM ILWR P VW RSWR
 * VSTATIONS_REFRESH_RATE = 21600
 * VSTATIONS_REFRESH_OFFSET = 3600
 * @endcode
 * 
 * @section grids_extract From gridded data
 * The meteorological time series are extracted from time series of user-provided grids. therefore a plugin for 2D grids must have been defined (with the GRID2D key in
 * the [Input] section). The following keys control this downscaling process:
 *    + RESAMPLING set to *true*;
 *    + RESAMPLING_STRATEGY set to either *GRID_EXTRACT* or *GRID_ALL*;
 *    + VSTATION# : provide the lat, lon and (optionally) the epsg code for a virtual station;
 *    + VIRTUAL_PARAMETERS: list of MeteoData::Parameters that have to be interpolated to populate the virtual stations;
 *
 * Currently, a DEM has to be provided in order to check the position of the stations and the consistency of the grids.
 * @code
 * DEM     = NETCDF
 * DEMFILE = ./input/grids/era-interim-dem.nc
 * 
 * GRID2D    = NETCDF
 * GRID2DFILE =  ./input/grids/era-interim.nc
 * NETCDF_SCHEMA = ECMWF
 * 
 * RESAMPLING = true
 * RESAMPLING_STRATEGY = GRID_EXTRACT
 * Virtual_parameters = TA RH PSUM ILWR P VW ISWR
 * 
 * #here the locations where the data will be generated. The caller will only see these stations!
 * VSTATION1 = latlon 43.359188 6.693612 150 ;great station
 * VSTATION2 = latlon 43.324887 6.629711 57 ;another great station
 * @endcode
 * @note If the temporal resolution of the gridded data is greater than the WINDOW_SIZE of the [Interpolations1D], then no data will be interpolated.
 * 
 * @section resampling_behind_the_scene Behind the scene
 * Behind the scene, this is a two stages setup: the IOManager uses either a TimeSeriesManager or a GridsManager object to retrieve the real data
 * and then injects this data as raw data into another TimeSeriesManager (so the temporal interpolations can be computed, the data generators can be called, etc).
 * Then the IOManager request the data as usual from the final TimeSeriesManager.
 *
 */

IOHandler::operation_mode IOManager::setMode(const Config& i_cfg)
{
	const std::string resampling_strategy_str = i_cfg.get("Resampling_strategy", "Input", IOUtils::nothrow);
	if (resampling_strategy_str.empty())
		return IOHandler::STD;
	if (resampling_strategy_str=="VSTATION")
		return IOHandler::VSTATION;
	if (resampling_strategy_str=="GRID_RESAMPLE")
		return IOHandler::GRID_RESAMPLE;
	if (resampling_strategy_str=="GRID_EXTRACT")
		return IOHandler::GRID_EXTRACT;
	if (resampling_strategy_str=="GRID_ALL")
		return IOHandler::GRID_ALL;
	if (resampling_strategy_str=="GRID_SMART")
		return IOHandler::GRID_SMART;
	//TODO TRAJECTORY
	
	throw InvalidArgumentException("The selected resampling_strategy is not supported", AT);
}

//TODO write an IOHandler that can directly tap into the buffers of a tsm or gdm
IOManager::IOManager(const std::string& filename_in) : cfg(filename_in), iohandler(cfg),
                                                       tsm1(iohandler, cfg), tsm2(iohandler, cfg), gdm1(iohandler, cfg), interpolator(cfg, tsm1, gdm1), resampling_dem(),
                                                       v_params(), v_stations(), vstations_refresh_rate(3600), vstations_refresh_offset(0), mode(setMode(cfg))
{
	initIOManager();
}

IOManager::IOManager(const Config& i_cfg) : cfg(i_cfg), iohandler(cfg),
                                            tsm1(iohandler, cfg), tsm2(iohandler, cfg), gdm1(iohandler, cfg), interpolator(cfg, tsm1, gdm1), resampling_dem(),
                                            v_params(), v_stations(), vstations_refresh_rate(3600), vstations_refresh_offset(0), mode(setMode(cfg))
{
	initIOManager();
}

void IOManager::initIOManager()
{
	//TODO support extra parameters by getting the param index from vecTrueMeteo[0]
	if (mode>IOHandler::GRID_EXTRACT) {
		std::vector<std::string> vecStr;
		cfg.getValue("Virtual_parameters", "Input", vecStr);
		for (size_t ii=0; ii<vecStr.size(); ii++) {
			const size_t param_idx = MeteoGrids::getParameterIndex( vecStr[ii] );
			if (param_idx==IOUtils::npos)
				throw InvalidArgumentException("Invalid parameter '" + vecStr[ii] + "', only standard parameters can be extracted from grids for virtual stations! ", AT);
			v_params.push_back( param_idx );
		}
	}
	
	if (mode==IOHandler::VSTATION) {
		std::vector<std::string> vecStr;
		cfg.getValue("Virtual_parameters", "Input", vecStr);
		for (size_t ii=0; ii<vecStr.size(); ii++) {
			const size_t param_idx = MeteoData::getStaticParameterIndex( vecStr[ii] );
			if (param_idx==IOUtils::npos)
				throw InvalidArgumentException("Invalid parameter '" + vecStr[ii] + "', only standard parameters can be extracted from grids for virtual stations! ", AT);
			v_params.push_back( param_idx );
		}
	}
	
	if (mode==IOHandler::VSTATION || mode==IOHandler::GRID_RESAMPLE) { //in these cases, we do not want to re-apply the filters
		tsm2.setProcessingLevel(IOUtils::resampled | IOUtils::generated);
	}
	//TODO: force filter pass2 on tsm2?
}

void IOManager::setProcessingLevel(const unsigned int& i_level)
{
	tsm1.setProcessingLevel(i_level);
	gdm1.setProcessingLevel(i_level);
}

void IOManager::clear_cache()
{
	tsm1.clear_cache();
	tsm2.clear_cache();
	gdm1.clear_cache();
}

void IOManager::initVirtualStations()
{
	gdm1.readDEM(resampling_dem);
	
	if (mode==IOHandler::GRID_ALL || mode==IOHandler::GRID_RESAMPLE) {
		v_stations = gdm1.initVirtualStationsAtAllGridPoints(resampling_dem);
	} else if (mode==IOHandler::GRID_EXTRACT || mode==IOHandler::GRID_SMART) {
		const bool adjust_coordinates = (mode==IOHandler::GRID_EXTRACT);
		const bool fourNeighbors = (mode==IOHandler::GRID_SMART);
		v_stations = gdm1.initVirtualStations(resampling_dem, adjust_coordinates, fourNeighbors);
	} else if (mode==IOHandler::VSTATION) {
		cfg.getValue("VSTATIONS_REFRESH_RATE", "Input", vstations_refresh_rate, IOUtils::nothrow);
		cfg.getValue("VSTATIONS_REFRESH_OFFSET", "Input", vstations_refresh_offset, IOUtils::nothrow);
		v_stations = gdm1.initVirtualStations(resampling_dem, false, false);
	} else
		throw InvalidArgumentException("Unsupported mode of operation", AT);
}

size_t IOManager::getStationData(const Date& date, STATIONS_SET& vecStation)
{
	vecStation.clear();

	if (mode==IOHandler::STD) return tsm1.getStationData(date, vecStation);
	
	if (v_stations.empty()) initVirtualStations();
	vecStation = v_stations;
	return vecStation.size();
}

size_t IOManager::getMeteoData(const Date& dateStart, const Date& dateEnd, std::vector< METEO_SET >& vecVecMeteo) 
{
	if (mode==IOHandler::STD) return tsm1.getMeteoData(dateStart, dateEnd, vecVecMeteo);
	
	if (mode==IOHandler::VSTATION) {
		if (v_stations.empty()) initVirtualStations();
		const Date bufferStart( tsm2.getRawBufferStart() );
		const Date bufferEnd( tsm2.getRawBufferEnd() );
		
		if (bufferStart.isUndef() || dateStart<bufferStart || dateEnd>bufferEnd) { //TODO: smarter rebuffer! (ie partial)
			tsm2.push_meteo_data(IOUtils::raw, dateStart, dateEnd, getVirtualStationsData(resampling_dem, dateStart, dateEnd));
		}
		
		return tsm2.getMeteoData(dateStart, dateEnd, vecVecMeteo);
	}
	
	if (mode>=IOHandler::GRID_EXTRACT) {
		if (v_stations.empty()) initVirtualStations();
		const Date bufferStart( tsm1.getRawBufferStart() );
		const Date bufferEnd( tsm1.getRawBufferEnd() );
		
		if (bufferStart.isUndef() || dateStart<bufferStart || dateEnd>bufferEnd) { //TODO: smarter rebuffer! (ie partial)
			vecVecMeteo = gdm1.getVirtualStationsFromGrid(resampling_dem, v_params, v_stations, dateStart, dateEnd);
			tsm1.push_meteo_data(IOUtils::raw, dateStart, dateEnd, vecVecMeteo);
		}
		
		return tsm1.getMeteoData(dateStart, dateEnd, vecVecMeteo);
	}
	//HACK GRID_RESAMPLE is not handled!
	
	throw InvalidArgumentException("Unsuppported operation_mode", AT);
}

//data can be raw or processed (filtered, resampled)
size_t IOManager::getMeteoData(const Date& i_date, METEO_SET& vecMeteo)
{
	vecMeteo.clear();

	if (mode==IOHandler::STD) {
		return tsm1.getMeteoData(i_date, vecMeteo);
	}
	
	if (mode==IOHandler::VSTATION) {
		if (v_stations.empty()) initVirtualStations();
		const Date bufferStart( tsm2.getRawBufferStart() );
		const Date bufferEnd( tsm2.getRawBufferEnd() );
		
		if (bufferStart.isUndef() || i_date<bufferStart || i_date>bufferEnd) { //TODO: smarter rebuffer! (ie partial)
			Duration buffer_size, buffer_before;
			tsm2.getBufferProperties(buffer_size, buffer_before);
		
			const Date dateStart( i_date - buffer_before );
			const Date dateEnd( i_date - buffer_before + buffer_size );
			tsm2.push_meteo_data(IOUtils::raw, dateStart, dateEnd, getVirtualStationsData(resampling_dem, dateStart, dateEnd));
		}
		
		return tsm2.getMeteoData(i_date, vecMeteo);
	}
	
	if (mode>=IOHandler::GRID_EXTRACT) {
		if (v_stations.empty()) initVirtualStations();
		const Date bufferStart( tsm1.getRawBufferStart() );
		const Date bufferEnd( tsm1.getRawBufferEnd() );
		
		if (bufferStart.isUndef() || i_date<bufferStart || i_date>bufferEnd) { //TODO: smarter rebuffer! (ie partial)
			Duration buffer_size, buffer_before;
			tsm1.getBufferProperties(buffer_size, buffer_before);
			
			const Date dateStart( i_date - buffer_before );
			const Date dateEnd( i_date - buffer_before + buffer_size );
			tsm1.push_meteo_data(IOUtils::raw, dateStart, dateEnd, gdm1.getVirtualStationsFromGrid(resampling_dem, v_params, v_stations, dateStart, dateEnd));
		}
		
		return tsm1.getMeteoData(i_date, vecMeteo);
	}
	//HACK GRID_RESAMPLE is not handled!
	
	throw InvalidArgumentException("Unsuppported operation_mode", AT);
}

bool IOManager::getMeteoData(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                  Grid2DObject& result)
{
	std::string info_string;
	const bool status = getMeteoData(date, dem, meteoparam, result, info_string);
	cerr << "[i] Interpolating " << MeteoData::getParameterName(meteoparam);
	cerr << " (" << info_string << ") " << endl;
	return status;
}

bool IOManager::getMeteoData(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                  Grid2DObject& result, std::string& info_string)
{
	interpolator.interpolate(date, dem, meteoparam, result, info_string);
	return (!result.empty());
}

void IOManager::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                            const std::vector<Coords>& in_coords, std::vector<double>& result)
{
	std::string info_string;
	interpolate(date, dem, meteoparam, in_coords, result, info_string);
	cerr << "[i] Interpolating " << MeteoData::getParameterName(meteoparam);
	cerr << " (" << info_string << ") " << endl;
}

void IOManager::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                            const std::vector<Coords>& in_coords, std::vector<double>& result, std::string& info_string)
{
	interpolator.interpolate(date, dem, meteoparam, in_coords, result, info_string);
}

void IOManager::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                            const std::vector<StationData>& in_stations, std::vector<double>& result, std::string& info_string)
{
	interpolator.interpolate(date, dem, meteoparam, in_stations, result, info_string);
}

double IOManager::getAvgSamplingRate() const 
{
	return tsm1.getAvgSamplingRate();
}

void IOManager::add_to_points_cache(const Date& i_date, const METEO_SET& vecMeteo) 
{
	tsm1.add_to_points_cache(i_date, vecMeteo);
}

std::vector<METEO_SET> IOManager::getVirtualStationsData(const DEMObject& dem, const Date& dateStart, const Date& dateEnd)
{
	const Date buff_start( Date::rnd(dateStart-vstations_refresh_offset, vstations_refresh_rate, Date::DOWN) + vstations_refresh_offset/(24.*3600.) );
	const Date buff_end( Date::rnd(dateEnd-vstations_refresh_offset, vstations_refresh_rate, Date::UP) + vstations_refresh_offset/(24.*3600.) );
	std::string info_string;
	
	std::vector<METEO_SET> vecvecMeteo(v_stations.size());
	const double date_inc = static_cast<double>(vstations_refresh_rate) / (24.*3600.);
	for (Date date=buff_start; date<=buff_end; date += date_inc) {
		for (size_t ii=0; ii<v_stations.size(); ii++) {
			MeteoData md(date, v_stations[ii]);
			vecvecMeteo[ii].push_back( md );
		}
		
		for (size_t param=0; param<v_params.size(); param++) {
			std::vector<double> result;
			interpolate(date, dem, static_cast<MeteoData::Parameters>(v_params[param]), v_stations, result, info_string);
			for (size_t ii=0; ii<v_stations.size(); ii++) {
				vecvecMeteo[ii].back()(v_params[param]) = result[ii];
			}
		}

	}
	
	return vecvecMeteo;
}

const std::string IOManager::toString() const {
	ostringstream os;
	os << "<IOManager>\n";
	os << "Config cfg = " << hex << &cfg << dec << "\n";
	os << iohandler.toString();
	os << tsm1.toString();
	os << gdm1.toString();
	os << interpolator.toString();
	os << "</IOManager>\n";
	return os.str();
}

} //namespace
