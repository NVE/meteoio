/***********************************************************************************/
/*  Copyright 2009-2012 WSL Institute for Snow and Avalanche Research  SLF-DAVOS   */
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
#include <meteoio/DataEditing.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/IOUtils.h>
#include <meteoio/FileUtils.h>
#include <meteoio/dataClasses/MeteoData.h> //needed for the merge strategies

#include <algorithm>
#include <fstream>

using namespace std;

namespace mio {

/**
 * @page data_editing Input Data Editing
 * Before any filters, resampling algorithms or data generators are applied, it is possible to edit the original data (in the following order):
 *     -# \ref data_move "rename certain parameters for all stations;"
 *     -# \ref data_exclusion "exclude/keep certain parameters on a per station basis;"
 *     -# \ref data_merging "merge stations together;"
 *     -# \ref data_copy "make a copy of a certain parameter under a new parameter name for all stations;"
 *     -# \ref data_creation "create certain parameters based on some parametrizations."
 *
 * @note Please note that the processing order is the following: the MOVE directives are processed first, then the EXCLUDE directives,
 * then the KEEP directives, then the MERGE directives and finally the COPY directives. The CREATE directives only come after all the raw data
 * has been edited.
 *
 * @section data_move 1. Data renaming (MOVE)
 * It is possible to rename a meteorological parameter thanks to the MOVE key. This key can take multiple source names that will be processed in the
 * order of declaration. The syntax is new_name::MOVE = {*space delimited list of original names*}. Original names that are not found in the current
 * dataset will silently be ignored, so it is safe to provide a list that contain many possible names:
 * @code
 * [InputEditing]
 * TA::MOVE = air_temp air_temperature temperature_air
 * @endcode
 * This can be used to rename non-standard parameter names into standard ones. In this example, if TA already had some values, it will keep
 * those and only points not having a value will be filled by either air_temp or air_temperature or temperature_air (the first one in
 * the list to have a value has the priority).
 *
 * @section data_exclusion 2. Data exclusion (EXCLUDE/KEEP)
 * It is possible to exclude specific parameters from given stations (on a per station basis). This is either done by using the station ID (or the '*' wildcard)
 * followed by "::exclude" as key with a space delimited list of \ref meteoparam "meteorological parameters" to exclude for the station as key.
 * Another possibility is to provide a file containing one station ID per line followed by a space delimited list of \ref meteoparam "meteorological parameters"
 * to exclude for the station (the path to the file can be a relative path and will be properly resolved).
 *
 * The exact opposite can also be done, excluding ALL parameters except the ones declared with the "::keep" statement (or a file containing one station ID
 * per line followed by a space delimited list of \ref meteoparam "meteorological parameters" to keep for the station).
 *
 * @code
 * [InputEditing]
 * WFJ2::EXCLUDE = HS PSUM                       ;inline declaration of parameters exclusion
 * KLO3::KEEP = TA RH VW DW                      ;inline declaration of parameters to keep
 *
 * EXCLUDE_FILE = ../input/meteo/excludes.csv    ;parameters exclusions defined in a separate file
 * KEEP_FILE = ../input/meteo/keeps.csv          ;parameters to keep defined in a separate file
 * @endcode
 *
 * In the second example (relying on a separate file), the file "../input/meteo/excludes.csv" could look like this:
 * @code
 * WFJ2 TA RH
 * KLO3 HS PSUM
 * @endcode
 *
 * Another example relying on wildcards (the kept/excluded parameters lists are <b>additive</b>):
 * @code
 * [InputEditing]
 * *::KEEP = TA RH                               ;all stations will keep TA and RH and reject the other parameters
 * WFJ2::KEEP = HS PSUM                          ;WFJ2 will keep TA and RH as defined above but also HS and PSUM
 * @endcode
 *
 * @note First the <i>exclude</i> directives are applied and then the <i>keep</i> directives.
 * 
 * @section data_merging 3. Data merging (MERGE)
 * @subsection stations_merging 3.1 Merging different stations (MERGE)
 * It is possible to merge different data sets together, with a syntax similar to the Exclude/Keep syntax. This merging occurs <b>after</b> any
 * EXCLUDE/KEEP commands. This is useful, for example, to provide measurements from different stations that actually share the
 * same measurement location or to build "composite" station from multiple real stations (in this case, using EXCLUDE and/or KEEP
 * commands to fine tune how the composite station(s) is/are built).
 * Please note that the order of declaration defines the priority (ie the first station that has a value for a given parameter has priority). Please also
 * note that only common timestamps will be merged! (ie if the stations have different sampling rates, it might end up that no merge gets performed)
 *
 * @code
 * [Input]
 * METEO = SMET
 * METEOPATH = ./input
 * STATION1 = STB
 * STATION2 = WFJ2
 * STATION3 = WFJ1
 * STATION4 = DAV1
 * [...]
 *
 * [InputEditing]
 * STB::EXCLUDE = ILWR PSUM
 * WFJ2::KEEP = PSUM ILWR RSWR
 *
 * STB::MERGE = WFJ2 WFJ1
 * DAV1::MERGE = WFJ2
 * @endcode
 * In order to avoid circular dependencies, a station can NOT receive data from a station AND contribute data to another station. Otherwise, a
 * station can be merged into multiple other stations. Moreover, the merging strategy can be controlled by setting the MERGE_STRATEGY key in
 * the [InputEditing] section (by default it is "EXPAND_MERGE", see MeteoData::Merge_Type).
 *
 * @note One limitation when handling "extra" parameters (ie parameters that are not in the default \ref meteoparam) is that these extra
 * parameters must be known from the beginning. So if station2 appears later in time with extra parameters, make sure that the buffer size
 * is large enough to reach all the way to this new station (by setting General::BUFFER_SIZE at least to the number of days from
 * the start of the first station to the start of the second station)
 *
 * @subsection automerge 3.2 Automerge
 * If the key \em AUTOMERGE is set to true in the [InputEditing] section, all stations that have identical IDs will be merged together. The first station
 * to come (usually, the first that was defined in the plugin) has the priority over the next ones.
 *
 * @section data_copy 4. Data copy (COPY)
 * It is also possible to duplicate a meteorological parameter as another meteorological parameter. This is done by specifying a COPY key, following the syntax
 * new_name::COPY = existing_parameter. For example:
 * @code
 * [InputEditing]
 * VW_avg::COPY = VW
 * @endcode
 * This creates a new parameter VW_avg that starts as an exact copy of the raw data of VW, for each station. This newly created parameter is
 * then processed as any other meteorological parameter (thus going through filtering, generic processing, spatial interpolations). This only current
 * limitation is that the parameter providing the raw data must be defined for all stations (even if filled with nodata, this is good enough).
 *
 * @section data_creation 5. Data creation (CREATE)
 * Finally, it is possible to create new data based on some parametrizations. If the requested parameter does not exists, it will be created. Otherwise,
 * any pre-existing data is kept and only missing values in the original data set are filled with the generated values, keeping the original sampling rate. As
 * with all raw data editing, this takes place *before* any filtering/resampling/data generators. As the available algorithms are the same as for the
 * data generators, they are listed in the \ref generators_keywords "data generators section" (but the data creators must be declared in the [InputEditing] section).
 * @code
 * [InputEditing]
 * P::create = STD_PRESS			#the pressure is filled with STD_PRESS if no measured values are available
 * ISWR_POT::create = clearSky_SW		#a new parameter "ISWR_POT" is created and filled with Clear Sky values
 * @endcode
 */

DataEditing::DataEditing(const Config& cfgreader)
           : timeproc(cfgreader), cfg(cfgreader), dataCreator(cfgreader), excluded_params(), kept_params(),
             merge_commands(), copy_commands(), move_commands(),
             merged_stations(), merge_strategy(MeteoData::EXPAND_MERGE),
             copy_ready(false), move_ready(false), excludes_ready(false), keeps_ready(false), merge_ready(false), automerge(false)
{
	cfg.getValue("AUTOMERGE", "InputEditing",automerge, IOUtils::nothrow);
	const std::string merge_strategy_str = cfg.get("MERGE_STRATEGY", "InputEditing", "");
	if (!merge_strategy_str.empty())
		merge_strategy = MeteoData::getMergeType(merge_strategy_str);
}

DataEditing& DataEditing::operator=(const DataEditing& source) 
{
	if (this != &source) {
		timeproc = source.timeproc;
		dataCreator = source.dataCreator;
		excluded_params = source.excluded_params;
		kept_params = source.kept_params;
		merge_commands = source.merge_commands;
		merged_stations = source.merged_stations;
		copy_commands = source.copy_commands;
		move_commands = source.move_commands;
		merge_strategy = source.merge_strategy;
		copy_ready = source.copy_ready;
		move_ready = source.move_ready;
		excludes_ready = source.excludes_ready;
		keeps_ready = source.keeps_ready;
		merge_ready = source.merge_ready;
	}
	return *this;
}

void DataEditing::editTimeSeries(STATIONS_SET& vecStation)
{
	if (automerge) automerge_stations(vecStation);

	if (!merge_ready) create_merge_map();
	merge_stations(vecStation);
}

void DataEditing::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	if (automerge) automerge_stations(vecMeteo);

	if (!move_ready) create_move_map();
	move_params(vecMeteo);

	if (!excludes_ready) create_exclude_map();
	exclude_params(vecMeteo);

	if (!keeps_ready) create_keep_map();
	keep_params(vecMeteo);

	if (!merge_ready) create_merge_map();
	merge_stations(vecMeteo);
	
	//remove trailing pure nodata MeteoData elements (if any)
	purgeTrailingNodata(vecMeteo);

	if (!copy_ready) create_copy_map();
	copy_params(vecMeteo);

	timeproc.process(vecMeteo);
	TimeProcStack::checkUniqueTimestamps(vecMeteo);

	dataCreator.createParameters(vecMeteo);
}

//some vector of MeteoData might have trailing elements that are purely nodata
void DataEditing::purgeTrailingNodata(std::vector<METEO_SET>& vecMeteo)
{
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		//purge trailing nodata
		for (size_t jj=vecMeteo[ii].size(); jj>0; jj--) {
			if (!vecMeteo[ii][jj-1].isNodata()) {
				if (jj!=vecMeteo[ii].size()) vecMeteo[ii].resize( jj );
				break;
			}
		}
	}
}

void DataEditing::create_merge_map()
{
	merge_ready = true;

	const std::vector<std::string> merge_keys( cfg.getKeys("::MERGE", "InputEditing", true) );
	const size_t nrOfStations = merge_keys.size();
	for (size_t ii=0; ii<nrOfStations; ++ii) {
		const size_t found = merge_keys[ii].find_first_of(":");
		if (found==std::string::npos) continue;

		const std::string station( IOUtils::strToUpper(merge_keys[ii].substr(0,found)) );
		std::vector<std::string> vecString;
		cfg.getValue(merge_keys[ii], "InputEditing", vecString);
		if (vecString.empty()) throw InvalidArgumentException("Empty value for key \""+merge_keys[ii]+"\"", AT);

		for (vector<string>::iterator it = vecString.begin(); it != vecString.end(); ++it) {
			IOUtils::toUpper( *it );
			const std::vector<std::string>::const_iterator vec_it = find (merged_stations.begin(), merged_stations.end(), *it);
			if (vec_it==merged_stations.end()) merged_stations.push_back( *it ); //this station will be merged into another one
		}
		merge_commands[ station ] = vecString;
	}

	//sort the merged_stations vector so searches will be faster
	std::sort(merged_stations.begin(), merged_stations.end());

	//make sure there is no "chain merge": station A merging station B and station C merging station A
	std::map< std::string, std::vector<std::string> >::iterator it_dest;
	for(it_dest=merge_commands.begin(); it_dest!=merge_commands.end(); ++it_dest) {
		const std::string stationID( it_dest->first );
		if (std::binary_search(merged_stations.begin(), merged_stations.end(), stationID))
			throw InvalidArgumentException("\'chain merge\' detected for station \'"+stationID+"\', this is not supported (see documentation)", AT);
	}
}

//merge stations that have identical names
void DataEditing::merge_stations(STATIONS_SET& vecStation) const
{
	if (merge_commands.empty()) return;

	for (size_t ii=0; ii<vecStation.size(); ii++) {
		const std::string toStationID( IOUtils::strToUpper( vecStation[ii].stationID ) );
		//we do not support "chain merge": station A merging station B and station C merging station A
		if ( std::find(merged_stations.begin(), merged_stations.end(), toStationID)!=merged_stations.end() ) continue;

		const std::map< string, vector<string> >::const_iterator it = merge_commands.find( toStationID );
		if (it == merge_commands.end()) continue; //no merge commands for this station

		const std::vector<std::string> merge_from( it->second );
		for (std::vector<std::string>::const_iterator it_set=merge_from.begin(); it_set != merge_from.end(); ++it_set) {
			const std::string fromStationID( *it_set );

			bool found = false;
			for (size_t jj=0; jj<vecStation.size(); jj++) {
				const std::string curr_station( IOUtils::strToUpper(vecStation[jj].stationID) );
				if (curr_station==fromStationID) {
					vecStation[ii].merge( vecStation[jj] );
					found = true;
				}
			}
			if (!found)
				throw InvalidArgumentException("Station ID '"+fromStationID+"' not found when merging toward station '"+toStationID+"'. Consider increasing BUFFER_SIZE!", AT);
		}
	}

	//remove the stations that have been merged into other ones
	for (size_t ii=0; ii<vecStation.size(); ii++) {
		const std::string stationID( IOUtils::strToUpper( vecStation[ii].stationID ) );
		const std::vector<std::string>::const_iterator it = std::find(merged_stations.begin(), merged_stations.end(), stationID);
		if ( it!=merged_stations.end() ) {
			std::swap( vecStation[ii], vecStation.back() );
			vecStation.pop_back();
			ii--; //in case we have multiple identical stations ID
		}
	}
}

//in this implementation, we consider that the station name does NOT change over time
void DataEditing::merge_stations(std::vector<METEO_SET>& vecVecMeteo) const
{
	if (merge_commands.empty()) return;

	for (size_t ii=0; ii<vecVecMeteo.size(); ii++) { //loop over the stations
		if (vecVecMeteo[ii].empty())  continue;
		const std::string toStationID( IOUtils::strToUpper(vecVecMeteo[ii][0].meta.stationID) );
		//we do not support "chain merge": station A merging station B and station C merging station A
		if ( std::find(merged_stations.begin(), merged_stations.end(), toStationID)!=merged_stations.end() ) continue;

		const std::map< std::string, std::vector<std::string> >::const_iterator it = merge_commands.find( toStationID );
		if (it == merge_commands.end()) continue; //no merge commands for this station

		const std::vector<std::string> merge_from( it->second );
		for (std::vector<std::string>::const_iterator it_set=merge_from.begin(); it_set != merge_from.end(); ++it_set) {
			const std::string fromStationID( *it_set );

			bool found = false;
			for (size_t jj=0; jj<vecVecMeteo.size(); jj++) { //loop over the available stations in the current dataset
				if (vecVecMeteo[jj].empty()) continue;
				const std::string curr_station( IOUtils::strToUpper(vecVecMeteo[jj][0].meta.stationID) );
				if (curr_station==fromStationID) {
					MeteoData::mergeTimeSeries(vecVecMeteo[ii], vecVecMeteo[jj], static_cast<MeteoData::Merge_Type>(merge_strategy)); //merge timeseries for the two stations
					found = true;
				}
			}
			if (!found)
				throw InvalidArgumentException("Station ID '"+fromStationID+"' not found when merging toward station '"+toStationID+"'. Consider increasing BUFFER_SIZE!", AT);
		}
	}

	//remove the stations that have been merged into other ones
	for (size_t ii=0; ii<vecVecMeteo.size(); ii++) {
		if (vecVecMeteo[ii].empty())  continue;
		const std::string stationID( IOUtils::strToUpper(vecVecMeteo[ii][0].meta.stationID) );
		const std::vector<std::string>::const_iterator it = std::find(merged_stations.begin(), merged_stations.end(), stationID);
		if ( it!=merged_stations.end() ) {
			std::swap( vecVecMeteo[ii], vecVecMeteo.back() );
			vecVecMeteo.pop_back();
			ii--; //in case we have multiple identical stations ID
		}
	}
}

//merge stations that have identical IDs
void DataEditing::automerge_stations(STATIONS_SET& vecStation) const
{
	for (size_t ii=0; ii<vecStation.size(); ii++) { //loop over the stations
		const std::string toStationID( IOUtils::strToUpper(vecStation[ii].stationID) );
		for (size_t jj=ii+1; jj<vecStation.size(); jj++) { //loop over the stations
			const std::string fromStationID( IOUtils::strToUpper(vecStation[jj].stationID) );
			if (fromStationID==toStationID) {
				vecStation[ii].merge( vecStation[jj] );
				std::swap( vecStation[jj], vecStation.back() );
				vecStation.pop_back();
				jj--; //we need to redo the current jj, because it contains another station
			}
		}
	}
}

//merge stations that have identical IDs
void DataEditing::automerge_stations(std::vector<METEO_SET>& vecVecMeteo) const
{
	for (size_t ii=0; ii<vecVecMeteo.size(); ii++) { //loop over the stations
		if (vecVecMeteo[ii].empty())  continue;
		const std::string toStationID( IOUtils::strToUpper(vecVecMeteo[ii][0].meta.stationID) );
		size_t nr_conflicts = 0;
		
		for (size_t jj=ii+1; jj<vecVecMeteo.size(); jj++) { //loop over the stations
			if (vecVecMeteo[jj].empty())  continue;
			const std::string fromStationID( IOUtils::strToUpper(vecVecMeteo[jj][0].meta.stationID) );
			if (fromStationID==toStationID) {
				nr_conflicts += MeteoData::mergeTimeSeries(vecVecMeteo[ii], vecVecMeteo[jj], static_cast<MeteoData::Merge_Type>(merge_strategy), MeteoData::CONFLICTS_AVERAGE); //merge timeseries for the two stations
				std::swap( vecVecMeteo[jj], vecVecMeteo.back() );
				vecVecMeteo.pop_back();
				jj--; //we need to redo the current jj, because it contains another station
			}
		}
		
		if (nr_conflicts>0) std::cerr << "[E] " << nr_conflicts << " automerge conflicts on station " <<  toStationID << "\n";
	}
}

void DataEditing::create_exclude_map()
{
	excludes_ready = true;
	const std::string exclude_file = cfg.get("EXCLUDE_FILE", "InputEditing", "");

	if (!exclude_file.empty()) {
		//if this is a relative path, prefix the path with the current path
		const std::string prefix = ( FileUtils::isAbsolutePath(exclude_file) )? "" : cfg.getConfigRootDir()+"/";
		const std::string path( FileUtils::getPath(prefix+exclude_file, true) );  //clean & resolve path
		const std::string filename( path + "/" + FileUtils::getFilename(exclude_file) );

		if (!FileUtils::fileExists(filename)) throw AccessException(filename, AT); //prevent invalid filenames
		std::ifstream fin(filename.c_str(), std::ifstream::in);
		if (fin.fail()) throw AccessException(filename, AT);

		try {
			const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file

			std::vector<std::string> tmpvec;
			std::string line;

			while (!fin.eof()) { //Go through file
				getline(fin, line, eoln); //read complete line meta information
				IOUtils::stripComments(line);
				const size_t ncols = IOUtils::readLineToVec(line, tmpvec, ' ');

				if (ncols > 1) {
					for (std::vector<std::string>::iterator it = tmpvec.begin()+1; it != tmpvec.end(); ++it) {
						IOUtils::toUpper( *it );
					}

					const std::set<std::string> tmpset(tmpvec.begin()+1, tmpvec.end());
					excluded_params[ IOUtils::strToUpper(tmpvec[0]) ] = tmpset;
				}
			}
		} catch (const std::exception&) {
			fin.close();
			throw;
		}

		fin.close();
	}

	const std::vector<std::string> exclude_keys( cfg.getKeys("::EXCLUDE", "InputEditing", true) );
	const size_t nrOfStations = exclude_keys.size();
	for (size_t ii=0; ii<nrOfStations; ++ii) {
		const size_t found = exclude_keys[ii].find_first_of(":");
		if (found==std::string::npos) continue;

		const std::string station( IOUtils::strToUpper(exclude_keys[ii].substr(0,found)) );
		std::vector<std::string> vecString;
		cfg.getValue(exclude_keys[ii], "InputEditing", vecString);
		if (vecString.empty()) throw InvalidArgumentException("Empty value for key \""+exclude_keys[ii]+"\"", AT);
		for (vector<string>::iterator it = vecString.begin(); it != vecString.end(); ++it) {
			IOUtils::toUpper( *it );
		}

		const std::set<std::string> tmpset(vecString.begin(), vecString.end());
		excluded_params[ station ] = tmpset;
	}

	//Handle "*" wildcard: add the params to all other declared stations
	std::map< std::string, std::set<std::string> >::const_iterator it_station = excluded_params.find("*");
	if (it_station!=excluded_params.end()) {
		const std::set<std::string> wildcard( excluded_params["*"] );
		for (it_station=excluded_params.begin(); it_station!=excluded_params.end(); ++it_station) {
			std::set<std::string> params( it_station->second );

			for (std::set<std::string>::const_iterator it=wildcard.begin(); it!=wildcard.end(); ++it)
				params.insert( *it ); //merging: keep in mind that a set can not contain duplicates

			excluded_params[ it_station->first ] = params;
		}
	}
}

void DataEditing::create_keep_map()
{
	keeps_ready = true;
	const std::string keep_file = cfg.get("KEEP_FILE", "InputEditing", "");

	if (!keep_file.empty()) {
		//if this is a relative path, prefix the path with the current path
		const std::string prefix = ( FileUtils::isAbsolutePath(keep_file) )? "" : cfg.getConfigRootDir()+"/";
		const std::string path( FileUtils::getPath(prefix+keep_file, true) );  //clean & resolve path
		const std::string filename( path + "/" + FileUtils::getFilename(keep_file) );

		if (!FileUtils::fileExists(filename)) throw AccessException(filename, AT); //prevent invalid filenames
		std::ifstream fin(filename.c_str(), std::ifstream::in);
		if (fin.fail()) throw AccessException(filename, AT);

		try {
			const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file

			std::vector<std::string> tmpvec;
			std::string line;

			while (!fin.eof()) { //Go through file
				getline(fin, line, eoln); //read complete line meta information
				IOUtils::stripComments(line);
				const size_t ncols = IOUtils::readLineToVec(line, tmpvec, ' ');

				if (ncols > 1) {
					for (vector<string>::iterator it = tmpvec.begin()+1; it != tmpvec.end(); ++it) {
						IOUtils::toUpper( *it );
					}

					const set<string> tmpset(tmpvec.begin()+1, tmpvec.end());
					kept_params[ IOUtils::strToUpper(tmpvec[0]) ] = tmpset;
				}
			}
		} catch (const std::exception&) {
			fin.close();
			throw;
		}

		fin.close();
	}

	const std::vector<std::string> keep_keys( cfg.getKeys("::KEEP", "InputEditing", true) );
	const size_t nrOfStations = keep_keys.size();
	for (size_t ii=0; ii<nrOfStations; ++ii) {
		const size_t found = keep_keys[ii].find_first_of(":");
		if (found==std::string::npos) continue;

		const std::string station( IOUtils::strToUpper(keep_keys[ii].substr(0,found)) );
		std::vector<std::string> vecString;
		cfg.getValue(keep_keys[ii], "InputEditing", vecString);
		if (vecString.empty()) throw InvalidArgumentException("Empty value for key \""+keep_keys[ii]+"\"", AT);
		for (vector<string>::iterator it = vecString.begin(); it != vecString.end(); ++it) {
			IOUtils::toUpper( *it );
		}

		const std::set<std::string> tmpset(vecString.begin(), vecString.end());
		kept_params[ station ] = tmpset;
	}

	//Handle "*" wildcard: add the params to all other declared stations
	std::map< std::string, std::set<std::string> >::const_iterator it_station = kept_params.find("*");
	if (it_station!=kept_params.end()) {
		const std::set<std::string> wildcard( kept_params["*"] );
		for (it_station=kept_params.begin(); it_station!=kept_params.end(); ++it_station) {
			std::set<std::string> params( it_station->second );

			for (std::set<std::string>::const_iterator it=wildcard.begin(); it!=wildcard.end(); ++it)
				params.insert( *it ); //merging: keep in mind that a set can not contain duplicates

			kept_params[ it_station->first ] = params;
		}
	}
}

/**
* @brief reset to nodata the parameters marked as EXCLUDE on a per station basis
*/
void DataEditing::exclude_params(std::vector<METEO_SET>& vecVecMeteo) const
{
	if (excluded_params.empty()) return;

	for (size_t station=0; station<vecVecMeteo.size(); ++station) { //loop over the stations
		if (vecVecMeteo[station].empty()) continue;
		const std::string stationID( IOUtils::strToUpper(vecVecMeteo[station][0].meta.stationID) );
		std::map< std::string, std::set<std::string> >::const_iterator it = excluded_params.find(stationID);
		if (it == excluded_params.end()) {
			it = excluded_params.find("*"); //fallback: is there a wildcard like "*::KEEP"?
			if (it == excluded_params.end()) continue;
		}

		const std::set<std::string> excluded( it->second );

		for (size_t ii=0; ii<vecVecMeteo[station].size(); ++ii) { //loop over the timesteps
			for (std::set<std::string>::const_iterator it_set=excluded.begin(); it_set != excluded.end(); ++it_set) {
				const std::string param( *it_set );
				if (vecVecMeteo[station][ii].param_exists(param))
					vecVecMeteo[station][ii](param) = IOUtils::nodata;
			}
		}
	}
}

/**
* @brief only keep the parameters marked as KEEP on a per station basis
*/
void DataEditing::keep_params(std::vector<METEO_SET>& vecVecMeteo) const
{
	if (kept_params.empty()) return;

	for (size_t station=0; station<vecVecMeteo.size(); ++station) { //loop over the stations
		if (vecVecMeteo[station].empty()) continue;

		const std::string stationID( IOUtils::strToUpper(vecVecMeteo[station][0].meta.stationID) );
		std::map< std::string, std::set<std::string> >::const_iterator it = kept_params.find(stationID);
		if (it == kept_params.end()) {
			it = kept_params.find("*"); //fallback: is there a wildcard like "*::KEEP"?
			if (it == kept_params.end()) continue;
		}

		const std::set<std::string> kept( it->second );

		for (size_t ii=0; ii<vecVecMeteo[station].size(); ++ii) {
			MeteoData& md_ref( vecVecMeteo[station][ii] );
			MeteoData md( md_ref );
			md.reset(); //delete all meteo fields

			for (std::set<std::string>::const_iterator it_set=kept.begin(); it_set != kept.end(); ++it_set) { //loop over the parameters to keep
				const std::string param( *it_set);
				if (!md.param_exists(param)) continue;
				 md(param) = md_ref(param);
			}

			//copy back the new object into vecVecMeteo
			md_ref = md;
		}
	}
}

/**
* Parse [InputEditing] section for potential parameters that the user wants
* renamed (as '%%::MOVE = %%')
*/
void DataEditing::create_move_map()
{
	const std::vector<std::string> move_keys( cfg.getKeys("::MOVE", "InputEditing", true) );
	const size_t nrOfMatches = move_keys.size();

	for (size_t ii=0; ii<nrOfMatches; ++ii) {
		const std::string dest_param( move_keys[ii].substr( 0, move_keys[ii].find_first_of(":") ) );
		std::vector<std::string> vecString;
		cfg.getValue(move_keys[ii], "InputEditing", vecString); //multiple source can be provided

		if (vecString.empty()) throw InvalidArgumentException("Empty value for key \""+move_keys[ii]+"\"", AT);
		for (vector<string>::iterator it = vecString.begin(); it != vecString.end(); ++it) {
			IOUtils::toUpper( *it );
		}

		const std::set<std::string> tmpset(vecString.begin(), vecString.end());
		move_commands[ dest_param ] = tmpset;
	}

	move_ready = true;
}

/**
* This procedure runs through the MeteoData objects in vecMeteo and according to user
* configuration renames a certain present meteo parameter to another one, named by the
* user in the [InputEditing] section of the io.ini, e.g.
* [InputEditing]
* TA::MOVE = air_temp air_temperature
* means that TA will be the name of a new parameter in MeteoData with the copied value
* of the original parameter air_temp or air_temperature
*
* If TA already had some values, it will keep those and only points not having a value will
* be filled by either air_temp or air_temperature (the first one in the list to have a value
* has the priority)
*/
void DataEditing::move_params(std::vector< METEO_SET >& vecMeteo) const
{
	if (move_commands.empty()) return; //Nothing configured

	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (vecMeteo[station].empty()) continue;

		std::map< std::string, std::set<std::string> >::const_iterator param = move_commands.begin();
		for (; param!=move_commands.end(); ++param) { //loop over all the MOVE commands
			const std::string dest_param( param->first );
			const std::set<std::string> src( param->second );

			for (std::set<std::string>::const_iterator it_set=src.begin(); it_set != src.end(); ++it_set) { //loop over the parameters to move
				const std::string src_param( *it_set );
				const size_t src_index = vecMeteo[station].front().getParameterIndex( src_param );
				if (src_index == IOUtils::npos) continue; //no such parameter for this station, skipping

				for (size_t jj=0; jj<vecMeteo[station].size(); ++jj) {
					const size_t dest_index = vecMeteo[station][jj].addParameter( dest_param ); //either add or just return the proper index
					if (vecMeteo[station][jj]( dest_index ) == IOUtils::nodata) {
						vecMeteo[station][jj]( dest_index ) = vecMeteo[station][jj]( src_index );
						vecMeteo[station][jj]( src_index ) = IOUtils::nodata;
					}
				}
			}
		}
	}
}


/**
* Parse [InputEditing] section for potential parameters that the user wants
* duplicated (as '%%::COPY = %%')
*/
void DataEditing::create_copy_map()
{
	const std::vector<std::string> copy_keys( cfg.getKeys("::COPY", "InputEditing", true) );
	const size_t nrOfMatches = copy_keys.size();

	for (size_t ii=0; ii<nrOfMatches; ++ii) {
		const std::string dest_param( copy_keys[ii].substr( 0, copy_keys[ii].find_first_of(":") ) );
		const std::string src_param = cfg.get(copy_keys[ii], "InputEditing");
		if (!dest_param.empty() && !src_param.empty())
			copy_commands[ dest_param ] = src_param;
	}

	copy_ready = true;
}

/**
* This procedure runs through the MeteoData objects in vecMeteo and according to user
* configuration copies a certain present meteo parameter to another one, named by the
* user in the [InputEditing] section of the io.ini, e.g.
* [InputEditing]
* TA2::COPY = TA
* means that TA2 will be the name of a new parameter in MeteoData with the copied value
* of the meteo parameter MeteoData::TA
*/
void DataEditing::copy_params(std::vector< METEO_SET >& vecMeteo) const
{
	if (copy_commands.empty()) return; //Nothing configured

	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (vecMeteo[station].empty()) continue;

		std::map< std::string, std::string >::const_iterator param = copy_commands.begin();
		for (; param!=copy_commands.end(); ++param) {
			const size_t src_index = vecMeteo[station].front().getParameterIndex(param->second);
			if (src_index == IOUtils::npos) {
				const std::string stationID( vecMeteo[station].front().meta.stationID );
				throw InvalidArgumentException("Station "+stationID+" has no parameter '"+param->second+"' to copy", AT);
			}

			const std::string dest( param->first );
			for (size_t jj=0; jj<vecMeteo[station].size(); ++jj) { //for each MeteoData object of one station
				const size_t dest_index = vecMeteo[station][jj].addParameter( dest );
				vecMeteo[station][jj]( dest_index ) = vecMeteo[station][jj]( src_index );
			}
		}
	}
}

const std::string DataEditing::toString() const
{
	std::ostringstream os;
	os << "<DataEditing>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";

	if (!excluded_params.empty()) {
		os << "<excluded_params>\n";
		std::map< std::string, std::set<std::string> >::const_iterator it_exc;
		for (it_exc=excluded_params.begin(); it_exc != excluded_params.end(); ++it_exc) {
			os << setw(10) << it_exc->first << " = ";
			std::set<std::string>::const_iterator it_set;
			for (it_set=(it_exc->second).begin(); it_set != (it_exc->second).end(); ++it_set)
				os << *it_set << " ";
			os << "\n";
		}
		os << "</excluded_params>\n";
	}

	if (!merge_commands.empty()) {
		os << "<merge_commands>\n";
		std::map< std::string, std::vector<std::string> >::const_iterator it_merge;
		for (it_merge=merge_commands.begin(); it_merge != merge_commands.end(); ++it_merge) {
			os << setw(10) << it_merge->first << " <- ";
			for (size_t ii=0; ii<it_merge->second.size(); ++ii)
				os << it_merge->second[ii] << " ";
			os << "\n";
		}
		os << "</merge_commands>\n";
	}

	os << dataCreator.toString();

	os << "</DataEditing>\n";
	return os.str();
}

} //end namespace
