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
#include <meteoio/DataEditingAlgorithms.h>
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
 * Before any filters, resampling algorithms or data generators are applied, it is possible to edit the original data. There are several
 * edition commands that can be stacked at will, per station ID. It is also possible to apply a stack of edits to all stations
 * by using the '*' wildcard instead of the station ID (in such a case, the wildcard stack will be applied before any other
 * stack). Currently, the data creation is applied after all stacks have been processed but this will change in the near future...
 *
 * The editing command that are available are the following:
 *     - \ref data_renaming "rename certain parameters;"
 *     - \ref data_exclusion "exclude/keep certain parameters on a per station basis;"
 *     - \ref data_merging "merge stations together;"
 *     - \ref data_copy "make a copy of a certain parameter under a new parameter name for all stations;"
 *     - \ref data_creation "create certain parameters based on some parametrizations."
 *
 * The general syntax is <i>{stationID}\:\:edit# = {command}</i> followed by the command's arguments as 
 * <i>{stationID}\:\:arg#::{argument name} = {values}</i> where '#' represent a number (so each key remains unique).
 * 
 * @section data_renaming 1. Data renaming (MOVE or SWAP)
 * @subsection data_move 1.1 Data renaming (MOVE)
 * It is possible to rename a meteorological parameter thanks to the MOVE key. This key can take 
 * multiple source names that will be processed in the order of declaration. Original names that are not found in the current
 * dataset will silently be ignored, so it is safe to provide a list that contain many possible names:
 * @code
 * [InputEditing]
 * SLF2::edit1 = MOVE
 * SLF2::arg1::dest = TA
 * SLF2::arg1::src = air_temp air_temperature temperature_air
 * @endcode
 * This can be used to rename non-standard parameter names into standard ones. In this example, if TA already had some values, it will keep
 * those and only points not having a value will be filled by either air_temp or air_temperature or temperature_air (the first one in
 * the list to have a value has the priority).
 * 
 * @subsection data_swap 1.2 Swapping parameters (SWAP)
 * It is possible to swap pairs of parameters with the SWAP key. This supports both standard \ref meteoparam "meteorological parameters" as well
 * as non-standard parameters (ie not in the list in the link). If a parameter does not exists, it will be transparently added with a nodata value.
 * 
 * @code
 * [InputEditing]
 * FLU2::edit2 = SWAP
 * FLU2::arg2::dest = ISWR
 * FLU2::arg2::src = RSWR
 * @endcode
 *
 * @section data_exclusion 2. Data exclusion (EXCLUDE/KEEP)
 * It is possible to exclude specific parameters with the "exclude" command. This is either done by providing a space delimited list of 
 * \ref meteoparam "meteorological parameters" to exclude for the station as key.
 *
 * The exact opposite can also be done, excluding ALL parameters except the ones declared with the "keep" command.
 *
 * @code
 * [InputEditing]
 * FLU2::edit3 = EXCLUDE
 * FLU2::arg3::exclude = TA RH TSS TSG
 * @endcode
 *
 * Another example relying on wildcards (the kept/excluded parameters lists are <b>additive</b>):
 * @code
 * [InputEditing]
 * *::edit1 = KEEP                               ;all stations will keep TA and RH and reject the other parameters
 * *::arg1::keep = TA RH
 * 
 * WFJ2::edit1 = KEEP                          ;WFJ2 will keep TA and RH as defined above but also HS and PSUM
 * WFJ2::arg1::keep = HS PSUM
 * @endcode
 * 
 * @section data_merging 3. Data merging (MERGE)
 * @subsection stations_merging 3.1 Merging different stations (MERGE)
 * It is possible to merge different data sets together, with the MERGE command. This is useful, for example, to 
 * provide measurements from different stations that actually share the same measurement location or to build 
 * "composite" station from multiple real stations (in this case, using EXCLUDE and/or KEEP commands to fine tune 
 * how the composite station(s) is/are built).
 * Please note that the order of declaration defines the priority (ie the first station that has a value for a given parameter has priority). Please also
 * note that which timestamps will be merged depends on the chosen merge strategy (see MeteoData::MERGE_TYPE).
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
 * STB::edit1 = EXCLUDE
 * STB::arg1::exclude = ILWR PSUM
 * 
 * WFJ2::edit1 = KEEP
 * WFJ2::arg1::keep = PSUM ILWR RSWR
 *
 * STB::edit2 = MERGE
 * STB::arg2::merge = WFJ2 WFJ1
 * STB::arg2::merge_strategy = FULL_MERGE
 * 
 * DAV1::edit1 = MERGE
 * DAV1::arg1::merge = WFJ2
 * @endcode
 * In order to avoid circular dependencies, a station can NOT receive data from a station AND contribute data to another station. Otherwise, a
 * station can be merged into multiple other stations. Moreover, the merging strategy can be controlled by setting the MERGE_STRATEGY 
 * optional argument (by default it is "EXPAND_MERGE", see MeteoData::Merge_Type).
 *
 * @note One limitation when handling "extra" parameters (ie parameters that are not in the default \ref meteoparam) is that these extra
 * parameters must be known from the beginning. So if station2 appears later in time with extra parameters, make sure that the buffer size
 * is large enough to reach all the way to this new station (by setting General::BUFFER_SIZE at least to the number of days from
 * the start of the first station to the start of the second station)
 *
 * @subsection automerge 3.2 Automerge
 * This is a special case of merge: only station's have the exact same ID will get merge together. This is useful when reading data
 * for the same station from multiple source in order to rebuild a consistent dataset. If merge conflicts are encountered (such as 
 * identical fields having different values at the same timestamp), warnings will be printed out.
 * 
 * @code
 * [InputEditing]
 * *::edit1 = AUTOMERGE                        ;all stations having the same ID will be merged together
 * @endcode
 *
 * @section data_copy 4. Data copy (COPY)
 * It is also possible to duplicate a meteorological parameter as another meteorological parameter. This is done with the COPY command, 
 * such as:
 * @code
 * [InputEditing]
 * DAV::edit1 = COPY
 * DAV::arg1::dest = TA_copy
 * DAV::arg1::src = TA
 * @endcode
 * This creates a new parameter TA_copy that starts as an exact copy of the raw data of TA, for the DAV station. This newly created parameter is
 * then processed as any other meteorological parameter (thus going through filtering, generic processing, spatial interpolations). 
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

const std::string EditingBlock::toString() const 
{
	std::ostringstream os;
	os << "[" << stationID << " - " << block_name << "]";
	return os.str();
}

EditingBlock* EditingBlockFactory::getBlock(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
{
	if (name == "SWAP"){
		return new EditingSwap(i_stationID, vecArgs, name);
	} else if (name == "MOVE"){
		return new EditingMove(i_stationID, vecArgs, name);
	} else if (name == "EXCLUDE"){
		return new EditingExclude(i_stationID, vecArgs, name);
	} else if (name == "KEEP"){
		return new EditingKeep(i_stationID, vecArgs, name);
	} else if (name == "MERGE"){
		return new EditingMerge(i_stationID, vecArgs, name);
	} else if (name == "AUTOMERGE"){
		return new EditingAutoMerge(i_stationID, vecArgs, name);
	} if (name == "COPY"){
		return new EditingCopy(i_stationID, vecArgs, name);
	} else {
		throw IOException("The input data editing block '"+name+"' does not exist! " , AT);
	}
}


////////////////////////////////////////////////// SWAP
EditingSwap::EditingSwap(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
            : EditingBlock(i_stationID, vecArgs, name), dest_param(), src_param()
{
	parse_args(vecArgs);
}

void EditingSwap::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );
	bool has_dest=false, has_src=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="DEST") {
			IOUtils::parseArg(vecArgs[ii], where, dest_param);
			IOUtils::toUpper( dest_param );
			has_dest = true;
		} else if (vecArgs[ii].first=="SRC") {
			IOUtils::parseArg(vecArgs[ii], where, src_param);
			IOUtils::toUpper( src_param ),
			has_src = true;
		}
	}

	if (!has_dest) throw InvalidArgumentException("Please provide a DEST value for "+where, AT);
	if (!has_src) throw InvalidArgumentException("Please provide an SRC value for "+where, AT);
}

void EditingSwap::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (vecMeteo[station].empty()) continue;
		if (stationID!="*" && stationID!=IOUtils::strToUpper(vecMeteo[station][0].meta.stationID)) continue;
		
		for (size_t ii=0; ii<vecMeteo[station].size(); ++ii) {
			const size_t src_index = vecMeteo[station][ii].addParameter( src_param ); //either add or just return the proper index
			const double src_value = vecMeteo[station][ii]( src_param );
			
			const size_t dest_index = vecMeteo[station][ii].addParameter( dest_param ); //either add or just return the proper index
			
			vecMeteo[station][ii]( src_index ) = vecMeteo[station][ii]( dest_index );
			vecMeteo[station][ii]( dest_index ) = src_value;
		}
	}
}


////////////////////////////////////////////////// MOVE
EditingMove::EditingMove(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
            : EditingBlock(i_stationID, vecArgs, name), src_params(), dest_param()
{
	parse_args(vecArgs);
}

void EditingMove::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );
	bool has_dest=false, has_src=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="DEST") {
			IOUtils::parseArg(vecArgs[ii], where, dest_param);
			IOUtils::toUpper( dest_param );
			has_dest = true;
		} else if (vecArgs[ii].first=="SRC") {
			IOUtils::readLineToSet( IOUtils::strToUpper(vecArgs[ii].second), src_params);
			has_src = true;
		}
	}

	if (!has_dest) throw InvalidArgumentException("Please provide a DEST value for "+where, AT);
	if (!has_src || src_params.empty()) throw InvalidArgumentException("Please provide a valid SRC value for "+where, AT);
}

void EditingMove::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (vecMeteo[station].empty()) continue; //no data for this station
		if (stationID!="*" && stationID!=IOUtils::strToUpper(vecMeteo[station][0].meta.stationID)) continue;
		
		for (std::set<std::string>::const_iterator it_set=src_params.begin(); it_set != src_params.end(); ++it_set) { //loop over the parameters to move
			const size_t src_index = vecMeteo[station].front().getParameterIndex( *it_set );
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


////////////////////////////////////////////////// EXCLUDE
EditingExclude::EditingExclude(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
            : EditingBlock(i_stationID, vecArgs, name), exclude_params()
{
	parse_args(vecArgs);
}

void EditingExclude::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );
	bool has_excludes=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="EXCLUDE") {
			IOUtils::readLineToSet( IOUtils::strToUpper(vecArgs[ii].second), exclude_params);
			has_excludes = true;
		}
	}

	if (!has_excludes || exclude_params.empty()) throw InvalidArgumentException("Please provide a valid EXCLUDE value for "+where, AT);
}

void EditingExclude::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (vecMeteo[station].empty()) continue; //no data for this station
		if (stationID!="*" && stationID!=IOUtils::strToUpper(vecMeteo[station][0].meta.stationID)) continue;
		
		for (size_t ii=0; ii<vecMeteo[station].size(); ++ii) { //loop over the timesteps
			for (std::set<std::string>::const_iterator it_set=exclude_params.begin(); it_set != exclude_params.end(); ++it_set) { //loop over the parameters to exclude
				const std::string param( *it_set );
				if (vecMeteo[station][ii].param_exists(param))
					vecMeteo[station][ii](param) = IOUtils::nodata;
			}
		}
	}
}


////////////////////////////////////////////////// KEEP
EditingKeep::EditingKeep(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
            : EditingBlock(i_stationID, vecArgs, name), keep_params()
{
	parse_args(vecArgs);
}

void EditingKeep::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );
	bool has_excludes=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="KEEP") {
			IOUtils::readLineToSet( IOUtils::strToUpper(vecArgs[ii].second), keep_params);
			has_excludes = true;
		}
	}

	if (!has_excludes || keep_params.empty()) throw InvalidArgumentException("Please provide a valid KEEP value for "+where, AT);
}

void EditingKeep::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (vecMeteo[station].empty()) continue; //no data for this station
		if (stationID!="*" && stationID!=IOUtils::strToUpper(vecMeteo[station][0].meta.stationID)) continue;
		
		for (size_t ii=0; ii<vecMeteo[station].size(); ++ii) {//loop over the timesteps
			MeteoData& md_ref( vecMeteo[station][ii] );
			MeteoData md( md_ref );
			md.reset(); //delete all meteo fields

			for (std::set<std::string>::const_iterator it_set=keep_params.begin(); it_set != keep_params.end(); ++it_set) { //loop over the parameters to keep
				const std::string param( *it_set);
				if (!md.param_exists(param)) continue;
				md(param) = md_ref(param);
			}

			//copy back the new object into vecMeteo
			md_ref = md;
		}
	}
}


////////////////////////////////////////////////// AUTOMERGE
EditingAutoMerge::EditingAutoMerge(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
            : EditingBlock(i_stationID, vecArgs, name), merge_strategy(MeteoData::EXPAND_MERGE)
{
	parse_args(vecArgs);
}

void EditingAutoMerge::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="MERGE_STRATEGY") {
			merge_strategy = MeteoData::getMergeType( vecArgs[ii].second );
		}
	}
}

void EditingAutoMerge::mergeStations(const size_t& toStationIdx, STATIONS_SET& vecStation)
{
	const std::string toStationID( vecStation[toStationIdx].stationID );
	
	//stations before toStationIdx are not == stationID both for the "*" station and for any specific stationID
	for (size_t jj=toStationIdx+1; jj<vecStation.size(); jj++) { //loop over the stations
		const std::string fromStationID( IOUtils::strToUpper(vecStation[jj].stationID) );
		if (fromStationID==toStationID) {
			vecStation[toStationIdx].merge( vecStation[jj] );
			std::swap( vecStation[jj], vecStation.back() );
			vecStation.pop_back();
			jj--; //we need to redo the current jj, because it contains another station
		}
	}
}

void EditingAutoMerge::editTimeSeries(STATIONS_SET& vecStation)
{
	if (stationID=="*") {
		for (size_t ii=0; ii<vecStation.size(); ii++) {
			mergeStations(ii, vecStation);
		}
	} else {
		//find our current station in vecStation
		size_t toStationIdx=IOUtils::npos;
		for (size_t ii=0; ii<vecStation.size(); ii++) {
			if (vecStation[ii].stationID==stationID) {
				toStationIdx = ii;
				break;
			}
		}
		
		if (toStationIdx==IOUtils::npos) return;
		mergeStations(toStationIdx, vecStation);
	}
}

void EditingAutoMerge::mergeMeteo(const size_t& toStationIdx, std::vector<METEO_SET>& vecMeteo) const
{
	const std::string toStationID( vecMeteo[toStationIdx].front().getStationID() );
	size_t nr_conflicts = 0;
	
	for (size_t jj=toStationIdx+1; jj<vecMeteo.size(); jj++) { //loop over the stations
		if (vecMeteo[jj].empty())  continue;
		const std::string fromStationID( IOUtils::strToUpper(vecMeteo[jj].front().getStationID()) );
		if (fromStationID==toStationID) {
			nr_conflicts += MeteoData::mergeTimeSeries(vecMeteo[toStationIdx], vecMeteo[jj], merge_strategy, MeteoData::CONFLICTS_AVERAGE); //merge timeseries for the two stations
			std::swap( vecMeteo[jj], vecMeteo.back() );
			vecMeteo.pop_back();
			jj--; //we need to redo the current jj, because it contains another station
		}
	}
	
	if (nr_conflicts>0) std::cerr << "[E] " << nr_conflicts << " automerge conflicts on station " <<  toStationID << "\n";
}

void EditingAutoMerge::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	if (stationID=="*") {
		for (size_t ii=0; ii<vecMeteo.size(); ii++) {
			if (vecMeteo[ii].empty()) continue;
			mergeMeteo(ii, vecMeteo);
		}
	} else {
		//find our current station in vecMeteo
		size_t toStationIdx=IOUtils::npos;
		for (size_t ii=0; ii<vecMeteo.size(); ii++) {
			if (vecMeteo[ii].empty()) continue;
			if (vecMeteo[ii].front().getStationID()==stationID) {
				toStationIdx = ii;
				break;
			}
		}
		
		if (toStationIdx==IOUtils::npos) return;
		//stations before toStationIdx are not == stationID, see above
		mergeMeteo(toStationIdx, vecMeteo);
	}
}


////////////////////////////////////////////////// MERGE
EditingMerge::EditingMerge(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
            : EditingBlock(i_stationID, vecArgs, name), merged_stations(), merge_strategy(MeteoData::EXPAND_MERGE)
{
	if (i_stationID=="*")
		throw InvalidArgumentException("It is not possible to do a MERGE on the '*' stationID", AT);
	
	parse_args(vecArgs);
}

void EditingMerge::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );
	bool has_merges=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="MERGE") {
			IOUtils::readLineToVec( IOUtils::strToUpper(vecArgs[ii].second), merged_stations);
			has_merges = true;
		} else if (vecArgs[ii].first=="MERGE_STRATEGY") {
			merge_strategy = MeteoData::getMergeType( vecArgs[ii].second );
		}
	}
	
	//check that each station ID to merge from is only included once
	const std::set<std::string> tmp(merged_stations.begin(), merged_stations.end());
	if (tmp.size()<merged_stations.size())
		throw InvalidArgumentException("Each station to merge from can only appear once in the list for "+where, AT);
	
	//check that the station does not merge with itself
	if (tmp.count(stationID)>0)
		throw InvalidArgumentException("A station can not merge with itself! Wrong argument in "+where, AT);
	
	if (!has_merges || merged_stations.empty()) throw InvalidArgumentException("Please provide a valid MERGE value for "+where, AT);
}

void EditingMerge::editTimeSeries(STATIONS_SET& vecStation)
{
	//find our current station in vecStation
	size_t toStationIdx=IOUtils::npos;
	for (size_t ii=0; ii<vecStation.size(); ii++) {
		if (vecStation[ii].stationID==stationID) {
			toStationIdx = ii;
			break;
		}
	}
	
	if (toStationIdx==IOUtils::npos) return;
	
	for (size_t jj=0; jj<merged_stations.size(); jj++) {
		const std::string fromStationID( IOUtils::strToUpper( merged_stations[jj] ) );
		
		for (size_t ii=0; ii<vecStation.size(); ii++) {
			if (vecStation[ii].stationID==fromStationID)
				vecStation[toStationIdx].merge( vecStation[ii] );
		}
	}
}

void EditingMerge::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	//find our current station in vecStation
	size_t toStationIdx=IOUtils::npos;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		if (vecMeteo[ii].empty()) continue;
		if (vecMeteo[ii].front().getStationID()==stationID) {
			toStationIdx = ii;
			break;
		}
	}
	
	if (toStationIdx==IOUtils::npos) return;
	
	for (size_t jj=0; jj<merged_stations.size(); jj++) {
		const std::string fromStationID( IOUtils::strToUpper( merged_stations[jj] ) );
		
		for (size_t ii=0; ii<vecMeteo.size(); ii++) {
			if (vecMeteo[ii].empty()) continue;
			if (vecMeteo[ii].front().getStationID()!=fromStationID) continue;
			
			MeteoData::mergeTimeSeries( vecMeteo[toStationIdx], vecMeteo[ii], merge_strategy );
		}
	}
}

std::set<std::string> EditingMerge::getDependencies() const
{
	const std::set<std::string> stations_to_purge(merged_stations.begin(), merged_stations.end());
	return stations_to_purge;
}


////////////////////////////////////////////////// COPY
EditingCopy::EditingCopy(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name)
            : EditingBlock(i_stationID, vecArgs, name), dest_param(), src_param()
{
	parse_args(vecArgs);
}

void EditingCopy::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "InputEditing::"+block_name+" for station "+stationID );
	bool has_dest=false, has_src=false;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="DEST") {
			IOUtils::parseArg(vecArgs[ii], where, dest_param);
			IOUtils::toUpper( dest_param ); //HACK not needed if vecArgs is prepared upper case
			has_dest = true;
		} else if (vecArgs[ii].first=="SRC") {
			IOUtils::parseArg(vecArgs[ii], where, src_param);
			IOUtils::toUpper( src_param ),
			has_src = true;
		}
	}

	if (!has_dest) throw InvalidArgumentException("Please provide a DEST value for "+where, AT);
	if (!has_src) throw InvalidArgumentException("Please provide an SRC value for "+where, AT);
}

void EditingCopy::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (vecMeteo[station].empty()) continue;
		if (stationID!="*" && stationID!=IOUtils::strToUpper(vecMeteo[station][0].meta.stationID)) continue;
		
		for (size_t ii=0; ii<vecMeteo[station].size(); ++ii) {
			const size_t dest_index = vecMeteo[station][ii].addParameter( dest_param ); //either add or just return the proper index
			
			if (vecMeteo[station][ii].param_exists(src_param))
				vecMeteo[station][ii](dest_index) = vecMeteo[station][ii].param_exists(src_param);
		}
	}
}

} //end namespace
