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
 * edition commands that can be stacked at will, per station ID. This is similar to the way that filters (processing elements) are
 * also stacked together. The general syntax is ('#' represent a number, so each key remains unique):
 * @code
 * {stationID}::edit# = {command}
 * {stationID}::arg#::{argument name} = {values}
 * 
 * #here is an example
 * WFJ2::edit1 = EXCLUDE
 * WFJ2::arg1::exclude = VW DW ISWR RSWR
 * @endcode
 * 
 * It is also possible to apply a stack of edits to all stations by using the '*' wildcard instead of the station ID 
 * (in such a case, the wildcard stack will be applied before any other stack). Currently, the data creation is applied 
 * after all stacks have been processed but this will change in the near future...
 * 
 * The following Input Data Editing commands are available:
 *     - SWAP: swap two parameters, see EditingSwap
 *     - MOVE: rename one or more parameters into a new name, see EditingMove
 *     - EXCLUDE: delete a list of parameters, see EditingExclude
 *     - KEEP: only keep a list of parameters and reject the others, see EditingKeep
 *     - AUTOMERGE: merge together stations sharing the same station ID, see EditingAutoMerge
 *     - MERGE: merge together one or more stations, see EditingMerge
 *     - COPY: make a copy of a given parameter under a new name, see EditingCopy
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
            : EditingBlock(i_stationID, vecArgs, name), merged_stations(), merge_strategy(MeteoData::EXPAND_MERGE), merge_conflicts(MeteoData::CONFLICTS_PRIORITY)
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
		} else if (vecArgs[ii].first=="MERGE_CONFLICTS") {
			merge_conflicts = MeteoData::getMergeConflicts( vecArgs[ii].second );
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
			
			MeteoData::mergeTimeSeries( vecMeteo[toStationIdx], vecMeteo[ii], merge_strategy, merge_conflicts );
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
