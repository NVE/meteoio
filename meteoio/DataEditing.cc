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

const char DataEditing::NUM[] = "0123456789";
const std::string DataEditing::command_key( "::EDIT" );
const std::string DataEditing::arg_key( "::ARG" );

DataEditing::DataEditing(const Config& cfgreader)
           : timeproc(cfgreader), dataCreator(cfgreader), editingStack()
{
	//cfgreader.getValue("AUTOMERGE", "InputEditing",automerge, IOUtils::nothrow);
	const std::set<std::string> editableStations( getEditableStations(cfgreader) );
	
	for (std::set<std::string>::const_iterator it = editableStations.begin(); it != editableStations.end(); ++it) {
		editingStack[ *it ] = buildStack(*it, cfgreader);
	}
}

DataEditing::~DataEditing() 
{
	std::map< std::string, std::vector< EditingBlock* > >::const_iterator it;
	for (it = editingStack.begin(); it != editingStack.end(); ++it) {
		for (size_t ii=0; ii<it->second.size(); ii++)
			delete it->second[ii];
	}
}

DataEditing& DataEditing::operator=(const DataEditing& source) 
{
	if (this != &source) {
		timeproc = source.timeproc;
		dataCreator = source.dataCreator;
		editingStack = source.editingStack;
	}
	
	return *this;
}

std::set<std::string> DataEditing::getEditableStations(const Config& cfg)
{
	const std::vector<std::string> vec_keys( cfg.getKeys(command_key, "INPUTEDITING", true) );

	std::set<std::string> set_stations;
	for (size_t ii=0; ii<vec_keys.size(); ++ii){
		const size_t found = vec_keys[ii].find_first_of(":");
		if (found != std::string::npos){
			if (vec_keys[ii].length()<=(found+2))
				throw InvalidFormatException("Invalid syntax: \""+vec_keys[ii]+"\"", AT);
			if (vec_keys[ii][found+1]!=':')
				throw InvalidFormatException("Missing ':' in \""+vec_keys[ii]+"\"", AT);
				
			const std::string tmp( vec_keys[ii].substr(0,found) );
			set_stations.insert(tmp);
		}
	}

	return set_stations;
}

std::vector< std::pair<std::string, std::string> > DataEditing::parseArgs(const Config& cfg, const std::string& cmd_key, const std::string& stationID)
{
	//extract the cmd number and perform basic checks on the syntax
	const size_t end_cmd = cmd_key.find(command_key); //we know this will be found since it has been matched in cfg.getValues()
	const size_t start_cmd_nr = cmd_key.find_first_of(NUM, end_cmd+command_key.length());
	const size_t end_cmd_nr = cmd_key.find_first_not_of(NUM, end_cmd+command_key.length());
	if (start_cmd_nr==std::string::npos || end_cmd_nr!=std::string::npos) throw InvalidArgumentException("Syntax error: "+cmd_key, AT);

	unsigned int cmd_nr;
	const std::string cmd_nr_str( cmd_key.substr(start_cmd_nr) );
	if ( !IOUtils::convertString(cmd_nr, cmd_nr_str) ) InvalidArgumentException("Can not parse command number in "+cmd_key, AT);

	//read the arguments and clean them up (ie remove the {stationID}::{args##}:: in front of the argument cmd_key itself)
	std::ostringstream arg_str;
	arg_str << stationID << arg_key << cmd_nr;
	std::vector< std::pair<std::string, std::string> > vecArgs( cfg.getValues(arg_str.str(), "INPUTEDITING") );
	for (size_t jj=0; jj<vecArgs.size(); jj++) {
		const size_t beg_arg_name = vecArgs[jj].first.find_first_not_of(":", arg_str.str().length());
		if (beg_arg_name==std::string::npos)
			throw InvalidFormatException("Wrong argument format for '"+vecArgs[jj].first+"'", AT);
		vecArgs[jj].first = vecArgs[jj].first.substr(beg_arg_name);
	}
	
	return vecArgs;
}

std::vector< EditingBlock* > DataEditing::buildStack(const std::string& station_ID, const Config& cfg)
{
	//extract each filter and its arguments, then build the filter stack
	const std::vector< std::pair<std::string, std::string> > vecCommands( cfg.getValues(station_ID+command_key, "INPUTEDITING") );
	std::vector< EditingBlock* > cmd_stack;
	cmd_stack.reserve( vecCommands.size() );
	
	for (size_t ii=0; ii<vecCommands.size(); ii++) {
		const std::string cmd_name( IOUtils::strToUpper( vecCommands[ii].second ) );
		if (cmd_name=="NONE") continue;
		
		const std::vector< std::pair<std::string, std::string> > vecArgs( parseArgs(cfg, vecCommands[ii].first, station_ID) );
		cmd_stack.push_back( EditingBlockFactory::getBlock(station_ID, vecArgs, cmd_name) );
	}
	
	return cmd_stack;
}

std::set<std::string> DataEditing::getMergedFromIDs() const
{
	//build the list of stations that are merged to and merged from
	std::set<std::string> mergedToIDs, mergedFromIDs;
	std::map< std::string, std::vector< EditingBlock* > >::const_iterator it_blocks;
	for (it_blocks = editingStack.begin(); it_blocks != editingStack.end(); ++it_blocks) {
		for (size_t jj=0; jj<it_blocks->second.size(); jj++) {
			const std::set<std::string> tmp_set( it_blocks->second[jj]->getPurgeIDs() );
			if (!tmp_set.empty()) {
				mergedToIDs.insert( it_blocks->first );
				mergedFromIDs.insert(tmp_set.begin(), tmp_set.end());
			}
		}
	}
	
	//make sure there is no "circular merge": station A merging station B and station C merging station A
	//we simply make sure that no destination station ID is also a source data for a merge
	for (std::set<std::string>::const_iterator it = mergedToIDs.begin(); it != mergedToIDs.end(); ++it) {
		if (mergedFromIDs.count( *it )>0)
			throw InvalidArgumentException("\'chain merge\' detected for station \'"+*it+"\', this is not supported (see documentation)", AT);
	}
	
	return mergedFromIDs;
}

void DataEditing::editTimeSeries(STATIONS_SET& vecStation)
{
	//check for "circular merge", build the list of stations that will be purged in the end 
	const std::set<std::string> mergedFromIDs( getMergedFromIDs() );
	
	//process widlcard commands first, knowing that '*' merges are prohibited
	if (editingStack.count("*")>0) {
		for (size_t jj=0; jj<editingStack["*"].size(); jj++)
			editingStack["*"][jj]->editTimeSeries(vecStation);
	}
	
	for (size_t ii=0; ii<vecStation.size(); ii++) {
		const std::string current_ID( vecStation[ii].getStationID() );
		if (editingStack.count(current_ID)>0) {
			for (size_t jj=0; jj<editingStack[current_ID].size(); jj++) {
				editingStack[current_ID][jj]->editTimeSeries(vecStation);
			}
		}
	}
	
	//remove the stations that have been merged into other ones, if necessary
	for (size_t ii=0; ii<vecStation.size(); ii++) {
		const std::string stationID( IOUtils::strToUpper(vecStation[ii].stationID) );
		if (mergedFromIDs.count( stationID ) >  0) {
			std::swap( vecStation[ii], vecStation.back() );
			vecStation.pop_back();
			ii--; //in case we have multiple identical stations ID
		}
	}
}

void DataEditing::editTimeSeries(std::vector<METEO_SET>& vecMeteo)
{
	//TODO handle CREATE command
	//TODO there is a problem: to exclude some parameters from one
	//station and then merge with another is not guaranteed to work
	//in this new architecture...
	//TODO: handle checking for circular dependencies
	
	//check for "circular merge", build the list of stations that will be purged in the end 
	const std::set<std::string> mergedFromIDs( getMergedFromIDs() );
	
	//process widlcard commands first, knowing that '*' merges are prohibited
	if (editingStack.count("*")>0) {
		for (size_t jj=0; jj<editingStack["*"].size(); jj++)
			editingStack["*"][jj]->editTimeSeries(vecMeteo);
	}
	
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		if (vecMeteo[ii].empty()) continue;
		
		const std::string current_ID( vecMeteo[ii].front().getStationID() );
		if (editingStack.count(current_ID)>0) {
			for (size_t jj=0; jj<editingStack[current_ID].size(); jj++) {
				editingStack[current_ID][jj]->editTimeSeries(vecMeteo);
			}
		}
	}
	
	//remove the stations that have been merged into other ones, if necessary
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		if (vecMeteo[ii].empty())  continue;
		const std::string stationID( IOUtils::strToUpper(vecMeteo[ii][0].meta.stationID) );
		if (mergedFromIDs.count( stationID ) >  0) {
			std::swap( vecMeteo[ii], vecMeteo.back() );
			vecMeteo.pop_back();
			ii--; //in case we have multiple identical stations ID
		}
	}
	
	//remove trailing pure nodata MeteoData elements (if any)
	purgeTrailingNodata(vecMeteo);

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

const std::string DataEditing::toString() const
{
	std::ostringstream os;
	os << "<DataEditing>\n";

	os << dataCreator.toString();

	os << "</DataEditing>\n";
	return os.str();
}

} //end namespace
