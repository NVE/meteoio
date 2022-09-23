// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2022 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/plugins/SyntheticIO.h>

#include <regex>

using namespace std;

namespace mio {
/**
 * @page synthio SynthIO
 * This plugin is quite special since it does not read any data but generates synthetic data. It is designed for
 * numerical experiments where controlled conditions are applied to the numerical setup.
 * 
 * @section synthio_keywords Keywords
 * This plugin uses the following keywords, all in the [Input] section:
 * - controlling the timestamps generation:
 *     - TIME_ZONE: the time zone of any dates that are provided; 
 *     - SYNTH_START: when to start generating timestamps (optional, by default it generates timestamps for any requested date);
 *     - SYNTH_END: when to stop generating timestamps (optional, by default it generates timestamps for any requested date);
 *     - SYNTH_SAMPLING: sampling rate in seconds, starting with SYNTH_START or the requested date if no SYNTH_START is provided (mandatory);
 * - providing the stations' metadata:
 *     - COORDSYS: coordinate system (see Coords);
 *     - COORDPARAM: extra coordinates parameters (see Coords);
 *     - STATION#: coordinates of the station (mandatory, see \link Coords::Coords(const std::string& in_coordinatesystem, const std::string& in_parameters, std::string coord_spec) Coords()\endlink for the syntax);
 *     - ID#: the (short) station id to use (optional but recommended, default is "ID_#")
 *     - NAME#: a descriptive station name to use (optional, default is "STATION_#");
 * 
 * @section synthio_examples
 * Example of use:
 * @code
 * [INPUT]
 * COORDSYS = CH1903
 * TIME_ZONE = 1.00
 * 
 * METEO = SYNTH
 * SYNTH_SAMPLING = 1800
 * SYNTH_START = 2022-09-02T12:40
 * STATION = latlon (46.8, 9.81, 1500)
 * @endcode
 */

SynthIO::SynthIO(const std::string& configfile) : cfg(configfile), mapSynthGenerators(), vecStations(), dt_start(), dt_end(), dt_step(0.), TZ(0.)
{
	init();
}

SynthIO::SynthIO(const Config& cfgreader) : cfg(cfgreader), mapSynthGenerators(), vecStations(), dt_start(), dt_end(), dt_step(0.), TZ(0.)
{
	init();
}

void SynthIO::init()
{
	//read start / end time as well as sampling rate
	cfg.get("TIME_ZONE", "INPUT", TZ);
	cfg.getValue("SYNTH_SAMPLING", "INPUT", dt_step);
	const std::string dt_start_spec = cfg.get("SYNTH_START", "INPUT", "");
	if (!dt_start_spec.empty() && !IOUtils::convertString(dt_start, dt_start_spec, TZ))
		throw InvalidFormatException("Could not process start date "+dt_start_spec+" for the SYNTH plugin", AT);
	const std::string dt_end_spec = cfg.get("SYNTH_END", "INPUT", "");
	if (!dt_end_spec.empty() && !IOUtils::convertString(dt_end, dt_end_spec, TZ))
		throw InvalidFormatException("Could not process start date "+dt_end_spec+" for the SYNTH plugin", AT);
	
	//read the stations' basic metadata
	std::string coordin, coordinparam; //projection parameters
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam);
	std::vector<std::string> vecIDs, vecNames;
	
	const std::vector< std::pair<std::string, std::string> > coords_specs( cfg.getValuesRegex("STATION[0-9]+", "INPUT") );
	//cfg.getValues("STATION", "INPUT", coords_specs);
	cfg.getValues("ID", "INPUT", vecIDs);
	cfg.getValues("NAME", "INPUT", vecNames);
	const bool has_ids = !vecIDs.empty();
	const bool has_names = !vecNames.empty();
	
	//check for consistency
	if (coords_specs.empty()) throw InvalidArgumentException("Please provide at least one STATION# key containing coordinates for the SYNTH plugin", AT);
	if (has_ids && vecIDs.size()!=coords_specs.size()) throw InvalidArgumentException("Please either provide no IDS or the exact same number as STATION for the plugin SYNTH", AT);
	if (has_names && vecNames.size()!=coords_specs.size()) throw InvalidArgumentException("Please either provide no NAMES or the exact same number as STATION for the plugin SYNTH", AT);
	
	//build the stations' metadata
	for (size_t ii=0; ii<coords_specs.size(); ii++) {
		const Coords loc(coordin, coordinparam, coords_specs[ii].second);
		const std::string id = (has_ids)? vecIDs[ii] : "ID_"+IOUtils::toString( ii+1 );
		const std::string name = (has_names)? vecNames[ii] : "STATION_"+IOUtils::toString( ii+1 );
		const StationData sd(loc, id, name);
		
		vecStations.push_back( sd );
		mapSynthGenerators[ id ] = getSynthGenerators( coords_specs[ii].first );
	}
}

//get all necessary generators for the current station ID identified by its STATION# user-provided key, for all declared MeteoParameters
std::map< std::string, SynthGenerator* > SynthIO::getSynthGenerators(const std::string& stationRoot) const
{
	const std::regex parname_regex(stationRoot+"::([^:]+)::([^:]+)"); //extract the meteo parameter name and its subkey
	const std::vector<std::string> vec_keys( cfg.getKeys(stationRoot, "Input") );
	//std::cout << "stationRoot=" << stationRoot << " vec_keys.size()=" << vec_keys.size() << "\n";
	
	std::map< std::string, std::vector< std::pair<std::string, std::string> > > mapArgs;
	std::map< std::string, std::string > mapTypes;
	
	for (auto& key : vec_keys) {
		std::smatch index_matches;
		if (std::regex_match(key, index_matches, parname_regex)) { //retrieve the parname index
			const std::string parname( index_matches.str(1) );
			const std::string subkey( index_matches.str(2) );
			const std::string value( cfg.get(key, "INPUT", "") );
			
			if (subkey=="TYPE") {
				mapTypes[parname] = value;
			} else {
				mapArgs[parname].push_back( std::make_pair(subkey, value) );
			}
		}
	}
	
	std::map< std::string, SynthGenerator* > resultsMap;
	for (const auto& item : mapArgs) {
		const std::string parname( item.first );
		resultsMap[parname] = SynthFactory::getSynth( mapTypes[parname], stationRoot, parname, item.second, TZ);
	}
	
	return resultsMap;
}

void SynthIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                             std::vector< std::vector<MeteoData> >& vecMeteo)
{
	static const double days_to_sec = 24.*3600.;
	const double dt_step_days = dt_step / days_to_sec;
	vecMeteo.clear();
	
	//check if there is anything to deliver or not
	if (!dt_start.isUndef() && dt_start>dateEnd) return;	//requested period is before the user-defined range
	if (!dt_end.isUndef() && dt_end<dateStart) return;	//requested period is after the user-defined range
	
	const Date true_end = (dt_end.isUndef())? dateEnd : std::min(dateEnd, dt_end);
	Date true_start( dateStart );
	
	const double julian_start = (dt_start.isUndef())? dateStart.getJulian(true) : std::max(dateStart.getJulian(true), dt_start.getJulian(true));
	if (!dt_start.isUndef()) {
		if (dt_start<dateStart) {
			const int nr_steps_to_start = static_cast<int>( std::ceil( (julian_start - dt_start.getJulian(true)) / dt_step_days ) );	//watch out, dt_step is in seconds
			const double julian_true_start = dt_start.getJulian(true) + nr_steps_to_start * dt_step_days;
			true_start.setDate( julian_true_start, 0. ); //the julian calculations were done in gmt
		} else if (dt_start>dateStart) {
			true_start.setDate( dt_start );
		}
	}
	
	for (auto &station : vecStations) {
		std::map< std::string, SynthGenerator* > *station_synth = &mapSynthGenerators[ station.getStationID() ];
		
		std::vector<MeteoData> vecM;
		for (Date dt=true_start; dt<=true_end; dt+=dt_step_days) {
			MeteoData md( dt, station );
			
			for (const auto& item : (*station_synth) ) {
				md( item.first ) = item.second->generate( dt );
			}
			
			vecM.push_back( md );
		}
		vecMeteo.push_back( vecM );
	}
}


///////////////////////////////////////////////////////
//below the constructors (and argument parsing) and generate()  methods for all SynthGenerator

CST_Synth::CST_Synth(const std::string& station, const std::string& parname, const std::vector< std::pair<std::string, std::string> >& vecArgs) 
          : value(IOUtils::nodata)
{
	const std::string where( "SYNTH::CST, " + station + "::" + parname );
	bool has_value = false;
	
	//parse the arguments
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="VALUE") {
			IOUtils::parseArg(vecArgs[ii], where, value);
			has_value=true;
		}
	}
	
	if (!has_value) throw InvalidArgumentException("Please provide the VALUE argument for the "+where, AT);
}

double CST_Synth::generate(const Date& /*dt*/) const
{
	return value;
}

STEP_Synth::STEP_Synth(const std::string& station, const std::string& parname, const std::vector< std::pair<std::string, std::string> >& vecArgs, const double& TZ) 
          : dt_step(), value_before(IOUtils::nodata), value_after(IOUtils::nodata)
{
	const std::string where( "SYNTH STEP, " + station + "::" + parname );
	bool has_step_date = false;
	bool has_value_before = false, has_value_after = false;
	
	//parse the arguments
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="STEP_DATE") {
			if (!IOUtils::convertString(dt_step, vecArgs[ii].second, TZ))
				throw InvalidArgumentException("Can not parse argument '"+vecArgs[ii].first+"' for " + where, AT);
			has_step_date=true;
		} else if (vecArgs[ii].first=="VALUE_BEFORE") {
			IOUtils::parseArg(vecArgs[ii], where, value_before);
			has_value_before=true;
		} else if (vecArgs[ii].first=="VALUE_AFTER") {
			IOUtils::parseArg(vecArgs[ii], where, value_after);
			has_value_after=true;
		}
	}
	
	if (!has_value_before) throw InvalidArgumentException("Please provide the VALUE_BEFORE argument for the "+where, AT);
	if (!has_value_after) throw InvalidArgumentException("Please provide the VALUE_AFTER argument for the "+where, AT);
	if (!has_step_date) throw InvalidArgumentException("Please provide the STEP_DATE argument for the "+where, AT);
}

double STEP_Synth::generate(const Date& dt) const
{
	if (dt<dt_step)
		return value_before;
	else
		return value_after;
}

SynthGenerator* SynthFactory::getSynth(std::string type, const std::string& station, const std::string& parname, const std::vector< std::pair<std::string, std::string> >& vecArgs, const double& TZ)
{
	IOUtils::toUpper( type );
	
	if (type == "CST"){
		return new CST_Synth(station, parname, vecArgs);
	} else if (type == "STEP"){
		return new STEP_Synth(station, parname, vecArgs, TZ);
	} /*else if (type == "RECTANGLE"){
		return new RECT_Synth(station, parname, vecArgs);
	}*/ else {
		throw IOException("The Synthetizer '"+type+"' does not exist!" , AT);
	}
}

} //namespace
