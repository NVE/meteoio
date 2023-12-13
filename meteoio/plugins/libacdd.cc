// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2020 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/plugins/libacdd.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/dataClasses/CoordsAlgorithms.h>
#include <meteoio/IOUtils.h>
#include <meteoio/FileUtils.h>
#include <meteoio/dataClasses/Date.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/IOExceptions.h>

#include <fstream>
#include <cerrno>
#include <cstring>
#include <algorithm>


using namespace std;

namespace mio {

/**
* @brief Read all config keys from the selected section and apply some special processing for some keys.
* @details This is used as some sort of caching, only keeping the section of interest.
* @param[in] cfg Config object to read the configuration keys from
* @param[in] section section to read the keys from (all keys from the section will be read at this point)
* @param[in] allow_multi_line should multi-line content be supported?
*/
void ACDD::setUserConfig(const mio::Config& cfg, const std::string& section, const bool& allow_multi_line)
{
#ifdef WEBSERVICE
	ENVIDAT = cfg.get("ENVIDAT", "Output", false);
#endif
	for (size_t ii=0; ii<name.size(); ii++) {
		cfg.getValue(cfg_key[ii], section, value[ii], mio::IOUtils::nothrow);
		
		if (cfg_key[ii]=="ACDD_SUMMARY") { //overwrite with the content of summary_file if available
			const std::string summary_file = cfg.get("ACDD_SUMMARY_FILE", section, "");
			if (!summary_file.empty()) {
				std::string buffer;
				std::ifstream fin( summary_file.c_str() );
				if (fin.fail())
					throw mio::AccessException("Error opening ACDD_SUMMARY_FILE \""+summary_file+"\", possible reason: "+std::strerror(errno), AT);

				const char eoln = mio::FileUtils::getEoln(fin); //get the end of line character for the file
				try {
					do {
						std::string line;
						getline(fin, line, eoln); //read complete line
						if (allow_multi_line) buffer.append(line+"\n");
						else buffer.append(line+" ");
					} while (!fin.eof());
					fin.close();
				} catch (const std::exception&){
					if (fin.is_open()) fin.close();
					throw;
				}
				
				value[ii] = buffer;
			}
		}
	}
	setContactInfo();
#ifdef WEBSERVICE
	//special handling for webservice
	setSLFAsPublisher();
#endif
}

size_t ACDD::countCommas(const std::string& str)
{
	const size_t count = std::count_if( str.begin(), str.end(), []( const char& c ){return c ==',';});
	return count;
}

void ACDD::setContactInfo(const std::string& category, const std::vector<std::string>& vecKeys, const std::vector<std::string>& default_values, const bool& setDefaults)
{
	bool unset = true;
	size_t min_count = IOUtils::npos, max_count = IOUtils::npos;
	
	for (std::string att_name : vecKeys) {
		const std::string att_value( getAttribute(att_name) );
		if (!att_value.empty()) {
			unset = false;
			const size_t num_elems = countCommas( att_value );
			if (max_count==IOUtils::npos || num_elems>max_count) max_count = num_elems;
			if (min_count==IOUtils::npos || num_elems<min_count) min_count = num_elems;
		}
	}
	
	if (min_count!=max_count) throw mio::InvalidFormatException("Please configure the same number of fields for each comma-delimited ACDD fields of type '"+category+"'", AT);
	
	if (unset && setDefaults) {
		for (size_t ii=0; ii<vecKeys.size(); ii++) {
			addAttribute(vecKeys[ii], default_values[ii]);
		}
	}
}

void ACDD::setContactInfo()
{
	static const std::vector<std::string> creator_keys{"creator_name", "creator_email", "creator_institution", "creator_url", "creator_type"};
	static const std::vector<std::string> defaults_creator{mio::IOUtils::getLogName(), "", mio::IOUtils::getDomainName(), mio::IOUtils::getDomainName(), "person"};
	setContactInfo("CREATOR", creator_keys, defaults_creator, true);

	static const std::vector<std::string> publisher_keys{"publisher_name", "publisher_email", "publisher_url", "publisher_type"};
	static const std::vector<std::string> defaults_publisher{mio::IOUtils::getLogName(), "", mio::IOUtils::getDomainName(), "person"};
	setContactInfo("PUBLISHER", publisher_keys, defaults_publisher, true);

	static const std::vector<std::string> contributor_keys{"contributor_name", "contributor_role"};
	static const std::vector<std::string> defaults_contributor{"", "technical"};
	setContactInfo("CONTRIBUTOR", contributor_keys, defaults_contributor,false);
}

void ACDD::defaultInit()
{
	mio::Date now; 
	now.setFromSys();
	addAttribute("date_created", now.toString(mio::Date::ISO_DATE));
	addAttribute("institution", mio::IOUtils::getDomainName(), "ACDD_INSTITUTION");

	//defaults are set separately in setContactInfo()
	addAttribute("creator_name", "", "ACDD_CREATOR");
	addAttribute("creator_email", "", "ACDD_CREATOR_EMAIL");
	addAttribute("creator_institution", "", "ACDD_CREATOR_INSTITUTION");
	addAttribute("creator_url", "", "ACDD_CREATOR_URL");
	addAttribute("creator_type", "", "ACDD_CREATOR_TYPE");
	addAttribute("contributor_name", "", "ACDD_CONTRIBUTOR");
	addAttribute("contributor_role", "", "ACDD_CONTRIBUTOR_ROLE");
	addAttribute("publisher_name", "", "ACDD_PUBLISHER");
	addAttribute("publisher_email", "", "ACDD_PUBLISHER_EMAIL");
	addAttribute("publisher_url", "", "ACDD_PUBLISHER_URL");
	addAttribute("publisher_type", "", "ACDD_PUBLISHER_TYPE");

	addAttribute("source", "MeteoIO-" + mio::getLibVersion(true), "ACDD_SOURCE");
	addAttribute("history", now.toString(mio::Date::ISO_Z) + ", " + mio::IOUtils::getLogName() + "@" + mio::IOUtils::getHostName() + ", MeteoIO-" + mio::getLibVersion(true));
	addAttribute("keywords_vocabulary", "AGU Index Terms", "ACDD_KEYWORDS_VOCABULARY");
	addAttribute("keywords", "Cryosphere, Mass Balance, Energy Balance, Atmosphere, Land/atmosphere interactions, Climatology", "ACDD_KEYWORDS");
	addAttribute("title", "", "ACDD_TITLE");
	addAttribute("project", "", "ACDD_PROJECT");
	addAttribute("program", "", "ACDD_PROGRAM");
	addAttribute("id", "", "ACDD_ID");
	addAttribute("references", "", "ACDD_REFERENCES");
	addAttribute("naming_authority", "", "ACDD_NAMING_AUTHORITY");
	addAttribute("processing_level", "", "ACDD_PROCESSING_LEVEL");
	addAttribute("summary", "", "ACDD_SUMMARY"); //special handling, see setUserConfig()
	addAttribute("comment", "", "ACDD_COMMENT");
	addAttribute("acknowledgement", "", "ACDD_ACKNOWLEDGEMENT");
	addAttribute("metadata_link", "", "ACDD_METADATA_LINK");
	addAttribute("license", "", "ACDD_LICENSE");
	addAttribute("product_version", "1.0", "ACDD_PRODUCT_VERSION");
	addAttribute("activity_type", "", "ACDD_ACTIVITY_TYPE");
	addAttribute("operational_status", "", "ACDD_OPERATIONAL_STATUS");
	addAttribute("wmo__wsi", "", "WIGOS_ID");
}

/**
* @brief Add an attribute and its content to the internal list
* @details This allows to create or edit attributes. For the MERGE or APPEND modes, if the attribute name is not found, it will be created.
* @param[in] att_name attribute name
* @param[in] att_value attribute value
* @param[in] att_cfg_key associated configuration key (to read user provided values from a mio::Config object)
* @param[in] mode write mode: MERGE (currently empty values will be replaced by the given arguments), APPEND (the value content will be expanded by
* what is provided in att_value, separated by ", ", REPLACE (the current attribute will be fully replaced by the provided arguments)
*/
void ACDD::addAttribute(const std::string& att_name, const std::string& att_value, const std::string& att_cfg_key, Mode mode)
{
	if (att_name.empty())
		throw mio::InvalidFormatException("The attribute name must be provided", AT);
	
	if (mode==MERGE) {
		const size_t pos = find( att_name );
		if (pos==mio::IOUtils::npos) {
			mode = REPLACE;
		} else {
			if (!att_value.empty()) value[pos] = att_value;
			if (!att_cfg_key.empty()) cfg_key[pos] = att_cfg_key;
			return;
		}
	} else if (mode==APPEND) {
		const size_t pos = find( att_name );
		if (pos==mio::IOUtils::npos) {
			mode = REPLACE;
		} else {
			value[pos] = value[pos] + ", " + att_value;
			return;
		}
	}
	
	if (mode==REPLACE) {
		name.push_back( att_name );
		value.push_back( att_value );
		cfg_key.push_back( att_cfg_key );
		return;
	}
	
	//we should not have come here -> throw
	throw mio::InvalidFormatException("The specified write mode does not exists", AT);
}

void ACDD::addAttribute(const std::string& att_name, const double& att_value, const std::string& att_cfg_key, const Mode& mode)
{
	std::ostringstream os;
	os << att_value;
	addAttribute(att_name, os.str(), att_cfg_key, mode);
}

void ACDD::getAttribute(const size_t ii, std::string &att_name, std::string & att_value) const
{
	if (ii<name.size()) {
		att_name=name[ii];
		att_value=value[ii];
	} else {
		att_name="";
		att_value="";
	}
}

/**
* @brief Given an attribute name, return its associated value (or an empty string if it does not exists)
* @param[in] att_name attribute name to get the value for
* @return attribute value or empty string
*/
std::string ACDD::getAttribute(std::string &att_name) const
{
	for (size_t ii=0; ii<name.size(); ii++) {
		if (name[ii]==att_name) return value[ii];
	}
	
	return "";
}

/**
* @brief Given an attribute name, return its associated index (or IOUtils::npos if it does not exists)
* @param[in] search_name attribute name to get the index for
* @return attribute index or IOUtils::npos
*/
size_t ACDD::find(const std::string& search_name) const
{
	for (size_t ii=0; ii<name.size(); ii++) {
		if (name[ii]==search_name) return ii;
	}
	
	return mio::IOUtils::npos;
}

void ACDD::setGeometry(const mio::Grid2DObject& grid, const bool& isLatLon)
{
	mio::Coords urcorner(grid.llcorner);
	urcorner.moveByXY(static_cast<double>(grid.getNx())*grid.cellsize, static_cast<double>(grid.getNy())*grid.cellsize);
	
	std::string epsg_str = "4326";
	std::string geometry;
	if (isLatLon) {
		std::ostringstream ss;
		ss << std::fixed << std::setprecision(10) << grid.llcorner.getLon() << " " << grid.llcorner.getLat() << ", ";
		ss << urcorner.getLon() << " " << grid.llcorner.getLat() << ", ";
		ss << urcorner.getLon() << " " << urcorner.getLat() << ", ";
		ss << grid.llcorner.getLon() << " " << urcorner.getLat();
		geometry = ss.str();
	}else {
		std::ostringstream os;
		os << grid.llcorner.getEPSG();
		epsg_str = os.str();
		
		std::ostringstream ss;
		ss << std::fixed << std::setprecision(10) << grid.llcorner.getEasting() << " " << grid.llcorner.getNorthing() << ", ";
		ss << urcorner.getEasting() << " " << grid.llcorner.getNorthing() << ", ";
		ss << urcorner.getEasting() << " " << urcorner.getNorthing() << ", ";
		ss << grid.llcorner.getEasting() << " " << urcorner.getNorthing();
		geometry = ss.str();
	}
	
	addAttribute("geospatial_bounds_crs", "EPSG:"+epsg_str);
	addAttribute("geospatial_bounds", "Polygon (("+geometry+"))");

	addAttribute("geospatial_vertical_positive", "up");
	addAttribute("geospatial_vertical_units", "m");
	const double min_val = grid.getMin();
	if (min_val!=IOUtils::nodata) {
		//if there is a min_val, there is also a max_val
		addAttribute("geospatial_vertical_min", min_val);
		const double max_val = grid.getMax();
		addAttribute("geospatial_vertical_max", max_val);
	}


}

void ACDD::setGeometry(const std::vector< std::vector<mio::MeteoData> >& vecMeteo, const bool& isLatLon)
{
	if (vecMeteo.empty()) return;
	
	std::string multiPts;
	short int epsg = -1;
	double lat_min=90., lat_max=-90., lon_min=360., lon_max=-360.;
	double alt_min=std::numeric_limits<double>::max(), alt_max=-std::numeric_limits<double>::max();
	bool found = false;
	for (const std::vector<mio::MeteoData>& timeseries : vecMeteo) {
		if (timeseries.empty()) continue;

		//create the strings for the MultiPoint property
		std::ostringstream ss;
		if (isLatLon) {
			ss  << std::fixed << std::setprecision(10) << "(" << timeseries.front().meta.position.getLon() << " " << timeseries.front().meta.position.getLat() << ")";
		} else {
			ss  << std::fixed << std::setprecision(0) << "(" << timeseries.front().meta.position.getEasting() << " " << timeseries.front().meta.position.getNorthing() << ")";
		}
		if (epsg==-1) { //first valid point
			epsg = (isLatLon)? 4326 : timeseries.front().meta.position.getEPSG();
			multiPts = ss.str();
		} else {
			if (!isLatLon && epsg!=timeseries.front().meta.position.getEPSG()) epsg = 0; //we use 0 as a marker for non-consistent epsg between points
			multiPts += ", "+ss.str();
		}

		const double curr_lat = timeseries.front().meta.position.getLat();
		const double curr_lon = timeseries.front().meta.position.getLon();
		const double curr_alt = timeseries.front().meta.position.getAltitude();
		found = true;
		
		if (lat_min>curr_lat) lat_min = curr_lat;
		if (lat_max<curr_lat) lat_max = curr_lat;
		if (lon_min>curr_lon) lon_min = curr_lon;
		if (lon_max<curr_lon) lon_max = curr_lon;
		if (alt_min>curr_alt) alt_min = curr_alt;
		if (alt_max<curr_alt) alt_max = curr_alt;
	}
	if (!found) return;
	
	if (epsg>0) { //ie there is at least one valid point and all further points use the same epsg
		std::ostringstream os;
		os << epsg;
		addAttribute("geospatial_bounds_crs", "EPSG:"+os.str());
		addAttribute("geospatial_bounds", "MultiPoint ("+multiPts+")");
	}
	addAttribute("geospatial_lat_min", lat_min);
	addAttribute("geospatial_lat_max", lat_max);
	addAttribute("geospatial_lon_min", lon_min);
	addAttribute("geospatial_lon_max", lon_max);

	addAttribute("geospatial_vertical_positive", "up");
	addAttribute("geospatial_vertical_units", "m");
	if (alt_min!=IOUtils::nodata) {
		addAttribute("geospatial_vertical_min", alt_min);
		addAttribute("geospatial_vertical_max", alt_max);
	}
}

void ACDD::setGeometry(const std::vector< mio::Coords >& vecLocation, const bool& isLatLon)
{
	if (vecLocation.empty()) return;
	
	std::string multiPts;
	short int epsg = -1;
	double lat_min=90., lat_max=-90., lon_min=360., lon_max=-360.;
	double alt_min=std::numeric_limits<double>::max(), alt_max=-std::numeric_limits<double>::max();
	
	for (const mio::Coords& location : vecLocation) {
		//create the strings for the MultiPoint property
		std::ostringstream ss;
		if (isLatLon) {
			ss  << std::fixed << std::setprecision(10) << "(" << location.getLon() << " " << location.getLat() << ")";
		} else {
			ss  << std::fixed << std::setprecision(0) << "(" << location.getEasting() << " " << location.getNorthing() << ")";
		}
		if (epsg==-1) { //first valid point
			epsg = (isLatLon)? 4326 : location.getEPSG();
			multiPts = ss.str();
		} else {
			if (!isLatLon && epsg!=location.getEPSG()) epsg = 0; //we use 0 as a marker for non-consistent epsg between points
			multiPts += ", "+ss.str();
		}

		const double curr_lat = location.getLat();
		const double curr_lon = location.getLon();
		const double curr_alt = location.getAltitude();
		
		if (lat_min>curr_lat) lat_min = curr_lat;
		if (lat_max<curr_lat) lat_max = curr_lat;
		if (lon_min>curr_lon) lon_min = curr_lon;
		if (lon_max<curr_lon) lon_max = curr_lon;
		if (alt_min>curr_alt) alt_min = curr_alt;
		if (alt_max<curr_alt) alt_max = curr_alt;
	}
	
	if (epsg>0) { //ie there is at least one valid point and all further points use the same epsg
		std::ostringstream os;
		os << epsg;
		addAttribute("geospatial_bounds_crs", "EPSG:"+os.str());
		
		const bool singlePoint = (lat_min==lat_max && lon_min==lon_max);
		if (singlePoint)
			addAttribute("geospatial_bounds", "Point ("+multiPts+")");
		else
			addAttribute("geospatial_bounds", "MultiPoint ("+multiPts+")");
	}
	
	addAttribute("geospatial_lat_min", lat_min);
	addAttribute("geospatial_lat_max", lat_max);
	addAttribute("geospatial_lon_min", lon_min);
	addAttribute("geospatial_lon_max", lon_max);
	
	addAttribute("geospatial_vertical_positive", "up");
	addAttribute("geospatial_vertical_units", "m");
	if (alt_min!=IOUtils::nodata) {
		addAttribute("geospatial_vertical_min", alt_min);
		addAttribute("geospatial_vertical_max", alt_max);
	}
}

void ACDD::setGeometry(const mio::Coords& location, const bool& isLatLon)
{
	std::string epsg_str = "4326";
	std::string geometry;
	if (isLatLon) {
		std::ostringstream ss;
		ss << std::fixed << std::setprecision(10) << location.getLon() << " " << location.getLat();
		geometry = ss.str();
	}else {
		std::ostringstream os;
		os << location.getEPSG();
		epsg_str = os.str();
		std::ostringstream ss;
		ss << std::fixed << std::setprecision(0) << location.getEasting() << " " << location.getNorthing();
	}
	addAttribute("geospatial_bounds_crs", "EPSG:"+epsg_str);
	addAttribute("geospatial_bounds", "Point ("+geometry+")");
	addAttribute("geospatial_lat_min", location.getLat());
	addAttribute("geospatial_lat_max", location.getLat());
	addAttribute("geospatial_lon_min", location.getLon());
	addAttribute("geospatial_lon_max", location.getLon());
	
	addAttribute("geospatial_vertical_positive", "up");
	addAttribute("geospatial_vertical_units", "m");
	addAttribute("geospatial_vertical_min", location.getAltitude());
	addAttribute("geospatial_vertical_max", location.getAltitude());
}

void ACDD::setTimeCoverage(const std::vector< std::vector<mio::MeteoData> >& vecMeteo)
{
	if (vecMeteo.empty()) return;
	
	mio::Date set_start( vecMeteo[0].front().date );
	mio::Date set_end( vecMeteo[0].back().date );
	int sampling_period = -1;
	for (const std::vector<mio::MeteoData>& timeseries : vecMeteo) { //we must redo station 0 in order to get sampling_period
		if (timeseries.empty()) continue;
		const mio::Date curr_start( timeseries.front().date );
		const mio::Date curr_end( timeseries.back().date );
		if (set_start>curr_start) set_start = curr_start;
		if (set_end<curr_end) set_end = curr_end;
		
		const size_t npts = timeseries.size();
		if (npts>1) {
			const int curr_sampling = static_cast<int>( (curr_end.getJulian() - curr_start.getJulian()) / static_cast<double>(npts-1) * 24.*3600. + .5);
			if (sampling_period<=0 || sampling_period>curr_sampling) sampling_period = curr_sampling;
		}
	}
	addAttribute( "time_coverage_start", set_start.toString(mio::Date::ISO_TZ));
	addAttribute("time_coverage_end", set_end.toString(mio::Date::ISO_TZ));
	
	if (sampling_period>0) {
		std::ostringstream os;
		os << "P" << sampling_period << "S"; //ISO8601 duration format
		addAttribute("time_coverage_resolution", os.str());
	}
}

void ACDD::setTimeCoverage(const std::vector<mio::MeteoData>& vecMeteo)
{
	if (vecMeteo.empty()) return;
	
	const mio::Date set_start( vecMeteo.front().date );
	const mio::Date set_end( vecMeteo.back().date );
	addAttribute( "time_coverage_start", set_start.toString(mio::Date::ISO_TZ));
	addAttribute("time_coverage_end", set_end.toString(mio::Date::ISO_TZ));
	
	const size_t npts = vecMeteo.size();
	if (npts>1) {
		const int sampling_period = static_cast<int>( (set_end.getJulian() - set_start.getJulian()) / static_cast<double>(npts-1) * 24.*3600. + .5);
		std::ostringstream os;
		os << "P" << sampling_period << "S"; //ISO8601 duration format
		addAttribute("time_coverage_resolution", os.str());
	}
}

void ACDD::setTimeCoverage(const std::vector<std::string>& vec_timestamp, const double& TZ)
{
	if (vec_timestamp.empty()) return;
	
	mio::Date set_start, set_end;
	mio::IOUtils::convertString(set_start, vec_timestamp.front(), TZ);
	mio::IOUtils::convertString(set_end, vec_timestamp.back(), TZ);
	addAttribute( "time_coverage_start", set_start.toString(mio::Date::ISO_TZ));
	addAttribute("time_coverage_end", set_end.toString(mio::Date::ISO_TZ));
	
	const size_t npts = vec_timestamp.size();
	if (npts>1) {
		const int sampling_period = static_cast<int>( (set_end.getJulian() - set_start.getJulian()) / static_cast<double>(npts-1) * 24.*3600. + .5);
		std::ostringstream os;
		os << "P" << sampling_period << "S"; //ISO8601 duration format
		addAttribute("time_coverage_resolution", os.str());
	}
}

void ACDD::setSLFAsPublisher()
{
	value[find("publisher_name")] = "WSL Institute for Snow and Avalanche Research SLF";
	value[find("publisher_type")] = "institution";
	if (!ENVIDAT) {
		value[find("publisher_email")] = "bavay@slf.ch, patrick.leibersperger@slf.ch";
		value[find("publisher_url")] = "https://service-meteoio.slf.ch";	}
	else {
		value[find("publisher_email")] = "";
		value[find("publisher_url")] = "https://www.envidat.ch";
	}
}

bool ACDD::ENVIDAT = false;

} //namespace

