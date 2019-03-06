/***********************************************************************************/
/*  Copyright 2019 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/plugins/Goes.h>
#include <meteoio/FileUtils.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <errno.h>
#include <cstdlib>

using namespace std;

namespace mio {
/**
 * @page goesio GoesIO
 * This plugin deals with data that has been transmitted through the
 * <A HREF="https://en.wikipedia.org/wiki/Geostationary_Operational_Environmental_Satellite">GOES</A> satellites
 * (see also https://www.goes-r.gov/resources/docs.html and https://www.rtl-sdr.com/tag/goes/). In order to reduce data
 * transfers, no headers are provided with the data and therefore all the metadata will have to be provided to the plugin.
 *
 * Currently, this plugin is only designed to read data in the GOES format.
 *
 * @section goes_format Format
 *
 * @section goes_units Units
 *
 * @section goes_keywords Keywords
 * This plugin uses the following keywords, all in the [Input] section:
 * - COORDSYS: coordinate system (see Coords);
 * - COORDPARAM: extra coordinates parameters (see Coords);
 * - TIMEZONE: timezone of the data;
 * - METEOPATH: directory where to read the data files from;
 * - GOES_NODATA: value used to represent nodata (default: -8190);
 * - GOES_DEBUG: should extra (ie very verbose) information be displayed? (default: false)
 */

//these are fixed by GOES
static const size_t nElems = 46;
static const size_t dataStartPos = 37;

GoesIO::GoesIO(const std::string& configfile) : vecFilenames(), fields_idx(), md_template(), meteopath(), coordin(), coordinparam(),
                                                in_TZ(0.), in_nodata(-8190.), stationID_idx(IOUtils::npos), year_idx(IOUtils::npos), jdn_idx(IOUtils::npos), hour_idx(IOUtils::npos), debug(false)
{
	parseInputOutputSection( Config(configfile) );
}

GoesIO::GoesIO(const Config& cfgreader) : vecFilenames(), fields_idx(), md_template(), meteopath(), coordin(), coordinparam(),
                                          in_TZ(0.), in_nodata(-8190.), stationID_idx(IOUtils::npos), year_idx(IOUtils::npos), jdn_idx(IOUtils::npos), hour_idx(IOUtils::npos), debug(false)
{
	parseInputOutputSection( cfgreader );
}

void GoesIO::parseInputOutputSection(const Config& cfg)
{
	cfg.getValue("TIME_ZONE", "Input", in_TZ);
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam);
	
	cfg.getValue("GOES_DEBUG", "Input", debug, IOUtils::nothrow);
	cfg.getValue("GOES_NODATA", "Input", in_nodata, IOUtils::nothrow);
	cfg.getValue("METEOPATH", "Input", meteopath);
	cfg.getValues("STATION", "Input", vecFilenames);
	std::vector<std::string> fields = cfg.get("FIELDS", "Input");
	if (fields.size()!=nElems) {
		ostringstream ss;
		ss << "Number of user-declared fields (" << fields.size() << ") don't match with GOES message ";
		ss << "number of fields (" << nElems << ")";
		throw InvalidArgumentException(ss.str(), AT);
	}
	parseFieldsSpecs(fields, md_template, fields_idx);
}

void GoesIO::readStationData(const Date& /*date*/, std::vector<StationData>& /*vecStation*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
	//implement a way to provide station metadata: location, timezone, units_offset & units factors, fields
}

void GoesIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                             std::vector< std::vector<MeteoData> >& vecMeteo)
{
	vecMeteo.clear();

	for (size_t ii=0; ii<vecFilenames.size(); ii++)
		readRawGoes( meteopath+"/"+vecFilenames[ii], dateStart, dateEnd, vecMeteo );
}

void GoesIO::readRawGoes(const std::string& file_and_path, const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo) const
{
	if (!FileUtils::validFileAndPath(file_and_path)) //Check whether filename is valid
		throw InvalidNameException(file_and_path, AT);

	errno = 0;
	std::ifstream fin(file_and_path.c_str(), ios::in);
	if (fin.fail()) {
		ostringstream ss;
		ss << "Error opening file \"" << file_and_path << "\" for reading, possible reason: " << strerror(errno);
		ss << " Please check file existence and permissions!";
		throw AccessException(ss.str(), AT);
	}

	std::map<std::string, size_t> mapStations;
	std::vector<float> raw_data(nElems);
	while (!fin.eof()){
		std::string line;
		getline(fin, line);
		if (line.length()<dataStartPos) continue;

		const std::string goes_id( line.substr(0, 8) ); //only first 8 characters of the station substring
		if (goes_id[0]!='8') continue; //invalid GOES ID, this must be an invalid line
		const std::string chain( line.substr(dataStartPos) );

		for (size_t ii=1; ii<=nElems; ii++) {
			raw_data[ii-1] = 0;
			const int A = chain[3*ii-2] & 15;
			const int B = chain[3*ii-1] & 63;
			const int C = chain[3*ii-0] & 63;

			if ((A*64+B) > 1008) {
				raw_data[ii-1] = static_cast<float>( (B-48)*64 + C + 9000 );
				continue;
			}

			float SF = ((A & 8) != 0)? -1 : 1;
			if ((A & 4) != 0) SF *= 0.01f;
			if ((A & 2) != 0) SF *= 0.1f;
			if ((A & 1) != 0) raw_data[ii-1] = 4096.;

			raw_data[ii-1] = (raw_data[ii-1] + static_cast<float>((B & 63)*64 + (C & 63))) * SF;
		}
		
		const Date dt( parseDate(raw_data) );
		if (dt.isUndef()) continue; //that was an invalid line
		if (dt<dateStart) continue;
		if (dt>dateEnd) {
			fin.close();
			return;
		}

		if (mapStations.count( goes_id )==0) { //create the station and get its index
			mapStations[ goes_id ] = vecMeteo.size();
			vecMeteo.push_back( std::vector<MeteoData>() );
		}
		const size_t st_idx = mapStations[ goes_id];

		const MeteoData md( parseDataLine("Goes::"+goes_id, dt, raw_data) );
		if (vecMeteo[st_idx].size() > 0 && md.date<=vecMeteo[st_idx].back().date) continue;
		vecMeteo[ st_idx ].push_back( md );
		
		if (debug) {
			std::cout << goes_id << " " << dt.toString(Date::ISO) << " -> ";
			for (size_t ii=0; ii<raw_data.size(); ii++) std::cout << "   " << raw_data[ii];
			std::cout << "\n";
		}
	}
	
	fin.close();
}

void GoesIO::parseFieldsSpecs(const std::vector<std::string>& fieldsNames, MeteoData &meteo_template, std::vector<size_t> &idx)
{
	idx.resize(fieldsNames.size(), IOUtils::npos);

	for (size_t ii=0; ii<fieldsNames.size(); ii++) {
		const std::string parname( IOUtils::strToUpper(fieldsNames[ii]) );
		if (parname=="SKIP") continue;

		if (parname=="STATIONID") {
			stationID_idx=ii;
			continue;
		}
		if (parname=="YEAR") {
			year_idx=ii;
			continue;
		}
		if (parname=="JDN") {
			jdn_idx=ii;
			continue;
		}
		if (parname=="HOUR") {
			hour_idx=ii;
			continue;
		}

		const size_t curr_idx = meteo_template.getParameterIndex(parname);
		if (curr_idx!=IOUtils::npos) idx[ii] = curr_idx;
		else idx[ii] = meteo_template.addParameter(parname);
	}
}

Date GoesIO::parseDate(const std::vector<float>& raw_data) const
{
	const double jdn = raw_data[ jdn_idx ] + raw_data[ hour_idx ]/24.;
	const int year = static_cast<int>(raw_data[ year_idx ]);
	return Date(year, jdn, in_TZ);
}

MeteoData GoesIO::parseDataLine(const std::string& stationName, const Date& dt, const std::vector<float>& raw_data) const
{
	MeteoData md( md_template );
	md.date.setDate(dt);
	md.meta.stationName = stationName;
	md.meta.stationID = IOUtils::toString( raw_data[stationID_idx] );

	for (size_t ii=0; ii<raw_data.size(); ii++) {
		const size_t idx = fields_idx[ii];
		if (idx==IOUtils::npos) continue;
		md( idx ) = IOUtils::standardizeNodata( raw_data[ii], in_nodata );
	}

	return md;
}

} //namespace
