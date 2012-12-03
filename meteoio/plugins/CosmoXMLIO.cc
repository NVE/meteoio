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
#include "CosmoXMLIO.h"
#include <meteoio/meteolaws/Atmosphere.h>
#include <sstream>

//To read a text file
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

namespace mio {
/**
 * @page cosmoxml COSMOXML
 * @section cosmoxml_format Format
 * This plugin reads the XML files as generated by the Cosmo system.
 * The files are outputted in Grib format and preprocessed by FieldExtra (MeteoSwiss)
 * to get XML files.
 * It requires libxml to compile and run.
 * When writing files, it creates one file per station nammed as {stationID}_{numerical date}.xml
 *
 * @section cosmoxml_units Units
 * The units are assumed to be the following:
 * - temperatures in celsius
 * - dew point in celsius (input)
 * - relative humidity in % (output)
 * - wind speed in m/s
 * - precipitations in mm/h
 * - radiation in W/m²
 * - maximal wind speed in m/s (not implemented yet...)
 *
 * @section cosmoxml_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS:  input coordinate system (see Coords) specified in the [Input] section
 * - METEOPATH: string containing the path to the xml files to be read, specified in the [Input] section
 * - METEOPATH: string containing the path where the XML files have to be written, specified in the [Output] section
 * - METEO:     specify COSMOXML for [Input] and/or [Output] section(s)
 *
 * Example:
 * @code
 * [Input]
 * COORDSYS	= CH1903
 * METEO	= COSMOXML
 * METEOPATH	= ./input/meteoXMLdata
 * [Output]
 * METEO	= COSMOXML
 * METEOPATH	= ./output/meteoXMLdata
 * @endcode
 */

const double CosmoXMLIO::in_tz = 0.; //Plugin specific timezone
const double CosmoXMLIO::out_tz = 0.; //Plugin specific time zone
const std::string CosmoXMLIO::dflt_extension = ".xml";

CosmoXMLIO::CosmoXMLIO(const std::string& configfile)
           : cfg(configfile), plugin_nodata(-999.), coordin(), coordinparam(), coordout(), coordoutparam()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

CosmoXMLIO::CosmoXMLIO(const Config& cfgreader)
           : cfg(cfgreader), plugin_nodata(-999.), coordin(), coordinparam(), coordout(), coordoutparam()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

CosmoXMLIO::~CosmoXMLIO() throw()
{

}

void CosmoXMLIO::read2DGrid(Grid2DObject& /*grid_out*/, const std::string& /*_name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void CosmoXMLIO::read2DGrid(Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void CosmoXMLIO::readDEM(DEMObject& /*dem_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void CosmoXMLIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void CosmoXMLIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void CosmoXMLIO::readStationData(const Date& station_date, std::vector<StationData>& vecStation)
{
	//Get all files from directory
	string meteopath, station_path;
 	cfg.getValue("METEOPATH", "Input", meteopath);
	if (meteopath == "")
		throw ConversionFailedException("Error while reading value for METEOPATH", AT);
	const string pattern = dflt_extension;
	list<string> dirlist;
	IOUtils::readDirectory(meteopath, dirlist, pattern);
	dirlist.sort();

	vecStation.clear(); //Initialize Station Vector

	for( list<string>::iterator itr = dirlist.begin(); itr != dirlist.end(); itr++ ) {
		//Loop over all stations in the meteopath directory
		station_path = meteopath + "/" + *itr;

		StationData sd;
		//Initialize variables
		double altitude=IOUtils::nodata, latitude=IOUtils::nodata, longitude=IOUtils::nodata;
		string station_name, station_ID;
		bool is_first=true;
		Date first_date=IOUtils::nodata, last_date=IOUtils::nodata, date_read=IOUtils::nodata;

		//Read station and meteo data
		xmlpp::TextReader reader(station_path);
		while(reader.read()) {
			if(reader.has_attributes()) {
				reader.move_to_first_attribute();
				const string key=reader.get_value();
				//StationData
				if(key=="identifier") station_name = getValue(reader);
				if(key=="station_abbreviation") station_ID = getValue(reader);
				if(key=="station.height") altitude = getDoubleValue(reader);
				if(key=="station.latitude") latitude = getDoubleValue(reader);
				if(key=="station.longitude") longitude = getDoubleValue(reader);
 				if(key=="missing_value_code") plugin_nodata = getDoubleValue(reader);
				//Test used to write station data only if the given date is in the station data interval
				if(key=="reference_ts") {
					if(!is_first) {
						//Check date
						date_read = getDateValue(reader);
						if (date_read < first_date) {
							first_date = date_read;
						}
						if (date_read > last_date) {
							last_date = date_read;
						}
					} else {
						//First date to be read -> is first and last date
						first_date = getDateValue(reader);
						last_date = first_date;
						is_first = false;
					}
				}
				//End of date test
			}
		}

		//Write station data if the given date is in the station data interval
		if ((station_date > first_date) && (station_date < last_date)) {
			sd.stationName = station_name;
			sd.stationID = station_ID;
			sd.position.setLatLon(latitude, longitude, altitude);
			vecStation.push_back(sd);	//Store results in vecStation
		}

		//Write station data for any date
// 		sd.stationName = station_name;
// 		sd.stationID = station_ID;
// 		sd.position.setLatLon(latitude, longitude, altitude);
// 		vecStation.push_back(sd);	//Store results in vecStation

	}
}

//Use parser to get value corresponding to a key
std::string CosmoXMLIO::getValue(xmlpp::TextReader& reader) {
	reader.move_to_element();
	reader.read();
	const string value = reader.get_value();
	reader.read();
	reader.read();
	return( value );
}

//Conversion from string to double of the value read by the parser
double CosmoXMLIO::getDoubleValue(xmlpp::TextReader& reader) {
	double value;
	const string Svalue = getValue(reader);
	if( IOUtils::convertString(value, Svalue) == false) {
		stringstream ss;
		ss << "Can not parse (double)" << Svalue;
		throw ConversionFailedException(ss.str(), AT);
	}
	return IOUtils::standardizeNodata(value, plugin_nodata);
}

//Conversion °C to °K
double CosmoXMLIO::c2k(const double& value) {
	if(value!=IOUtils::nodata)
		return C_TO_K(value);
	else
		return value;
}

//Conversion °K to °C
double CosmoXMLIO::k2c(const double& value) {
	if(value!=IOUtils::nodata)
		return K_TO_C(value);
	else
		return value;
}

//Conversion from string to Date of the value read by the parser
Date CosmoXMLIO::getDateValue(xmlpp::TextReader& reader) {
	Date date;
	const string Svalue = getValue(reader);
	if( IOUtils::convertString(date, Svalue, in_tz) == false) {
		stringstream ss;
		ss << "Can not parse (date)" << Svalue;
		throw ConversionFailedException(ss.str(), AT);
	}
	return date;
}

//Converts dew point to humidity and writes station position
void CosmoXMLIO::finishMeteo(const double& latitude, const double& longitude, const double& altitude,
                             double& dew_point, MeteoData& meteo) {
	// Set position
	if(altitude!=IOUtils::nodata && latitude!=IOUtils::nodata && longitude!=IOUtils::nodata) {
		meteo.meta.position.setLatLon(latitude, longitude, altitude);
	} else {
		throw InvalidFormatException("Meteo data found, but no position information", AT);
	}

	//Set RH
	if ((meteo(MeteoData::TA) != IOUtils::nodata) && (dew_point != IOUtils::nodata)) {
		meteo(MeteoData::RH) = Atmosphere::DewPointtoRh(dew_point, meteo(MeteoData::TA), TRUE);
	}

	//Reset dew_point
	dew_point = IOUtils::nodata;
}

//----------> Read Station and Meteo data, uses all stations in the "meteopath" directory <----------
void CosmoXMLIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                               std::vector< std::vector<MeteoData> >& vecMeteo,
                               const size_t&)
{
	//Get all files from directory
	string meteopath, station_path;
 	cfg.getValue("METEOPATH", "Input", meteopath);
	if (meteopath == "")
		throw ConversionFailedException("Error while reading value for METEOPATH", AT);
	const string pattern = "xml";
	list<string> dirlist;
	IOUtils::readDirectory(meteopath, dirlist, pattern);
	dirlist.sort();

	vecMeteo.clear();	//Initialize Meteo Vector

	vecMeteo.insert(vecMeteo.begin(), dirlist.size(), std::vector<MeteoData>());	//Allocation for the vectors
	list<string>::iterator itr;	//To loop in the stations list

	unsigned int ii=0;	//Declare and initialize counter (to know which station we are dealing with)
	for( itr = dirlist.begin(); itr != dirlist.end(); itr++ ) {	//Loop over all stations in the meteopath directory
		station_path = meteopath + "/" + *itr;

		MeteoData meteo;
		//Initialize variables
		double altitude=IOUtils::nodata, latitude=IOUtils::nodata, longitude=IOUtils::nodata;
		double dew_point=IOUtils::nodata;
		bool is_first=true;
		bool next_station = false;

		//Read station and meteo data
		xmlpp::TextReader reader(station_path);
		while(reader.read()) {
			if(reader.has_attributes()) {
				reader.move_to_first_attribute();
				const string key=reader.get_value();
				//StationData
				if(key=="identifier") meteo.meta.stationName = getValue(reader);
				if(key=="station_abbreviation") meteo.meta.stationID = getValue(reader);
				if(key=="station.height") altitude = getDoubleValue(reader);
				if(key=="station.latitude") latitude = getDoubleValue(reader);
				if(key=="station.longitude") longitude = getDoubleValue(reader);
	// 			if(key=="model_station_height") altitude = getDoubleValue(reader);
	// 			if(key=="model_station_latitude") latitude = getDoubleValue(reader);
	// 			if(key=="model_station_longitude") longitude = getDoubleValue(reader);
 				if(key=="missing_value_code") plugin_nodata = getDoubleValue(reader);

				//MeteoData
				if(key=="T_2M") meteo(MeteoData::TA) = c2k( getDoubleValue(reader) );
				if(key=="TD_2M") dew_point = c2k( getDoubleValue(reader) );
				if(key=="GLOB") meteo(MeteoData::ISWR) = getDoubleValue(reader);
				if(key=="TOT_PREC") meteo(MeteoData::HNW) = getDoubleValue(reader);
				if(key=="FF_10M") meteo(MeteoData::VW) = getDoubleValue(reader);
				if(key=="VMAX_10M") meteo(MeteoData::VW_MAX) = getDoubleValue(reader);
				if(key=="reference_ts") {
					if(!is_first) {
						//We can finish our object and push it.
						//Before we check if it is in the right time interval
						if (meteo.date < dateStart) {
							meteo.reset();
							meteo.date = getDateValue(reader);
							//Go to the next entry
							continue;
						}
						if (meteo.date > dateEnd) {
							meteo.reset();
							//Go to the next station
							next_station = true;
							break;
						}
						finishMeteo(latitude, longitude, altitude, dew_point, meteo);
						vecMeteo[ii].push_back( meteo );
					} else {
						//Nothing in memory to be written when we meet the first date entry
						is_first = false;
					}
					meteo.reset();
					meteo.date = getDateValue(reader);
				}
			}
		}
		if (next_station==true) continue;
		finishMeteo(latitude, longitude, altitude, dew_point, meteo);
		vecMeteo[ii].push_back( meteo );
		reader.close();
		ii++;
	}
}

//Writes header and description of station data
void CosmoXMLIO::writeHeader(std::stringstream& XMLdata)
{
	XMLdata << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	XMLdata << "<!DOCTYPE tree-table-xml SYSTEM \"tree-table-xml.dtd\">\n<tree-table-xml>\n";
	XMLdata << "<ttable id=\"MeteoIO->GribXML_save\">\n<row>\n<col id=\"datainformation-tables\">\n<ttable id=\"metadata\">\n";
	XMLdata << "<row>\n<col id=\"id\">identifier</col>\n<col id=\"element\">unique composed identifier</col>\n<col id=\"unit\">character</col>\n<col id=\"class\">java.lang.String</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">0</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">station_abbreviation</col>\n<col id=\"element\">station abbrevation</col>\n<col id=\"unit\">character</col>\n<col id=\"class\">java.lang.String</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">0</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">station_height</col>\n<col id=\"element\">station height</col>\n<col id=\"unit\">m</col>\n<col id=\"class\">java.lang.Double</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">0</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">station_latitude</col>\n<col id=\"element\">station latitude</col>\n<col id=\"unit\">degree</col>\n<col id=\"class\">java.lang.Double</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">4</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">station_longitude</col>\n<col id=\"element\">station longitude</col>\n<col id=\"unit\">degree</col>\n<col id=\"class\">java.lang.Double</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">4</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">model_station_height</col>\n<col id=\"element\">height of model grid point associated to station</col>\n<col id=\"unit\">m</col>\n<col id=\"class\">java.lang.Double</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">0</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">model_station_latitude</col>\n<col id=\"element\">latitude of model grid point associated to station</col>\n<col id=\"unit\">degree</col>\n<col id=\"class\">java.lang.Double</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">4</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">model_station_longitude</col>\n<col id=\"element\">longitude of model grid point associated to station</col>\n<col id=\"unit\">degree</col>\n<col id=\"class\">java.lang.Double</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">4</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">data_cat_number</col>\n<col id=\"element\">data cat number</col>\n<col id=\"unit\">decimal code</col>\n<col id=\"class\">java.lang.Double</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">0</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">data_cat_description</col>\n<col id=\"element\">data cat description</col>\n<col id=\"unit\">character</col>\n<col id=\"class\">java.lang.String</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">0</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">0</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">model_configuration</col>\n<col id=\"element\">model configuration</col>\n<col id=\"unit\">decimal code</col>\n<col id=\"class\">java.lang.Double</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">0</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">missing_value_code</col>\n<col id=\"element\">missing value code</col>\n<col id=\"unit\">code</col>\n<col id=\"class\">java.lang.Double</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">1</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">model_run_dt</col>\n<col id=\"element\">model run date (UTC)</col>\n<col id=\"unit\">YYYYMMDDHHmm</col>\n<col id=\"class\">java.lang.String</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">1</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "</ttable>\n<!-- metadata datainformation section -->\n</col>\n";
}

//Writes station data
void CosmoXMLIO::writeLocationHeader(const StationData& station, std::stringstream& XMLdata)
{
		XMLdata << "<col>\n<ttable id=\"data\">\n<row>" << "\n";
		XMLdata << "<col id=\"identifier\">" << station.getStationName() << "</col>\n";
		XMLdata << "<col id=\"station_abbreviation\">" << station.getStationID() << "</col>\n";

		//Same altitude, latitude and longitude data for station and model.
		XMLdata << "<col id=\"station.height\">" << station.position.getAltitude() << "</col>\n";
		XMLdata << "<col id=\"station.latitude\">" << station.position.getLat() << "</col>\n";
		XMLdata << "<col id=\"station.longitude\">" << station.position.getLon() << "</col>\n";
 		XMLdata << "<col id=\"model_station_height\">" << station.position.getAltitude() << "</col>\n";
		XMLdata << "<col id=\"model_station_latitude\">" << station.position.getLat() << "</col>\n";
		XMLdata << "<col id=\"model_station_longitude\">" << station.position.getLon() << "</col>\n";
// 		XMLdata << "<col id=\"data_cat_number\">" << "TO DO!" << "</col>\n";
// 		XMLdata << "<col id=\"data_cat_description\">" << "TO DO!" << "</col>\n";
// 		XMLdata << "<col id=\"model_configuration\">" << "TO DO!" << "</col>\n";
 		XMLdata << "<col id=\"missing_value_code\">" << IOUtils::nodata << "</col>\n";
// 		XMLdata << "<col id=\"model_run_dt\">" << "TO DO!" << "</col>\n";
}

//Writes description of meteo data
void CosmoXMLIO::writeMeteoDataDescription(std::stringstream& XMLdata)
{
	XMLdata << "</row>\n</ttable>\n<!-- data datainformation section -->\n</col></row>\n";
	XMLdata << "<row>\n<col id=\"values-tables\">\n<ttable id=\"metadata\">\n";
	XMLdata << "<row>\n<col id=\"id\">identifier</col>\n<col id=\"element\">unique composed identifier</col>\n<col id=\"unit\">character</col>\n<col id=\"class\">java.lang.String</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">0</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">reference_ts</col>\n<col id=\"element\">timestamp (UTC)</col>\n<col id=\"unit\">YYYYMMDDHHmm</col>\n<col id=\"class\">java.lang.String</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">0</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">T_2M</col>\n<col id=\"element\">Temperature 2m above ground, height corrected</col>\n<col id=\"unit\">degree Celsius</col>\n<col id=\"class\">java.lang.Double</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">1</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">TD_2M</col>\n<col id=\"element\">Dew point temperature 2 m above ground, height corrected</col>\n<col id=\"unit\">degree Celsius</col>\n<col id=\"class\">java.lang.Double</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">1</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">GLOB</col>\n<col id=\"element\">Global radiation</col>\n<col id=\"unit\">W/m^2</col>\n<col id=\"class\">java.lang.Double</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">1</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">TOT_PREC</col>\n<col id=\"element\">Total precipitation, hourly sum</col>\n<col id=\"unit\">mm</col>\n<col id=\"class\">java.lang.Double</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">1</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">FF_10M</col>\n<col id=\"element\">Wind speed 10m above ground</col>\n<col id=\"unit\">m/s</col>\n<col id=\"class\">java.lang.Double</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">1</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "<row>\n<col id=\"id\">VMAX_10M</col>\n<col id=\"element\">Wind gusts 10m above ground</col>\n<col id=\"unit\">m/s</col>\n<col id=\"class\">java.lang.Double</col>\n<col id=\"decimal-precision\" class=\"java.lang.Integer\">1</col>\n<col id=\"scale\" class=\"java.lang.Integer\">1</col>\n<col id=\"mandatory\" class=\"java.lang.Boolean\">1</col>\n</row>\n";
	XMLdata << "</ttable>\n<!-- metadata values section -->\n</col>\n";
}

//Writes meteo data
void CosmoXMLIO::writeMeteo(const std::vector<MeteoData>& vecMeteo, std::stringstream& XMLdata)
{
	XMLdata << "<col>\n<ttable id=\"data\">\n";
	for(unsigned int jj=0; jj<vecMeteo.size(); jj++) {
		XMLdata << "<row>\n";
		XMLdata << "<col id=\"identifier\">" << vecMeteo[jj].meta.getStationName() << "</col>\n";
		Date tmp_date(vecMeteo[jj].date);
		tmp_date.setTimeZone(out_tz);
		XMLdata << "<col id=\"reference_ts\">" << tmp_date.toString(Date::NUM) << "</col>\n";
		XMLdata << "<col id=\"T_2M\">" << k2c(vecMeteo[jj](MeteoData::TA)) << "</col>\n";

		if ((vecMeteo[jj](MeteoData::RH) == IOUtils::nodata) || (vecMeteo[jj](MeteoData::TA) == IOUtils::nodata)) {
			XMLdata << "<col id=\"TD_2M\">" << IOUtils::nodata << "</col>\n";
		} else {
			const double dew_point=Atmosphere::RhtoDewPoint(vecMeteo[jj](MeteoData::RH), vecMeteo[jj](MeteoData::TA), TRUE);
			XMLdata << "<col id=\"TD_2M\">" << k2c(dew_point) << "</col>\n";
		}

		XMLdata << "<col id=\"GLOB\">" << vecMeteo[jj](MeteoData::ISWR) << "</col>\n";
		XMLdata << "<col id=\"TOT_PREC\">" << vecMeteo[jj](MeteoData::HNW) << "</col>\n";
		XMLdata << "<col id=\"FF_10M\">" << vecMeteo[jj](MeteoData::VW) << "</col>\n";
		XMLdata << "<col id=\"VMAX_10M\">" << vecMeteo[jj](MeteoData::VW_MAX) << "</col>\n";
		XMLdata << "</row>\n";
	}
}

//Writes last lines of the file
void CosmoXMLIO::writeFooter(std::stringstream& XMLdata)
{
	XMLdata << "</ttable>\n<!-- data values section -->\n</col>\n";
	XMLdata << "</row>\n</ttable>\n</tree-table-xml>\n";
}

void CosmoXMLIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo,
                                const std::string&)
{
	ofstream fout;
	string filename="", meteopath_out = "", line;
	stringstream XMLdata;
	cfg.getValue("METEOPATH", "Output", meteopath_out);
	for (unsigned int ii=0; ii < vecMeteo.size(); ii++) {
		//Test existence of element 0
		if(vecMeteo[ii].size()==0) {
			cerr << "[E] Station " << ii+1 << " exists but contains no data! Skip writing it...\n";
			continue;
		}

		XMLdata.str("");	//Initialize XMLdata stringstream

		//Write the first part of the XML file (header and station data description)
		writeHeader(XMLdata);

		//Insert station data
		writeLocationHeader( vecMeteo[ii][0].meta, XMLdata );

		//Write the second part of the XML file (meteo data description)
		writeMeteoDataDescription(XMLdata);

		//Insert meteo data
		writeMeteo(vecMeteo[ii], XMLdata);

		//Write the last part of the XML file (after meteo data, only a few lines)
		writeFooter(XMLdata);

		//Save file
		std::string stat_id = vecMeteo[ii][0].meta.getStationID();
		if(stat_id=="") {
			stringstream ss;
			ss << "station" << ii;
			stat_id = ss.str();
		}
		filename = meteopath_out + "/" + stat_id + "_" + vecMeteo[ii][0].date.toString(Date::NUM) + ".xml";
		fout.open(filename.c_str());
		fout << XMLdata.str();
		fout.close();
	}
}

void CosmoXMLIO::readSpecialPoints(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void CosmoXMLIO::write2DGrid(const Grid2DObject& /*grid_in*/, const std::string& /*name*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void CosmoXMLIO::write2DGrid(const Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void CosmoXMLIO::cleanup() throw()
{

}

} //namespace
