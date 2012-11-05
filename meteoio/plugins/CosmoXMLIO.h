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
#ifndef __CosmoXMLIO_H__
#define __CosmoXMLIO_H__

#include <meteoio/IOInterface.h>
#include <meteoio/Config.h>
#include <meteoio/IOUtils.h>
#include <meteoio/Coords.h>
#include <meteoio/IOExceptions.h>

#include <libxml++/libxml++.h>
#include <libxml++/parsers/textreader.h>

#include <string>
#include <sstream>
#include <iostream>

namespace mio {

/**
* @class CosmoXMLIO
* @brief Reading of FieldExtra XML meteorological data.
* This is the plugin for reading the XML data genereated by FieldExtra, the post-processor
* of the MeteoSwiss COSMO meteorological model.
*
* @author Marc Diebold
* @date   January/February 2011
*/
class CosmoXMLIO : public IOInterface {
	public:
		CosmoXMLIO(const std::string& configfile);
		CosmoXMLIO(const CosmoXMLIO&);
		CosmoXMLIO(const Config& cfgreader);
		~CosmoXMLIO() throw();

		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);

		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo, const size_t& stationindex=IOUtils::npos);

		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string& name="");

		virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
		virtual void readSpecialPoints(std::vector<Coords>& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& filename);
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);

	private:
		void cleanup() throw();
		std::string getValue(xmlpp::TextReader& reader);
		double getDoubleValue(xmlpp::TextReader& reader);
		Date getDateValue(xmlpp::TextReader& reader);
		double c2k(const double& value);
		double k2c(const double& value);

		void finishMeteo(const double& latitude, const double& longitude, const double& altitude,
					double& dew_point, MeteoData& meteo);
		void writeHeader(std::stringstream& XMLdata);	//Write the first lines of the XML output file
		void writeLocationHeader(const StationData& station, std::stringstream& XMLdata);
		void writeMeteoDataDescription(std::stringstream& XMLdata);	//Write the middle of the XML output file
		void writeMeteo(const std::vector<MeteoData>& vecMeteo, std::stringstream& XMLdata);
		void writeFooter(std::stringstream& XMLdata);	//Write the last lines of the XML output file

		const Config& cfg;
		static const std::string dflt_extension;
		double plugin_nodata; //plugin specific no data value
		static const double in_tz, out_tz; //plugin specific time zones
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
};

} //namespace
#endif
