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
#ifndef __BORMAIO_H__
#define __BORMAIO_H__

#include "IOInterface.h"
#include "ConfigReader.h"
#include "IOUtils.h"

#include <string>
#include <sstream>
#include <libxml++/libxml++.h>
#include <iostream>

#include "IOExceptions.h"
#include "DynamicLibrary.h"


/**
 * @class BormaIO
 * @brief This class enables the access meteo data in Borma's XML format
 *
 * @author Thomas Egger
 * @date   2008-11-20
 */
class BormaIO : public IOInterface {
	public:
		//virtual BormaIO* clone() const;

		BormaIO(void (*delObj)(void*), const std::string& filename);

		BormaIO(const std::string& configfile);
		BormaIO(const BormaIO&);
		BormaIO(const ConfigReader&);
		~BormaIO() throw();

		virtual void read2DGrid(Grid2DObject& dem_out, const std::string& parameter="");

		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);

		virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  const unsigned int& stationindex=IOUtils::npos);

		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, 
							   const std::vector< std::vector<StationData> >& vecStation,
							   const std::string& name="");

		virtual void readAssimilationData(const Date_IO&, Grid2DObject& da_out);
		virtual void readSpecialPoints(std::vector<Coords>& pts);

		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& name);

		void read2DMeteo(const Date_IO&, std::vector<MeteoData>&); ///< No buffering
		void read2DMeteo(const Date_IO&, std::vector<MeteoData>&, std::vector<StationData>&); ///<No buffering

	private:
		void convertUnits(MeteoData& meteo);
		void checkForMeteoFiles(const std::string& xmlpath, const std::string& stationname, const Date_IO& date_in,
						    std::string& filename_out, Date_IO& date_out);
		void xmlParseStringToDouble(const std::string& str_in, double& d_out, const std::string& parname);
		std::string xmlGetNodeContent(xmlpp::Node* pNode, const std::string& nodename);
		void xmlExtractData(const std::string& filename, const Date_IO& date_in, MeteoData& md, StationData& sd);
		std::string xmlGetNodeName(xmlpp::Node* pNode);
		xmlpp::Node* xmlGetNode(xmlpp::Node* parentNode, const std::string& nodename);
		void stringToDate_IO(const std::string& tmp, Date_IO& date_out) const;
		bool validFilename(const std::string& tmp) const;
		void cleanup() throw();
		void getFiles(const std::string& stationsname, const Date_IO& start_date, const Date_IO& end_date,
				std::vector<std::string>& vecFiles, std::vector<Date_IO>& vecDate_IO);
		void readStationNames(void);
		bool bufferData(const Date_IO& dateStart, const Date_IO& dateEnd, 
					 std::vector< std::vector<MeteoData> >& vecMeteo, 
					 std::vector< std::vector<StationData> >& vecStation, 
					 const unsigned int& stationnr);
		void getProjectionParameters();

		ConfigReader cfg;
		std::ifstream fin; //Input file streams
		std::vector<std::string> vecStationName;
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		std::string coordsys, coordparam; //projection parameters
};

#endif
