/***********************************************************************************/
/*  Copyright 2009 EPFL                                                            */
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
#ifndef __GEOTOPIO_H__
#define __GEOTOPIO_H__

#include "MeteoIO.h"

#include <string>
#include <sstream>
#include <iostream>

namespace mio {

/**
 * @class GeotopIO
 * @brief This class enables the access meteo data in legacy Geotop format
 *
 * @author Thomas Egger
 * @date   2009-07-02
 */
class GeotopIO : public IOInterface {
	public:
		//virtual GeotopIO* clone() const;

		GeotopIO(void (*delObj)(void*), const std::string& filename);

		GeotopIO(const std::string& configfile);
		GeotopIO(const GeotopIO&);
		GeotopIO(const ConfigReader&);
		~GeotopIO() throw();

		virtual void read2DGrid(Grid2DObject& dem_out, const std::string& parameter="");

		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  const unsigned int& stationindex=IOUtils::npos);

		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, 
							   const std::vector< std::vector<StationData> >& vecStation,
							   const std::string& name="");

		virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
		virtual void readSpecialPoints(std::vector<Coords>& pts);

		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& name);
		void read2DMeteo(const Date&, std::vector<MeteoData>&, std::vector<StationData>&); ///<No buffering

	private:
		void initParamNames(std::map<std::string, unsigned int>& mapParam);
		void readMetaData(std::vector<StationData>& vecStation, std::vector<std::string>& vecColumnNames,
					   const std::string& metafile);
		void makeColumnMap(const std::vector<std::string>& tmpvec, 
					    const std::vector<std::string>& vecColumnNames, 
					    std::map<std::string, unsigned int>& mapHeader);
		void convertUnits(MeteoData& meteo);
		void convertUnitsBack(MeteoData& meteo);
		void cleanup() throw();
		void parseDate(const std::string& datestring, const std::string& fileandline, Date& _date);

		ConfigReader cfg;
		std::ifstream fin; //Input file streams
		std::ofstream fout; //Output file streams
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
};

} //end namespace

#endif
