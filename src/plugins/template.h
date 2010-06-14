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
#ifndef __TEMPLATE_H__
#define __TEMPLATE_H__

#include "IOInterface.h"
#include "Configreader.h"

#include <string>

namespace mio {

/**
 * @class TEMPLATE
 * @brief This (empty) class is to be used as a template for developing new plugins
 *
 * @author Mathias Bavay
 * @date   2010-06-14
 */
class TEMPLATE : public IOInterface {
	public:
		TEMPLATE(void (*delObj)(void*), const std::string& filename);

		TEMPLATE(const std::string& configfile);
		TEMPLATE(const ARPSIO&);
		TEMPLATE(const ConfigReader& cfgreader);
		~TEMPLATE() throw();

		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");

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
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& filename);

		void read3DGrid(Grid3DObject& grid_out, const std::string& _name);

	private:
		void cleanup() throw();

		ConfigReader cfg;
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
};

} //namespace
#endif
