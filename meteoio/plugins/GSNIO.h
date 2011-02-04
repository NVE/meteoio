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
#ifndef __GSNIO_H__
#define __GSNIO_H__

#include "gsn/soapA3DWebServiceSoap12BindingProxy.h"
#include "gsn/A3DWebServiceSoap12Binding.nsmap"

#include <meteoio/Config.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOUtils.h>
#include <meteoio/Coords.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/DynamicLibrary.h>

#include <string>
#include <sstream>
#include <iostream>
#include <vector>

namespace mio {

/**
 * @class GSNIO
 * @brief This class enables the access to the GSN web service
 *
 * @ingroup plugins
 * @author Thomas Egger
 * @date   2009-09-25
 */

class GSNIO : public IOInterface {
	public:
		GSNIO(void (*delObj)(void*), const Config& i_cfg);

		GSNIO(const std::string& configfile);
		GSNIO(const GSNIO&);
		GSNIO(const Config&);
		~GSNIO() throw();

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  const unsigned int& stationindex=IOUtils::npos);

		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, 
							   const std::string& name="");

		virtual void read2DGrid(Grid2DObject& dem_out, const std::string& parameter="");
		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);
		virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
		virtual void readSpecialPoints(std::vector<Coords>& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& name);

	private:
		void parseString(const std::string& in_string, std::vector<std::string>& vecString, MeteoData& md);
		void convertStringToDouble(double& d, const std::string& in_string, const std::string& in_parname);
		void convertUnits(MeteoData& meteo);
		void initGSNConnection();
		void readStationNames();
		void readStationMetaData(StationData& sd, const unsigned int& stationindex);
		void readData(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo, 
				    const StationData& sd, const unsigned int& stationindex);

		A3DWebServiceSoap12BindingProxy gsn;
		Config cfg;
		std::vector<std::string> vecStationName;
		std::string endpoint, hostname, port, userid, passwd; ///< Variables for proxy configuration
		int proxyport;                              ///< Variable for proxy configuration
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		double in_tz, out_tz; //default time zones
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
};

} //end namespace mio

#endif
