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
#include "IOInterface.h"
#include "ConfigReader.h"
#include "IOUtils.h"
#include "IOExceptions.h"
#include "DynamicLibrary.h"

#include <string>
#include <sstream>
#include <iostream>
#include <vector>

/**
 * @class GSNIO
 * @brief This class enables the access to the GSN web service
 *
 * @author Thomas Egger
 * @date   2009-09-25
 */
class GSNIO : public IOInterface {
	public:
		GSNIO(void (*delObj)(void*), const string& filename);

		GSNIO(const std::string& configfile);
		GSNIO(const GSNIO&);
		GSNIO(const ConfigReader&);
		~GSNIO() throw();

		virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  const unsigned int& stationindex=IOUtils::npos);

		virtual void read2DGrid(Grid2DObject& dem_out, const string& parameter="");
		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);
		virtual void readAssimilationData(const Date_IO&, Grid2DObject& da_out);
		virtual void readSpecialPoints(CSpecialPTSArray& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const string& name);

	private:
		void parseString(const std::string& _string, std::vector<std::string>& vecString, MeteoData& md);
		void convertStringToDouble(double& d, const std::string& _string, const std::string& _parname);
		void convertUnits(MeteoData& meteo);
		void initGSNConnection();
		void readStationNames();
		void readStationMetaData(StationData& sd, const unsigned int& stationindex);
		void readData(const Date_IO& dateStart, const Date_IO& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo, 
				    std::vector< std::vector<StationData> >& vecStation, const StationData& sd, const unsigned int& stationindex);

		A3DWebServiceSoap12BindingProxy gsn;
		ConfigReader cfg;
		vector<string> vecStationName;
		std::string hostname, port, userid, passwd; ///< Variables for proxy configuration
		int proxyport;                              ///< Variable for proxy configuration
};

#endif
