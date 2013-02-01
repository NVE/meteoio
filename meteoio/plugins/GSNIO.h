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

#include <math.h>
#define _isnan isnan

//This is big hack to allow older versions of gcc to turn warnings off in soap files
//hopefully we would soon get better soap files and would scrap the whole mess...
#if defined(__GNUC__) && defined(__DEBUG)
	//#pragma GCC diagnostic push_options
	#pragma GCC diagnostic ignored "-Wall"
	#pragma GCC diagnostic ignored "-Weffc++"
	#pragma GCC diagnostic ignored "-Wlong-long"
	#pragma GCC diagnostic ignored "-Wshadow"
	#include "gsn/soapGSNWebServiceSoap12BindingProxy.h"
	//#include "gsn/GSNWebServiceSoap12Binding.nsmap"
	//#pragma GCC diagnostic pop_options
	#pragma GCC diagnostic warning "-Wall"
	#pragma GCC diagnostic warning "-Weffc++"
	#pragma GCC diagnostic warning "-Wlong-long"
	#pragma GCC diagnostic warning "-Wshadow"
#else
	#include "gsn/soapGSNWebServiceSoap12BindingProxy.h"
	//#include "gsn/GSNWebServiceSoap12Binding.nsmap"
#endif

#ifdef WIN32 //because we collected c**p from windows.h that was included by gsoap
	#undef max
	#undef min
#endif

#include <meteoio/Config.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOUtils.h>
#include <meteoio/Coords.h>
#include <meteoio/IOExceptions.h>

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>

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
		GSNIO(const std::string& configfile);
		GSNIO(const GSNIO&);
		GSNIO(const Config&);
		~GSNIO() throw();

		virtual void read2DGrid(Grid2DObject& dem_out, const std::string& parameter="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo,
		                           const size_t& stationindex=IOUtils::npos);

		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo,
		                            const std::string& name="");


		virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
		virtual void readSpecialPoints(std::vector<Coords>& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& name);
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);

	private:
		void listSensors(std::vector<std::string>& vec_names);
		void convertUnits(MeteoData& meteo);
		void initGSNConnection();
		void readStationNames();
		void readMetaData();
		void readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex);
		void map_parameters(const std::vector<ns2__GSNWebService_USCOREDataField*>& field, MeteoData& md,
		                    std::vector<size_t>& index);
		double olwr_to_tss(const double& olwr);
		void parse_streamElement(const std::vector<size_t>& index, const bool& olwr_present,
				  std::vector<MeteoData>& vecMeteo, MeteoData& tmpmeteo, ns2__GSNWebService_USCOREStreamElement* streamElement);

		GSNWebServiceSoap12BindingProxy gsn;
		const Config& cfg;
		std::vector<std::string> vecStationName;
		std::vector<StationData> vecMeta;
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		std::string endpoint, hostname, port, userid, passwd; ///< Variables for proxy configuration
		int proxyport;                                        ///< Variable for proxy configuration
		double default_timezone;
};

} //end namespace mio

#endif
