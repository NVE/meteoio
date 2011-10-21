/***********************************************************************************/
/*  Copyright 2010 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef __SMETIO_H__
#define __SMETIO_H__

#include <meteoio/IOInterface.h>
#include <meteoio/Config.h>
#include <meteoio/libsmet.h>

#include <string>

namespace mio {

/**
 * @class SMETIO
 * @brief This (empty) class is to be used as a template for developing new plugins
 *
 * @ingroup plugins
 * @author Mathias Bavay
 * @date   2010-06-14
 */
class SMETIO : public IOInterface {
	public:
		SMETIO(void (*delObj)(void*), const Config& i_cfg);

		SMETIO(const std::string& configfile);
		SMETIO(const SMETIO&);
		SMETIO(const Config& cfgreader);
		~SMETIO() throw();

		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");

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
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& filename);

	private:
		void read_meta_data(const smet::SMETReader& myreader, StationData& meta);
		void identify_fields(const std::vector<std::string>& fields, std::vector<size_t>& indexes,
		                     bool& julian_present, MeteoData& md);
		void copy_data(const smet::SMETReader& myreader, const std::vector<std::string>& timestamps,
 		               const std::vector<double>& mydata, std::vector<MeteoData>& vecMeteo);

		void parseInputOutputSection();
		bool checkConsistency(const std::vector<MeteoData>& vecMeteo, StationData& sd);
		size_t getNrOfParameters(const std::string& stationname, const std::vector<MeteoData>& vecMeteo);
		void checkForUsedParameters(const std::vector<MeteoData>& vecMeteo, const size_t& nr_parameters, double& tz,
		                            std::vector<bool>& vecParamInUse, std::vector<std::string>& vecColumnName);
		void getFormatting(const size_t& param, int& prec, int& width);
		double olwr_to_tss(const double& olwr);
		void generateHeaderInfo(const StationData& sd, const bool& outputIsAscii, const bool& isConsistent,
		                        const double& timezone, const size_t& nr_of_parameters,
		                        const std::vector<bool>& vecParamInUse,
		                        const std::vector<std::string>& vecColumnName,
		                        smet::SMETWriter& mywriter);

		size_t nr_stations; //number of stations to read from
		std::vector<std::string> vecFiles;  //read from the Config [Input] section
		std::string outpath;                //read from the Config [Output] section
		bool outputIsAscii, outputIsGzipped;//read from the Config [Output] section
		double in_dflt_TZ, out_dflt_TZ;     //default time zones
		double plugin_nodata;

		Config cfg;
		std::string coordin, coordinparam, coordout, coordoutparam; //default projection parameters
		std::vector<smet::SMETReader> vec_smet_reader;
};

} //namespace
#endif
