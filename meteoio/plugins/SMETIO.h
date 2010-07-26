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

#include <string>

namespace mio {

/**
 * @class SMETIO
 * @brief This (empty) class is to be used as a template for developing new plugins
 *
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
		                           std::vector< std::vector<StationData> >& vecStation,
		                           const unsigned int& stationindex=IOUtils::npos);

		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, 
		                            const std::vector< std::vector<StationData> >& vecStation,
		                            const std::string& name="");

		virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
		virtual void readSpecialPoints(std::vector<Coords>& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& filename);

	private:
		static const std::string smet_version;
		static std::map<std::string, MeteoData::Parameters> mapParameterByName; ///<Associate name and meteo parameter
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map meteoparamname 
		static double& getParameter(const std::string& columnName, MeteoData& md);
		static void checkColumnNames(const std::vector<std::string>& vecColumns, const bool& locationInHeader);

		void cleanup() throw();
		void parseInputOutputSection();
		bool checkConsistency(const std::vector<StationData>& vecStation, StationData& sd);
		void checkForUsedParameters(const std::vector<MeteoData>& vecMeteo, double& timezone, 
                                      std::vector<bool>& vecParamInUse);
		void writeHeaderSection(const bool& writeLocationInHeader, const StationData& sd, 
                                  const double& timezone, const std::vector<bool>& vecParamInUse);
		void writeDataAscii(const bool& writeLocationInHeader, const std::vector<MeteoData>& vecMeteo,
                              const std::vector<StationData>& vecStation, const std::vector<bool>& vecParamInUse);
		void writeDataBinary(const bool& writeLocationInHeader, const std::vector<MeteoData>& vecMeteo,
                               const std::vector<StationData>& vecStation, const std::vector<bool>& vecParamInUse);

		void readHeader(const char& eoln, const std::string& filename, bool& locationInHeader,
                          double& timezone, StationData& sd, std::vector<std::string>& vecDataSequence);
		void readDataAscii(const char& eoln, const std::string& filename, const double& timezone,
					    const StationData& sd, const std::vector<std::string>& vecDataSequence,
					    const Date& dateStart, const Date& dateEnd,
					    std::vector<MeteoData>& vecMeteo, std::vector<StationData>& vecStation);
		void readDataBinary(const char& eoln, const std::string& filename, const double& timezone,
					     const StationData& sd, const std::vector<std::string>& vecDataSequence,
					     const Date& dateStart, const Date& dateEnd,
					     std::vector<MeteoData>& vecMeteo, std::vector<StationData>& vecStation);

		void checkSignature(const std::vector<std::string>& vecSignature, const std::string& filename, bool& isAscii);
		void setFormatting(const MeteoData::Parameters& paramindex);

		std::vector<std::string> vecFiles;  //read from the Config [Input] section
		std::string outpath;                //read from the Config [Output] section
		bool outputIsAscii, outputIsGzipped;//read from the Config [Output] section

		Config cfg;
		std::ifstream fin; //Input file streams
		std::ofstream fout; //Output file streams

};

} //namespace
#endif
