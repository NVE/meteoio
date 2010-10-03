/***********************************************************************************/
/*  Copyright 2009, 2010 WSL Institute for Snow and Avalanche Research   SLF-DAVOS */
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
#ifndef __IMISIO_H__
#define __IMISIO_H__

#include <meteoio/Config.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOUtils.h>
#include <meteoio/Coords.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/DynamicLibrary.h>
#include <meteoio/BufferedIOHandler.h>

#include <string>
#include <sstream>
#include <set>
#include <map>
#include <iostream>
#include <vector>
#include <ctime>
#include <occi.h>

namespace mio {

	class AnetzData{
		public:
			AnetzData(): nrOfAnetzStations(0), nrOfCoefficients(0)
			{
				anetzstations.push_back("");
				anetzstations.push_back("");
				anetzstations.push_back("");
				coeffs.push_back(IOUtils::nodata);
				coeffs.push_back(IOUtils::nodata);
				coeffs.push_back(IOUtils::nodata);
			}

			AnetzData(const unsigned int& nr_anetz,
					const std::string& i_anetz1, const std::string& i_anetz2, const std::string& i_anetz3, 
					const unsigned int& nr_coeffs, 
					const double& coeff1, const double& coeff2, const double& coeff3) 
				: nrOfAnetzStations(nr_anetz), nrOfCoefficients(nr_coeffs)
			{
				anetzstations.push_back(i_anetz1);
				anetzstations.push_back(i_anetz2);
				anetzstations.push_back(i_anetz3);
				coeffs.push_back(coeff1);
				coeffs.push_back(coeff2);
				coeffs.push_back(coeff3);
			}

			std::string anetzstation1, anetzstation2, anetzstation3;
			unsigned int nrOfAnetzStations, nrOfCoefficients;
			std::vector<double> coeffs;
			std::vector<std::string> anetzstations;
	};

/**
 * @class ImisIO
 * @brief The class with-in the data from the database are treated. The MeteoData and the StationData will be set in.
 * This class also herited to IOInterface class which is abstract.
 * @author Moustapha Mbengue and Thomas Egger
 * @date 2009-05-12
 */
class ImisIO : public IOInterface {
	public:
		ImisIO(void (*delObj)(void*), const Config& i_cfg);

		ImisIO(const std::string& configfile);
		ImisIO(const ImisIO&);
		ImisIO(const Config&);
		~ImisIO() throw();

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

	private:
		void cleanup() throw();
		void getDBParameters();
		void getStationData(const std::string& stat_abk, const std::string& stao_nr, std::vector<std::string>& data2S);
		void getImisData(const std::string& stat_abk, const std::string& stao_nr,
		                 const std::vector<int>& datestart, const std::vector<int>& dateend,
		                 std::vector< std::vector<std::string> >& dataImis);
		void parseDataSet(const std::vector<std::string>& meteo_in, MeteoData& md);
		void readData(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo,
		              std::vector< std::vector<StationData> >& vecStation, const unsigned int& stationindex);
		void readStationNames(std::vector<std::string>& vecStationName);
		void parseStationName(const std::string& stationName, std::string& stName, std::string& stNumber);
		void readStationMetaData();
		void convertUnits(MeteoData& meteo);
		void initializeAnetzBuffer(const unsigned int& indexStart, const unsigned int& indexEnd,
							  std::map<std::string, unsigned int>& mapAnetzNames, Config& anetzcfg);
		void assimilateAnetzData(const unsigned int& indexStart, const unsigned int& indexEnd,
							std::vector< std::vector<MeteoData> >& vecMeteo,
							std::vector< std::vector<StationData> >& vecStation,		   
							std::map<std::string, unsigned int>& mapAnetzNames, Config& anetzcfg);
		double getHNW(const std::vector<MeteoData>& vecAnetz, const AnetzData& ad, const unsigned int& index, 
				    const std::map<std::string, unsigned int>& mapAnetzNames);

		static const double in_tz; //timezone
		Config cfg;
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		std::vector<StationData> vecMyStation;
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		static const std::string sqlQueryMeteoData;
		static const std::string sqlQueryStationData;
		std::string oracleUserName_in;
		std::string oraclePassword_in;
		std::string oracleDBName_in;
		bool useAnetz;
		/*std::string oracleUserName_out;
		std::string oraclePassword_out;
		std::string oracleDBName_out;*/

		static std::map<std::string, AnetzData> mapAnetz;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map meteoparamname 
};

} //end namespace mio

#endif

