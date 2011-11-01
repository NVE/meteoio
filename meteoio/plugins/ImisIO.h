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
		AnetzData(): nrOfAnetzStations(0), nrOfCoefficients(0) {
			anetzstations.push_back("");
			anetzstations.push_back("");
			anetzstations.push_back("");
			coeffs.push_back(IOUtils::nodata);
			coeffs.push_back(IOUtils::nodata);
			coeffs.push_back(IOUtils::nodata);
		}

		AnetzData(const size_t& nr_anetz,
		          const std::string& i_anetz1, const std::string& i_anetz2, const std::string& i_anetz3,
		          const size_t& nr_coeffs,
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
		size_t nrOfAnetzStations, nrOfCoefficients;
		std::vector<double> coeffs;
		std::vector<std::string> anetzstations;
};

/**
 * @class ImisIO
 * @brief The class with-in the data from the database are treated. The MeteoData and the StationData will be set in.
 * This class also herited to IOInterface class which is abstract.
 *
 * @ingroup plugins
 * @author Thomas Egger, first version by Moustapha Mbengue
 * @date 2010-11-02
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
		                           const size_t& stationindex=IOUtils::npos);

		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo,
		                            const std::string& name="");

		virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
		virtual void readSpecialPoints(std::vector<Coords>& pts);

		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& name);

	private:
		void cleanup() throw();

		void openDBConnection(oracle::occi::Environment*& env, oracle::occi::Connection*& conn);
		void closeDBConnection(oracle::occi::Environment*& env, oracle::occi::Connection*& conn);
		void getDBParameters();

		size_t getStationIDs(const std::string& stat_code,
		                     const std::string& sqlQuery, std::vector<std::string>& vecStationMetaData,
		                     oracle::occi::Connection*& conn);
		size_t getStationMetaData(const std::string& stat_abk, const std::string& stao_nr,
		                          const std::string& sqlQuery, std::vector<std::string>& vecStationMetaData,
		                          oracle::occi::Connection*& conn);
		size_t getSensorDepths(const std::string& stat_abk, const std::string& stao_nr,
		                       const std::string& sqlQuery, std::vector<std::string>& vecHTS1,
		                       oracle::occi::Connection*& conn);
		bool getStationData(const std::string& stat_abk, const std::string& stao_nr,
		                    const std::vector<int>& datestart, const std::vector<int>& dateend,
		                    const std::vector<std::string>& i_vecHTS1,
		                    std::vector< std::vector<std::string> >& vecMeteoData,
		                    oracle::occi::Environment*& env, oracle::occi::Connection*& conn);

		void parseDataSet(const std::vector<std::string>& meteo_in, MeteoData& md, bool& _fullStation);
		void readData(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo,
		              const size_t& stationindex, const std::vector<StationData>& vecStationID,
		              oracle::occi::Environment*& env, oracle::occi::Connection*& conn);
		void readStationIDs(std::vector<std::string>& vecStationID);
		void parseStationID(const std::string& stationID, std::string& stnAbbrev, std::string& stnNumber);

		void readStationMetaData(oracle::occi::Connection*& conn);
		void convertSnowTemperature(MeteoData& meteo, const std::string& parameter);
		void convertSensorDepth(MeteoData& meteo, const std::string& parameter);
		void convertUnits(MeteoData& meteo);

		//helper functions for the Anetz coefficient mangling:
		void findAnetzStations(const size_t& indexStart, const size_t& indexEnd,
		                       std::map<std::string, size_t>& mapAnetzNames, std::vector<StationData>& vecAnetzStation);
		void getAnetzHNW(const AnetzData& ad, const std::map<std::string, size_t>& mapAnetzNames,
		                 const std::vector< std::vector<double> >& vec_of_psums, std::vector<double>& psum);
		void assimilateAnetzData(const Date& dateStart, const AnetzData& ad,
		                         const std::vector< std::vector<double> > vec_of_psums,
		                         const std::map<std::string, size_t>& mapAnetzNames, const size_t& stationindex,
		                         std::vector< std::vector<MeteoData> >& vecMeteo);
		void calculatePsum(const Date& dateStart, const Date& dateEnd,
		                   const std::vector< std::vector<MeteoData> >& vecMeteoAnetz,
		                   std::vector< std::vector<double> >& vec_of_psums);

		static const double in_tz; //timezone
		Config cfg;
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		std::vector<StationData> vecStationMetaData;
		std::map<std::string, std::string> mapDriftStation;
		std::vector<double> vecHTS1; // vector of sensor depths
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		static const std::string sqlQueryStationIDs;
		static const std::string sqlQueryStationMetaData;
		static const std::string sqlQuerySensorDepths;
		static const std::string sqlQueryMeteoDataDrift; // combined snow_drift query from two stations (ams.v_ams_raw)
		static const std::string sqlQueryMeteoData; // single station query (ams.v_ams_raw)
		std::string oracleUserName_in;
		std::string oraclePassword_in;
		std::string oracleDBName_in;
		bool useAnetz;

		static std::map<std::string, AnetzData> mapAnetz;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map meteoparamname
};

} //end namespace mio

#endif

