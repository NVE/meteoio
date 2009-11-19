#ifndef __IMISIO_H__
#define __IMISIO_H__


#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <occi.h>

#include "IOInterface.h"
#include "ConfigReader.h"
#include "IOUtils.h"
#include "DynamicLibrary.h"
#include "MeteoData.h"
#include "StationData.h"
#include "IOExceptions.h"

#define IMIS_BUFF_SIZE 5000

using namespace std;
using namespace IOUtils;

class ImisIO : public IOInterface {
	public:
		ImisIO(const string& configfile);
		ImisIO(void (*delObj)(void*), const string& filename);
		//ImisIO(const ImisIO&);
		virtual ~ImisIO() throw();
				
		virtual void read2DGrid(Grid2DObject& grid_out, const string& parameter="");
		
		virtual void readDEM(DEMObject& dem_out);
		
		virtual void readLanduse(Grid2DObject& landuse_out);
		
		virtual void readAssimilationData(const Date_IO& date_in, Grid2DObject& da_out);
		
		virtual void readSpecialPoints(CSpecialPTSArray& pts);
		
		virtual void write2DGrid(const Grid2DObject& grid_in, const string& name="");
		
		virtual void readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo);
		
	  	virtual void readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation);
		
		/*virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  const unsigned int& stationindex=IOUtils::npos);*/
		
	private:
		void createData(vector< vector<string> >& meteo_in, vector<string>& station_in, MeteoBuffer& mb);

		void cleanup() throw();
			
		void getStationName();
		
		ConfigReader getCfg();
		
		vector<string> getVecStationName();
		
		/**
		 * @brief Returns mbImis
		 */
		vector<MeteoBuffer> getMbImis();
		
		/**
		 * @brief Method which allows putting meteodata and stationdata into mbImis (vector<MeteoBuffer>)
		 * @param date_in Date_IO: recording date
		 * @param stationName const string: string keys for database's queries
		 * @param buffer MeteoBuffer: container in which data're filled
		 */
		void setMbImis(Date_IO date_in, const string& stationName, MeteoBuffer& buffer);
		
		/**
		 * @brief Selects the data corresponding to date_in. But whether we're looking for a date which is not in the database
		 * but that is between two recording date, this method will interpolate data for this given date.
		 * @param meteo MeteoData:  a meteodata corresponding to date_in
		 * @param station StationData:  a stationdata corresponding to date_in
		 * @param date_in const Date_IO: recording date
		 * @param mb MeteoBuffer: data's container
		 */
		void resampleMbImis(MeteoData& meteo, StationData& station, const Date_IO& date_in, MeteoBuffer& mb);
		
		/**
		 * @brief Defines the maximum size of MeteoBuffer (mbImis)
		 */
		void createBuffer();
		
		void getStation2Data(const string stat_abk, unsigned int stao_nr, vector<string>& data2S);
		
		void getImisData(const string &stat_abk, const unsigned int &stao_nr, vector<int> date_in, vector< vector<string> >& dataImis);
		
		/**
		 * @brief Variable in which the configFile will be saved
		 */
		ConfigReader cfg;
				
		/**
		 * @brief A vector for StationName's saving
		 */
		vector<string> vecStationName;
		
		/**
		 * @brief A vector in which MeteoData and StationData are saved for each station in "io.ini".
		 * In MeteoBuffer class, there're two (2) deque :
		 * - meteobuffer : in which MeteoData are saved
		 * - stationbuffer : in which StationData are saved
		 */		
		vector<MeteoBuffer> mbImis;
};

#endif

