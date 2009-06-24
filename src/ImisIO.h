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
#include "MeteoBuffer.h"
#include "Meteo1DResampler.h"
#include "IOExceptions.h"


using namespace std;
using namespace IOUtils;

/**
 * @class ImisIO 
 * @brief The class with-in the data from the database are treated. The MeteoData and the StationData will be set in.
 * This class also herited to IOInterface class which is abstract.
 * @author Moustapha Mbengue
 * @date 2009-05-12
 */

class ImisIO : public IOInterface {
	public:
		ImisIO(const string& configfile);
		ImisIO(void (*delObj)(void*), const string& filename);
		//ImisIO(const ImisIO&);
		virtual ~ImisIO() throw();
			
		/**
		 * @brief Once the data's retrieved form the database, this method is used to create a MeteoBuffer which contain 
		 * the meteo data and the station data of each single station in the configfile.
		 * @param meteo_in <vector<string>> : meteo data from the database.
		 * @param station_in <string> : station data from the database.
		 * @param mb MeteoBuffer : variable in which stationdata and meteodata are filled.
		 */
		void createData(vector< vector<string> >& meteo_in, vector<string>& station_in, MeteoBuffer& mb);
		
		virtual void get2DGridSize(int& nx, int& ny);
		
		virtual void read2DGrid(Grid2DObject& grid_out, const string& parameter="");
		
		virtual void readDEM(Grid2DObject& dem_out);
		
		virtual void readLanduse(Grid2DObject& landuse_out);
		
		virtual void readAssimilationData(const Date_IO& date_in, Grid2DObject& da_out);
		
		virtual void readSpecialPoints(CSpecialPTSArray& pts);
		
		virtual void write2DGrid(const Grid2DObject& grid_in, const string& options="");
	  	
		/**
		 * @brief refer to void readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation)
		 */
		virtual void readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo);
		
		/**
		 * @brief this method is used to get back meteo and station data for a given day
		 * into the respective vectors MeteoData and StationData
		 * @param date_in Date_IO: recording date
		 * @param vecMeteo <MeteoData>: vector in which meteo data will be filled
		 * @param vecStation <StationData>: vector in which station data will be filled
		 */
	  	virtual void readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation);
		
		
		/**
		 * @brief ImisIO cleaner function
		 */
		void cleanup() throw();
		
		/**
		 * @brief Get back station's name from the configfile
		 */		
		void getStationName();
		
		/**
		 * @brief Returns the ConfigFile
		 */
		ConfigReader getCfg();
		
		/**
		 * @brief Returns vecStationName, a string vector in which StationName are saved
		 */
		vector<string> getVecStationName();
		
		/**
		 * @brief Returns mbImis
		 */
		vector<MeteoBuffer> getMbImis();
		
		/**
		 * @brief Method which allows putting meteodata and stationdata into mbImis (vector<MeteoBuffer>)
		 * @param date_in Date_IO: recording date 
		 */
		void setMbImis(Date_IO date_in);
		
		/**
		 * @brief Selects the data corresponding to date_in. But whether we're looking for a date which is not in the database
		 * but that is between two recording date, this method will interpolate data for this given date.
		 * @param vecMeteo vector<MeteoData>: vector of meteodata corresponding to date_in
		 * @param vecStation vector<StationData>: vector of stationdata corresponding to date_in
		 * @param date_in const Date_IO: recording date
		 */
		void resampleMbImis(vector<MeteoData>& vecMeteo, vector<StationData>& vecStation, const Date_IO& date_in);
		
		/**
		 * @brief Defines the maximum size of MeteoBuffer (mbImis)
		 */
		void createBuffer();
		
 	private:
		
 		/**
 		 * @brief This is a private function. It gets back data from station2 which is a table of the database and fill them in a string vector
 		 * @param stat_abk : a string key of station2 
 		 * @param stao_nr : an integer key of station2
 		 * @param data2S : string vector in which data will be filled
 		 */
		void getStation2Data(const string stat_abk, unsigned int stao_nr, vector<string>& data2S);
		
		/**
		 * @brief This is a private function. It gets back data from ams.v_imis which is a table of the database
		 * and fill them in a vector of vector of string. It seems that each record is a string vector
		 * @param stat_abk : a string key of ams.v_imis
		 * @param stao_nr : an integer key of ams.v_imis
		 * @param date_in : a vector of five(5) integer corresponding to the recording date
		 * @param datatImis : a vector of vector of string in which data will be filled
		 */
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

