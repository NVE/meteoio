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
 * This class also herited to IOhandler class which is abstract.
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
		 * @brief Once the data're got form the database, this method is used to fill the meteo data and the station data
		 * into the respective vectors MeteoData and StationData
		 * @param vecMeteo <MeteoData>: vector in which meteo data will be filled
		 * @param vecStation <StationData>: vector in which station data will be filled
		 * @param meteo_in <vector<string>>: meteo data from the database
		 * @param station_in <string>: station data from the database
		 */
		void createData(vector<MeteoData>& vecMeteo, vector<StationData>& vecStation,
				vector< vector<string> >& meteo_in, vector<string>& station_in, MeteoBuffer& mb);
		
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
		 * @brief Get back station's name from "io.ini"
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

		void test(vector<int> date);	
		
		/**
		 * @brief It allows to convert a date which is in string format to SLFIO date/time format
		 * @param instr : string to convert
		 * @param date_out : date in Date_IO format (SFLIO date/time)
		 */
		static void stringToDate(const string& instr, Date_IO& date_out);
		
		/**
		 * @brief Data from the db are all string, this function convert the double in string format to real double
		 * @param str : string to convert into double
		 * @return : double which corresponds to the double in string format
		 */
		static double strToDouble(const string &str);
		
		/**
		 * @brief function which display meteo and station data
		 */
		void displayData(); 
		
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
		//ifstream fin; //Input file streams		
				
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

