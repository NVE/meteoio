#ifndef __BUFFEREDIOHANDLER_H__
#define __BUFFEREDIOHANDLER_H__

#include "IOInterface.h"
#include "ConfigReader.h"
#include "Meteo1DResampler.h"

#include <map>
#include <vector>
#include <string>

/**
 * @class BufferedIOHandler
 * @brief This class serves as a wrapper around all children of IOInterface. It internally handles 
 *        the buffering of the data and introduces a few convenient functions to access meteo data 
 *
 * @author Thomas Egger
 * @date   2009-07-25
 */
#ifdef _POPC_
class BufferedIOHandler {
#else
class BufferedIOHandler : public IOInterface {
#endif
	public:
	
		/**
		 * @brief The constructor accepts an already initialized child of IOInterface (e.g. A3DIO, BoschungIO, ImisIO)
		 *        and a ConfigReader object
		 *    
		 * Example Usage:
		 * @code
		 * IOInterface *io1;
		 * ConfigReader cfg("io.ini");
		 * io1 = new A3DIO(cfg);
		 * BufferedIOHandler bio(*io1, cfg);
		 * @endcode
		 */
		BufferedIOHandler(IOInterface& _iohandler, const ConfigReader& _cfg);
	#ifdef _POPC_
		~BufferedIOHandler();
	#else
		~BufferedIOHandler() throw();
	#endif

		/**
		 * @brief The function returns the next MeteoData object for each station with a 
		 *        date >= to the parameter _date. vecMeteo and vecStation will be empty if there 
		 *        is no MeteoData to be found.
		 *        NOTE: only the real measured data is looked at: no resampled values are taken into account
		 *
		 * @param _date start date of the data search for each station
		 * @param vecMeteo   A vector of MeteoData objects to be filled with data
		 * @param vecStation A vector of StationData objects to be filled with meta data
		 */
		void getNextMeteoData(const Date_IO& _date, std::vector<MeteoData>& vecMeteo, std::vector<StationData>& vecStation);
		
		/**
		 * @brief See BufferedIOHandler::readMeteoData(const Date& date_in, 
		 *                                             vector<MeteoData>& vecMeteo, 
		 *                                             vector<StationData>& vecStation).
		 */
		void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo);

		/**
		 * @brief Fill vector<MeteoData> and vector<StationData> objects with multiple datasets 
		 * corresponding to the time indicated by the Date object.
		 * Matching rule: Find first data set for every station which has an event time (measurement time) 
		 * that is greater (newer) or equal to the time represented by the Date object parameter. The 
		 * vector<StationData> object holds multiple StationData objects representing meta information 
		 * about the meteo stations that recorded the meteo data.
		 *
		 * NOTE: 
		 * - vecMeteo and vecStation will contain nodata objects if an exact time match is impossible
		 *   and resampling is turned off. If resampling is turned on a resampled value is returned 
		 *   if resampling is possible (enough measurements), otherwise nodata objects will be returned
		 * - is there is absolutely no data to be found, and hence not even station data, vecMeteo and vecStation
		 *   will be filled with only one nodata obejct of MeteoData and StationData respectively
		 *    
		 * Example Usage:
		 * @code
		 * vector<MeteoData> vecMeteo;      //empty vector
		 * vector<StationData> vecStation;  //empty vector
		 * BufferedIOHandler bio(A3DIO("io.ini"), ConfigReader("io.ini"));
		 * bio.readMeteoData(Date(2008,06,21,11,00), vecMeteo, vecStation); //21.6.2008 11:00
		 * @endcode
		 * @param _date       A Date object representing the date/time for the sought MeteoData objects
		 * @param vecMeteo    A vector of MeteoData objects to be filled with data
		 * @param vecStation  A vector of StationData objects to be filled with data
		 */
		void readMeteoData(const Date_IO& _date, std::vector<MeteoData>& vecMeteo, std::vector<StationData>& vecStation);

		/**
		 * @brief Clear all buffers in BufferedIOHandler and hence force rebuffering
		 */
		void clearBuffer();

		/**
		 * @brief Turn on/off resampling. It is turned on by default
		 * @param _enable true turns on resampling, false turns it off
		 */
		void enableResampling(const bool& _enable);

		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		virtual void readDEM(Grid2DObject& dem_out);
		virtual void readAssimilationData(const Date_IO& date_in, Grid2DObject& da_out);
		virtual void readLanduse(Grid2DObject& landuse_out);
		virtual void readSpecialPoints(CSpecialPTSArray& pts);
		virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  const unsigned int& stationindex=IOUtils::npos);

		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& options="");


		static const unsigned int npos = (unsigned int)-1;             ///<npos is the out-of-range value

	private:
		unsigned int seek(const Date_IO& date_in, std::vector<MeteoData>& vecM);
		bool bufferData(const Date_IO& _date, const unsigned int& stationindex);
		void bufferAllData(const Date_IO& _date);
		
		IOInterface& iohandler;
		ConfigReader cfg;
		std::vector< std::vector<MeteoData> > meteoBuffer;
		std::vector< std::vector<StationData> > stationBuffer;
		std::map<std::string, Grid2DObject> mapBufferedGrids;
		bool resample;
};
#endif
