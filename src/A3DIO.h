#ifndef __A3DIO_H__
#define __A3DIO_H__

#include "IOInterface.h"

#include "ConfigReader.h"
#include "IOExceptions.h"
#include "IOUtils.h"
#include "MeteoBuffer.h"

#include <string>
#include <vector>
#include <map>

class A3DIO : public IOInterface {
	public:
		//virtual A3DIO* clone() const;

		//A3DIO(void (*delObj)(void*), const std::string& filename);

		A3DIO(const std::string& configfile);
		A3DIO(const A3DIO&);
		A3DIO(const ConfigReader&);
		~A3DIO() throw();

		virtual void read2DGrid(Grid2DObject& dem_out, const string& parameter="");

		virtual void readDEM(Grid2DObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);
		virtual void readAssimilationData(const Date_IO&, Grid2DObject& da_out);

		virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo);

		virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  const unsigned int& stationindex=IOUtils::npos);

		virtual void readSpecialPoints(CSpecialPTSArray& pts);

		virtual void write2DGrid(const Grid2DObject& grid_in, const string& name);

	private:
		void read1DMeteo(const Date_IO& dateStart, const Date_IO& dateEnd, 
					  std::vector< std::vector<MeteoData> >&, std::vector< std::vector<StationData> >&); ///< No buffering
		void read2DMeteo(std::vector< std::vector<MeteoData> >&, std::vector< std::vector<StationData> >&); ///< No buffering

		void constructMeteo2DFilenames(const Date_IO& _date, std::vector<string>& _filenames);
		bool readMeteoDataLine(std::string& line, const Date_IO& date_in, MeteoData& tmpdata, std::string filename);
		void convertUnits(MeteoData& meteo);
		void cleanup() throw();
		void read2DMeteoData(const std::string&, const std::string&, std::map<std::string,unsigned int>& hashStations, 
						 std::vector< std::vector<MeteoData> >&, unsigned int& bufferindex);
		void read2DMeteoHeader(const string& filename, std::map<std::string, unsigned int>& hashStations, 
						   std::vector<StationData>&);
		unsigned int getNrOfStations(std::vector<std::string>& filenames, 
							    std::map<std::string, unsigned int>& hashStations);

		ConfigReader cfg;
		ifstream fin; //Input file streams
		ofstream fout;//Output file streams
};

#endif
