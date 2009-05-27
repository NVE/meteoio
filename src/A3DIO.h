#ifndef __A3DIO_H__
#define __A3DIO_H__

#include "IOInterface.h"

#include "FilterFacade.h"
#include "ConfigReader.h"
#include "IOExceptions.h"
#include "IOUtils.h"
#include "MeteoBuffer.h"
#include "Meteo1DResampler.h"

#include <string>
#include <vector>
#include <map>

using namespace IOUtils;

class A3DIO : public IOInterface {
 public:
	//virtual A3DIO* clone() const;

	//A3DIO(void (*delObj)(void*), const std::string& filename);

	A3DIO(const std::string& configfile);
	A3DIO(const A3DIO&);
	A3DIO(const ConfigReader&);
	~A3DIO() throw();

	virtual void get2DGridSize(int& nx, int& ny);
	virtual void read2DGrid(Grid2DObject& dem_out, const string& parameter="");

	virtual void readDEM(Grid2DObject& dem_out);
	virtual void readLanduse(Grid2DObject& landuse_out);

	virtual void readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo);
	virtual void readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation);

	virtual void readAssimilationData(const Date_IO&, Grid2DObject& da_out);
	virtual void readSpecialPoints(CSpecialPTSArray& pts);

	virtual void write2DGrid(const Grid2DObject& grid_in, const string& filename);

 private:
	void read1DMeteo(const Date_IO&, MeteoData&); ///< No buffering
	void read1DMeteo(const Date_IO&, MeteoData&, StationData&); ///< No buffering
	void read2DMeteo(const Date_IO&, vector<MeteoData>&); ///< No buffering
	void read2DMeteo(const Date_IO&, vector<MeteoData>&, vector<StationData>&); ///< No buffering

	void constructMeteo2DFilenames(const Date_IO& date_in, vector<string>& filenames);
	void readMeteoDataLine(std::string& line, const Date_IO& date_in, MeteoData& tmpdata, std::string filename);
	void convertUnits(MeteoData& meteo);
	void cleanup() throw();
	void read2DMeteoData(const string&, const string&, const Date_IO&, 
					 map<string,unsigned int>& hashStations, vector<MeteoData>&, unsigned int& bufferindex);
	void read2DMeteoHeader(const string& filename, map<string, unsigned int>& hashStations, vector<StationData>&);
	//unsigned int getNrOfStations(const string& filename);
	unsigned int getNrOfStations(vector<string>& filenames, map<string, unsigned int>& hashStations);

	ConfigReader cfg;
	ifstream fin; //Input file streams
	ofstream fout;//Output file streams
	vector<MeteoBuffer> unfilteredMeteoBuffer;
	vector<MeteoBuffer> filteredMeteoBuffer;
};

#endif
