#ifndef __ASCIIFILEIO_H__
#define __ASCIIFILEIO_H__

#include "IOHandler.h"

#include "FilterFacade.h"
#include "ConfigReader.h"
#include "slfexceptions.h"
#include "slfutils.h"
#include "MeteoBuffer.h"
#include "Meteo1DResampler.h"

#include <string>
#include <vector>
#include <map>

class ASCIIFileIO : public IOHandler {
 public:
	//virtual ASCIIFileIO* clone() const;

	//ASCIIFileIO(void (*delObj)(void*), const std::string& filename);

	ASCIIFileIO(const std::string& configfile);
	ASCIIFileIO(const ASCIIFileIO&);
	ASCIIFileIO(const ConfigReader&);
	~ASCIIFileIO() throw();

	virtual void get2DGridSize(int& nx, int& ny);
	virtual void read2DGrid(Grid2DObject& dem_out, const string& parameter="");

	virtual void readDEM(Grid2DObject& dem_out);
	virtual void readLanduse(Grid2DObject& landuse_out);

	virtual void readMeteoData(const Date& date_in, vector<MeteoData>& vecMeteo);
	virtual void readMeteoData(const Date& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation);

	virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
	virtual void readSpecialPoints(CSpecialPTSArray& pts);

	virtual void write2DGrid(const Grid2DObject& grid_in, const string& filename);

 private:
	void read1DMeteo(const Date&, MeteoData&); ///< No buffering
	void read1DMeteo(const Date&, MeteoData&, StationData&); ///< No buffering
	void read2DMeteo(const Date&, vector<MeteoData>&); ///< No buffering
	void read2DMeteo(const Date&, vector<MeteoData>&, vector<StationData>&); ///< No buffering

	void constructMeteo2DFilenames(const Date& date_in, vector<string>& filenames);
	void readMeteoDataLine(std::string& line, const Date& date_in, MeteoData& tmpdata, std::string filename);
	void convertUnits(MeteoData& meteo);
	void cleanup() throw();
	void read2DMeteoData(const string&, const string&, const Date&, 
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
