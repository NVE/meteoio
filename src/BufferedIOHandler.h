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
 * @brief 
 *
 * @author Thomas Egger
 * @date   2009-07-25
 */
class BufferedIOHandler : public IOInterface {
 public:

	BufferedIOHandler(IOInterface& _iohandler, const ConfigReader& _cfg);
	~BufferedIOHandler() throw();

	void getNextMeteoData(const Date_IO& _date, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation);
	void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, vector< vector<MeteoData> >& vecMeteo);
	void readMeteoData(const Date_IO& _date, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation);
	void clearBuffer();
	void enableResampling(const bool& _enable);

	virtual void read2DGrid(Grid2DObject& grid_out, const string& parameter="");
	virtual void readDEM(Grid2DObject& dem_out);
	virtual void readAssimilationData(const Date_IO& date_in, Grid2DObject& da_out);
	virtual void readLanduse(Grid2DObject& landuse_out);
	virtual void readSpecialPoints(CSpecialPTSArray& pts);
	virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
						  std::vector< std::vector<MeteoData> >& vecMeteo, 
						  std::vector< std::vector<StationData> >& vecStation,
						  unsigned int stationindex=IOUtils::npos);

	virtual void write2DGrid(const Grid2DObject& grid_in, const string& options="");


	static const unsigned int npos = (unsigned int)-1;             ///<npos is the out-of-range value

 private:
	unsigned int seek(const Date_IO& date_in, std::vector<MeteoData>& vecM);
	bool bufferData(const Date_IO& _date, const unsigned int& stationindex);
	void bufferAllData(const Date_IO& _date);

	IOInterface& iohandler;
	ConfigReader cfg;
	vector< vector<MeteoData> > meteoBuffer;
	vector< vector<StationData> > stationBuffer;
	std::map<std::string, Grid2DObject> mapBufferedGrids;
	bool resample;
};

#endif
