#ifndef __BOSCHUNGIO_H__
#define __BOSCHUNGIO_H__

#include "IOHandler.h"
#include "ConfigReader.h"
#include "slfutils.h"
#include "MeteoBuffer.h"
#include "Meteo1DResampler.h"

#include <string>
#include <sstream>
#include <libxml++/libxml++.h>
#include <iostream>

#include "slfexceptions.h"
#include "DynamicLibrary.h"

class BoschungIO : public IOHandler {
 public:
	//virtual BoschungIO* clone() const;

	BoschungIO(void (*delObj)(void*), const string& filename);

	BoschungIO(const std::string& configfile);
	BoschungIO(const BoschungIO&);
	BoschungIO(const ConfigReader&);
	~BoschungIO() throw();

	virtual void get2DGridSize(int& nx, int& ny);
	virtual void read2DGrid(Grid2DObject& dem_out, const string& parameter="");

	virtual void readDEM(Grid2DObject& dem_out);
	virtual void readLanduse(Grid2DObject& landuse_out);

	virtual void readMeteoData(const Date& date_in, vector<MeteoData>& vecMeteo);
	virtual void readMeteoData(const Date& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation);

	virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
	virtual void readSpecialPoints(CSpecialPTSArray& pts);

	virtual void write2DGrid(const Grid2DObject& grid_in, const string& filename);

	void read2DMeteo(const Date&, vector<MeteoData>&); ///< No buffering
	void read2DMeteo(const Date&, vector<MeteoData>&, vector<StationData>&); ///<No buffering

 private:
	void convertUnits(MeteoData& meteo);
	void createBuffer(void);
	void checkForMeteoFiles(const string& xmlpath, const string& stationname, const Date& date_in,
					    string& filename_out, Date& date_out);
	void xmlParseStringToDouble(const string& str_in, double& d_out, const string& parname);
	std::string xmlGetNodeContent(xmlpp::Node* pNode, const string& nodename);
	void xmlExtractData(const string& filename, const Date& date_in, MeteoData& md, StationData& sd);
	std::string xmlGetNodeName(xmlpp::Node* pNode);
	xmlpp::Node* xmlGetNode(xmlpp::Node* parentNode, const string& nodename);
	void stringToDate(const string& tmp, Date& date_out) const;
	bool validFilename(const string& tmp) const;
	void cleanup() throw();
	void getFiles(const string& stationsname, const Date& start_date, const Date& end_date, vector<string>& vecFiles, vector<Date>& vecDate);
	void readStationNames(void);
	bool bufferData(const Date& date_in, const unsigned int& stationnr);

	ConfigReader cfg;
	ifstream fin; //Input file streams
	vector<MeteoBuffer> unfilteredMeteoBuffer;
	vector<MeteoBuffer> filteredMeteoBuffer;
	vector<string> vecStationName;
};

#endif
