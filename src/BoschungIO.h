#ifndef __BOSCHUNGIO_H__
#define __BOSCHUNGIO_H__

#include "IOInterface.h"
#include "ConfigReader.h"
#include "IOUtils.h"
#include "MeteoBuffer.h"
#include "Meteo1DResampler.h"

#include <string>
#include <sstream>
#include <libxml++/libxml++.h>
#include <iostream>

#include "IOExceptions.h"
#include "DynamicLibrary.h"

using namespace IOUtils;

class BoschungIO : public IOInterface {
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

		virtual void readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo);
		virtual void readMeteoData(const Date_IO& date_in, vector<MeteoData>& vecMeteo, vector<StationData>& vecStation);

		virtual void readAssimilationData(const Date_IO&, Grid2DObject& da_out);
		virtual void readSpecialPoints(CSpecialPTSArray& pts);

		virtual void write2DGrid(const Array2D<double>& grid_in, const double& xllcorner, const double& yllcorner, const double& cellsize, const string& name);
		virtual void write2DGrid(const Grid2DObject& grid_in, const string& name);

		void read2DMeteo(const Date_IO&, vector<MeteoData>&); ///< No buffering
		void read2DMeteo(const Date_IO&, vector<MeteoData>&, vector<StationData>&); ///<No buffering

	private:
		void convertUnits(MeteoData& meteo);
		void createBuffer(void);
		void checkForMeteoFiles(const string& xmlpath, const string& stationname, const Date_IO& date_in,
						    string& filename_out, Date_IO& date_out);
		void xmlParseStringToDouble(const string& str_in, double& d_out, const string& parname);
		std::string xmlGetNodeContent(xmlpp::Node* pNode, const string& nodename);
		void xmlExtractData(const string& filename, const Date_IO& date_in, MeteoData& md, StationData& sd);
		std::string xmlGetNodeName(xmlpp::Node* pNode);
		xmlpp::Node* xmlGetNode(xmlpp::Node* parentNode, const string& nodename);
		void stringToDate_IO(const string& tmp, Date_IO& date_out) const;
		bool validFilename(const string& tmp) const;
		void cleanup() throw();
		void getFiles(const string& stationsname, const Date_IO& start_date, const Date_IO& end_date,
				vector<string>& vecFiles, vector<Date_IO>& vecDate_IO);
		void readStationNames(void);
		bool bufferData(const Date_IO& date_in, const unsigned int& stationnr);

		ConfigReader cfg;
		ifstream fin; //Input file streams
		vector<MeteoBuffer> unfilteredMeteoBuffer;
		vector<MeteoBuffer> filteredMeteoBuffer;
		vector<string> vecStationName;
};

#endif
