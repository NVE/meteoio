#ifndef __ARCIO_H__
#define __ARCIO_H__

#include "IOInterface.h"
#include "ConfigReader.h"
#include "IOUtils.h"

#include <string>
#include <sstream>
#include <iostream>

#include "IOExceptions.h"
#include "DynamicLibrary.h"

using namespace IOUtils;

class ARCIO : public IOInterface {
	public:
		ARCIO(void (*delObj)(void*), const string& filename);

		ARCIO(const std::string& configfile);
		ARCIO(const ARCIO&);
		ARCIO(const ConfigReader&);
		~ARCIO() throw();

		virtual void read2DGrid(Grid2DObject& dem_out, const string& parameter="");

		virtual void readDEM(Grid2DObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);

		virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  const unsigned int& stationindex=IOUtils::npos);

		virtual void readAssimilationData(const Date_IO&, Grid2DObject& da_out);
		virtual void readSpecialPoints(CSpecialPTSArray& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const string& filename);

	private:
		void cleanup() throw();
		ConfigReader cfg;
		ifstream fin; //Input file streams
		ofstream fout;//Output file streams
};

#endif
