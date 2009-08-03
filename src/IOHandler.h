#ifndef __IOHANDLER_H__
#define __IOHANDLER_H__

#ifdef _POPC_
#error
#endif

#include "IOInterface.h"
#include "A3DIO.h"
#include "IOExceptions.h"

class IOHandler : public IOInterface {
	public:
		// virtual IOHandler* clone() const; // lwk : not used yet

		IOHandler(const std::string& configfile);
		IOHandler(const IOHandler&);
		IOHandler(const ConfigReader&);
		~IOHandler() throw();

		virtual void read2DGrid(Grid2DObject& dem_out, const string& parameter="");

		virtual void readDEM(Grid2DObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);

		virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  unsigned int stationindex=IOUtils::npos);

		virtual void readAssimilationData(const Date_IO&, Grid2DObject& da_out);
		virtual void readSpecialPoints(CSpecialPTSArray& pts);

		virtual void write2DGrid(const Grid2DObject& grid_in, const string& name);

	private:
		string ascii_src;
		string boschung_src;
		string imis_src;
		string geotop_src;

	private:
		void cleanup() throw();
		void loadDynamicPlugins();
		void loadPlugin(const std::string& libname, const std::string& classname, 
					 DynamicLibrary*& dynLibrary, IOInterface*& io);
		void deletePlugin(DynamicLibrary*& dynLibrary, IOInterface*& io) throw();

		ConfigReader cfg;
		A3DIO fileio;
		DynamicLibrary* dynLibraryBoschung;
		IOInterface* boschungio;
		DynamicLibrary* dynLibraryImis;
		IOInterface* imisio;
		DynamicLibrary* dynLibraryGeoTOP;
		IOInterface* geotopio;
};

#endif
