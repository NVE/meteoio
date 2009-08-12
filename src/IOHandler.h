#ifndef __IOHANDLER_H__
#define __IOHANDLER_H__

#ifdef _POPC_
#error
#endif

#include "IOInterface.h"
#include "A3DIO.h"
#include "IOExceptions.h"

#include <map>

typedef std::map<string, string> MAP_STR_STR;

class IOPlugin {
 public:
	std::string libname;
	std::string classname;
	IOInterface *io;
	DynamicLibrary *dynLibrary;
	
	IOPlugin(std::string _s1, std::string _s2, IOInterface *p1, DynamicLibrary *p2) : libname(_s1), classname(_s2), io(p1), dynLibrary(p2){}
	IOPlugin() : libname(""), classname(""), io(NULL), dynLibrary(NULL){}
};

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
							  std::vector<METEO_DATASET>& vecMeteo, 
							  std::vector<STATION_DATASET>& vecStation,
							  const unsigned int& stationindex=IOUtils::npos);

		virtual void readAssimilationData(const Date_IO&, Grid2DObject& da_out);
		virtual void readSpecialPoints(CSpecialPTSArray& pts);

		virtual void write2DGrid(const Grid2DObject& grid_in, const string& name);

	private:
		string ascii_src;
		string boschung_src;
		string imis_src;
		string geotop_src;

	private:
		void loadDynamicPlugins();
		void loadPlugin(const std::string& libname, const std::string& classname, 
					 DynamicLibrary*& dynLibrary, IOInterface*& io);
		void deletePlugin(DynamicLibrary*& dynLibrary, IOInterface*& io) throw();
		void registerPlugins();
		IOInterface *getPlugin(const std::string&);

		ConfigReader cfg;
		map<std::string, IOPlugin> mapPlugins;
		map<std::string, IOPlugin>::iterator mapit;
		A3DIO fileio;
};

#endif
