#ifndef __IOHANDLER_H__
#define __IOHANDLER_H__

#include "IOInterface.h"
#include "A3DIO.h"
#include "IOExceptions.h"

#include <map>

#include "marshal_meteoio.h"

parclass IOHandler;

parclass IOPlugin {
	public:
		std::string libname;
		std::string classname;
		IOInterface *io;
		DynamicLibrary *dynLibrary;
		
		IOPlugin(std::string _s1, std::string _s2, IOInterface *p1, DynamicLibrary *p2) : libname(_s1), classname(_s2), io(p1), dynLibrary(p2){}
		IOPlugin() : libname(""), classname(""), io(NULL), dynLibrary(NULL){}
};

parclass IOHandler {
// Note : No heritage here for POPC++ : a parclass cannot herit from a class
		classuid(1003);
	public:
		// virtual IOHandler* clone() const; // lwk : not used yet
		IOHandler(const std::string& configfile) @{ power=100 ?: 50; };
		IOHandler(const IOHandler&);
		IOHandler(const ConfigReader&); //Pbl with stream in ConfigReader
		//~IOHandler() throw(); //pbl with throw
		~IOHandler();

		virtual void read2DGrid([out]Grid2DObject& dem_out, const string& parameter="");

		virtual void readDEM([out]Grid2DObject& dem_out);
		virtual void readLanduse([out]Grid2DObject& landuse_out);

		//Pbl with vector of vector...
		virtual void readMeteoData([in]const Date_IO& dateStart, 
			     			     [in]const Date_IO& dateEnd,
			     			     std::vector< std::vector<MeteoData> >& vecMeteo, 
						     std::vector< std::vector<StationData> >& vecStation,
						     const unsigned int& stationindex=IOUtils::npos);

		virtual void readAssimilationData([in] const Date_IO&,[out] Grid2DObject& da_out);
		virtual void readSpecialPoints([out,proc=marshal_CSpecialPTSArray]CSpecialPTSArray& pts);

		virtual void write2DGrid([in]const Grid2DObject& grid_in, [in]const string& name);

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
