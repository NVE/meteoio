#ifndef __IOHANDLER_H__
#define __IOHANDLER_H__

#include "IOInterface.h"
#include "A3DIO.h"
//#include "LegacyIO.ph"
#include "IOExceptions.h"

/*#include "ConfigReader.h"
#include "StationData.h"
#include "MeteoData.h"
#include "Grid2DObject.h"
#include "Date.h"
#include "Alpine3D.h"*/

#include "marshal_meteoio.h"

parclass IOHandler;

parclass IOHandler{ // Note : No heritage here for POPC++ : a parclass cannot herit from a class
		classuid(1003);
	public:
		IOHandler(const string& configfile) @{ power=100 ?: 50; };
		~IOHandler();

		virtual void get2DGridSize(int& nx, int& ny);
		virtual void read2DGrid([out]Grid2DObject& dem_out, const string& parameter="");

		virtual void readDEM([out]Grid2DObject& dem_out);
		virtual void readLanduse([out]Grid2DObject& landuse_out);

		virtual void readMeteoData([in]const Date_IO& date_in,
				     [out, proc=marshal_vector_MeteoData]vector<MeteoData>& vecMeteo);

		virtual void readMeteoData([in]const Date_IO& date_in,
				     [out, proc=marshal_vector_MeteoData]vector<MeteoData>& vecMeteo,
				     [out, proc=marshal_vector_StationData]vector<StationData>& vecStation);

		virtual void readAssimilationData([in] const Date_IO&,[out] Grid2DObject& da_out);
		virtual void readSpecialPoints([out,proc=marshal_CSpecialPTSArray]CSpecialPTSArray& pts);

		virtual void write2DGrid([in]const Grid2DObject& grid_in, const string& filename);

	private:
		string ascii_src;
		string boschung_src;
		string imis_src;

	private:
		void cleanup();// throw();
		void loadDynamicPlugins();

		ConfigReader cfg;
		A3DIO fileio;
		DynamicLibrary* dynLibraryBoschung;
		IOInterface* boschungio;
		DynamicLibrary* dynLibraryImis;
		IOInterface* imisio;
};

#endif
