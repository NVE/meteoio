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

		virtual void read2DGrid([out]Grid2DObject& dem_out, const string& parameter="");

		virtual void readDEM([out]Grid2DObject& dem_out);
		virtual void readLanduse([out]Grid2DObject& landuse_out);

		virtual void readMeteoData([in]const Date_IO& dateStart, 
			     			     [in]const Date_IO& dateEnd,
			     			     std::vector< std::vector<MeteoData> >& vecMeteo, 
						     std::vector< std::vector<StationData> >& vecStation,
						     const unsigned int& stationindex=IOUtils::npos);
		/*				     
		virtual void readMeteoData([in]const Date_IO& date_in,
				     [out, proc=marshal_vector_MeteoData]vector<MeteoData>& vecMeteo,
				     [out, proc=marshal_vector_StationData]vector<StationData>& vecStation);
		*/
		virtual void readAssimilationData([in] const Date_IO&,[out] Grid2DObject& da_out);
		virtual void readSpecialPoints([out,proc=marshal_CSpecialPTSArray]CSpecialPTSArray& pts);

		//BUG: popc does not see that TYPE_DOUBLE2D and Array2D are the same... and it does not like the Array2D<double> either...
		//virtual void write2DGrid([in, proc=marshal_TYPE_DOUBLE2D]const Array2D<double>& grid_in, [in]const double& xllcorner, [in]const double& yllcorner, [in]const double& cellsize, [in]const string& name);
		virtual void write2DGrid([in]const Grid2DObject& grid_in, [in]const string& name);

	private:
		string ascii_src;
		string boschung_src;
		string imis_src;

	private:
		void loadDynamicPlugins();

		ConfigReader cfg;
		A3DIO fileio;
		DynamicLibrary* dynLibraryBoschung;
		IOInterface* boschungio;
		DynamicLibrary* dynLibraryImis;
		IOInterface* imisio;
};

#endif
