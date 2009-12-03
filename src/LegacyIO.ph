#ifndef LEGACYIO_H
#define LEGACYIO_H


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "Array.h"
#include "Array2D.h"


#define ERROR_GRIDSIZE long(1)
#define ERROR_STATION_INDEX long(2)
#define ERROR_INIT_STATION long(3)
#define ERROR_SIZE long(4)
#define ERROR_SNOWPACK_PARAM long(5)
#define ERROR_SNOWPACK_COMPUTE long(6)
#define ERROR_SNOWDRIFT_COMPUTE long(7)
#define ERROR_EB_SUN  long(8)
#define ERROR_EB_COMPUTE long(9)

#define ERROR_INPUT_GRIDSIZE long(10)
#define ERROR_INPUT_GRIDPOINTS long(11)
#define ERROR_INPUT_METEODATA long(12)
#define ERROR_INPUT_LANDUSE long(13)
#define ERROR_INPUT_SNOWCOVER long(14)
#define ERROR_INPUT_GRIDDATA long(15)
#define ERROR_INPUT_DEM long(16)
#define ERROR_INPUT_METEO2D long(17)

#define ERROR_DATA_CONSISTENCY long(18)
#define ERROR_INPUT_RADIATION long(19)

#define MAX_STRING_LENGTH 256
#define MAX_LINE_LENGTH 6000

typedef struct {
	double x;
	double y;
	double z;
	double u;
	double v;
	double w;

	double slope;	/* In wind direction [-PI/2,PI/2]*/
	double sl;	/* General (maximum) slope angle */
	double tet;
	double p;
	double Km;
	double lm;
	double wstar;
	double e;
	double c;
	double azi;	/* Slope Azimut */
	double sx;	/* x -component of normal on Surface element */
	double sy;	/* y -component of normal on Surface element */
	
	double rh;		//subl	
	double wnd;		//subl
	double repRadius; //subl
	double initMassChange; //subl
	double initSubl;//subl
	double subl;
	/*double c_new;*/
	double specHumidity;//subl
	double iniSpecHum;
} NODE;

//Enumerate the slope shapes
enum {Flat, Luff, Lee, N_SLOPE};

typedef Array<NODE> CNodeArray;

typedef Array2D<int> CElementArray;

typedef Array<double> CDoubleArray;

typedef Array<int> CIntArray;

// marshal_meteoio must be included after the structures'definitions
#include "marshal_meteoio.h"

parclass LegacyIO
{
	public:
		LegacyIO( [in] const std::string &meteopath) @{od.url("localhost");};
		~LegacyIO();

		void GetGridSize([in] const std::string& grid_name, [out] int &nx, [out] int &ny, [out] int &nz);
		void GetGridPoints([in] const std::string& grid_name, [out, proc=marshal_CDoubleArray] CDoubleArray &x, [out, proc=marshal_CDoubleArray]  CDoubleArray &y, [out, proc=marshal_CDoubleArray]  CDoubleArray &z);
		void GetGridData([out, proc=marshal_input_CNodeArray] CNodeArray &data, [in] const std::string& hour);

		classuid(1002);

	private:
		void moveToMarker(FILE *fp, const std::string& file_name, const std::string& marker);
		std::string meteopathname;
		int dimx, dimy, dimz;

};

#endif
