#ifndef LEGACYIO_H
#define LEGACYIO_H

#ifdef _POPC_
#error
#endif

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


/*---------------------------------------------------------------+                                                                          
 | Define Data Structures                                        |                                                                          
 +---------------------------------------------------------------*/
/*---------------------------------------------------------------+                                                                          
 | Nodes                                                         |                                                                          
 +---------------------------------------------------------------*/
typedef struct {
	double x;
	double y;
	double z;
	double u;
	double v;
	double w;

	double slope; /* In wind direction [-PI/2,PI/2]*/
	double sl;    /* General (maximum) slope angle */
	double tet;
	double p;
	double Km;
	double lm;
	double wstar;
	double e;
	double c;
	double azi;  /* Slope Azimut */
	double sx;   /* x -component of normal on Surface element */
	double sy;   /* y -component of normal on Surface element */
} NODE;

/*---------------------------------------------------------------+                                                                          
 | Enumerate the slope shapes                                    |                                                                          
 +---------------------------------------------------------------*/
enum {Flat, Luff, Lee, N_SLOPE};

typedef CArray<NODE> CNodeArray;

//typedef CArray<int[8]> CElementArray;                                                                                                     
typedef Array2D<int> CElementArray;

typedef CArray<double> CDoubleArray;

typedef CArray<int> CIntArray;


class LegacyIO
{
	public:
		LegacyIO(char *meteopath);
		~LegacyIO();

		void GetGridSize(int &nx, int &ny, int &nz);
		void GetGridPoints(CDoubleArray &x, CDoubleArray &y, CDoubleArray &z);
		void GetGridData(CNodeArray &data, char *hour);
		void PrepareNextWindField(char *hour);

	private:
		char meteopathname[MAX_STRING_LENGTH];
		int dimx, dimy, dimz;
};


#endif
