/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
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

namespace mio {

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
	double iniTet;
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
	double iniSpecHum; //subl
	double fraction; //subl
	
} NODE;

//Enumerate the slope shapes
enum {Flat, Luff, Lee, N_SLOPE};

typedef Array<NODE> CNodeArray;

typedef Array2D<int> CElementArray;

typedef Array<double> CDoubleArray;

typedef Array<int> CIntArray;

class LegacyIO
{
	public:
		LegacyIO(const std::string &meteopath);
		~LegacyIO();

		void GetGridSize(const std::string& grid_name, int &nx, int &ny, int &nz);
		void GetGridPoints(const std::string& grid_name, CDoubleArray &x, CDoubleArray &y, CDoubleArray &z);
		void GetGridData(CNodeArray &data, const std::string& hour);

	private:
		void moveToMarker(FILE *fp, const std::string& file_name, const std::string& marker);
		std::string meteopathname;
		int dimx, dimy, dimz;
};
} //end namespace

#endif
