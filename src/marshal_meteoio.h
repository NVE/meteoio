#ifndef MARSHAL_METEOIO_H
#define MARSHAL_METEOIO_H

#include "Grid2DObject.h"
#include "StationData.h"
#include "MeteoData.h"
#include <vector>
#include "LegacyIO.ph"

typedef Array2D<double> TYPE_DOUBLE2D;
typedef Array2D<int> TYPE_INT2D;

void marshal_TYPE_DOUBLE2D(POPBuffer &buf, TYPE_DOUBLE2D &data,int maxsize, int flag, POPMemspool *temp);

void marshal_TYPE_INT2D(POPBuffer &buf, TYPE_INT2D &data,int maxsize, int flag, POPMemspool *temp);

void marshal_CDoubleArray(POPBuffer &buf, CDoubleArray &data,int maxsize, int flag, POPMemspool *temp);

void marshal_CNodeArray(POPBuffer &buf,CNodeArray &data,int maxsize, int flag, POPMemspool *temp);

void marshal_update_CNodeArray(POPBuffer &buf,CNodeArray &data, int maxsize, int flag, POPMemspool *temp);

void marshal_input_CNodeArray(POPBuffer &buf,CNodeArray &data, int maxsize, int flag, POPMemspool *temp);

void marshal_CSpecialPTSArray(POPBuffer &buf,CSpecialPTSArray &data, int maxsize, int flag, POPMemspool *temp);

void marshal_METEO_DATASET(POPBuffer &buf, METEO_DATASET &data, int maxsize, int flag, POPMemspool *temp);

void marshal_vector_METEO_DATASET(POPBuffer &buf, std::vector<METEO_DATASET> &data, int maxsize, int flag, POPMemspool *temp);

void marshal_STATION_DATASET(POPBuffer &buf, STATION_DATASET &data, int maxsize, int flag, POPMemspool *temp);

void marshal_vector_STATION_DATASET(POPBuffer &buf, std::vector<STATION_DATASET> &data, int maxsize, int flag, POPMemspool *temp);

void marshal_vector_Grid2DObject(POPBuffer &buf, vector<Grid2DObject> &data, int maxsize, int flag, POPMemspool *temp);

#endif
