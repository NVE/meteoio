#ifndef MARSHAL_ALPINE3D_H
#define MARSHAL_ALPINE3D_H
//#include "paroc_buffer.h"

//#include "Alpine3D.h"
//#include "Constants.h"
//#include "Snowpack.h"
//#include "DriftData.h"

#include "Grid2DObject.h"
#include "StationData.h"
#include "MeteoData.h"
#include <vector>
#include "LegacyIO.ph"

typedef CArray2D<double> TYPE_DOUBLE2D;

void marshal_TYPE_DOUBLE2D(paroc_buffer &buf, TYPE_DOUBLE2D &data,int maxsize, int flag, paroc_memspool *temp);

void marshal_CDoubleArray(paroc_buffer &buf, CDoubleArray &data,int maxsize, int flag, paroc_memspool *temp);

void marshal_CNodeArray(paroc_buffer &buf,CNodeArray &data,int maxsize, int flag, paroc_memspool *temp);

void marshal_update_CNodeArray(paroc_buffer &buf,CNodeArray &data, int maxsize, int flag, paroc_memspool *temp);

void marshal_input_CNodeArray(paroc_buffer &buf,CNodeArray &data, int maxsize, int flag, paroc_memspool *temp);

void marshal_CSpecialPTSArray(paroc_buffer &buf,CSpecialPTSArray &data, int maxsize, int flag, paroc_memspool *temp);

void marshal_vector_MeteoData(paroc_buffer &buf, vector<MeteoData> &data, int maxsize, int flag, paroc_memspool *temp);

void marshal_vector_StationData(paroc_buffer &buf, vector<StationData> &data, int maxsize, int flag, paroc_memspool *temp);

void marshal_vector_Grid2DObject(paroc_buffer &buf, vector<Grid2DObject> &data, int maxsize, int flag, paroc_memspool *temp);

#endif
