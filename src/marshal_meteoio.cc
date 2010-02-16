/***********************************************************************************/
/*  Copyright 2009 HES-SO Fribourg                                                 */
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
#ifdef _POPC_

#include "marshal_meteoio.h"

void marshal_uint(POPBuffer &buf, unsigned int &data, int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if (flag & FLAG_MARSHAL) {
		int n=(int)data;
		buf.Pack(&n,1);
	} else {
		int n;
		buf.UnPack(&n,1);
		data=(unsigned int)n;
	}
}

void marshal_slope_type(POPBuffer &buf, DEMObject::slope_type &data, int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if (flag & FLAG_MARSHAL) {
		int n=(int)data;
		buf.Pack(&n,1);
	} else {
		int n;
		buf.UnPack(&n,1);
		data=(DEMObject::slope_type)n;
	}
}

void marshal_geo_distances(POPBuffer &buf, Coords::geo_distances &data, int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if (flag & FLAG_MARSHAL) {
		int n=(int)data;
		buf.Pack(&n,1);
	} else {
		int n;
		buf.UnPack(&n,1);
		data=(Coords::geo_distances)n;
	}
}

void marshal_POINTSArray(POPBuffer &buf,POINTSArray &data, int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if (flag & FLAG_MARSHAL) {
		int n=data.size();
		buf.Pack(&n,1);
		if (n) buf.Pack((int *)((POINT *)&data[0]), 2*n);
	} else {
		int n;

		buf.UnPack(&n,1);
		data.resize(n);
		if (n) buf.UnPack((int *)((POINT *)&data[0]), 2*n);
	}
}

void marshal_METEO_DATASET(POPBuffer &buf, METEO_DATASET &data, int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if(flag&FLAG_MARSHAL) {
		int n=data.size();
		buf.Pack(&n,1);
		for(int i=0;i<n;i++) {
			data[i].Serialize(buf,true);
		}
	} else {
		int n=0;
		buf.UnPack(&n,1);
		data.clear();  
		for(int i=0;i<n;i++) {
			MeteoData obj;
			obj.Serialize(buf,false);
			data.push_back(obj);
		}
	}
}

void marshal_vector_METEO_DATASET(POPBuffer &buf, std::vector<METEO_DATASET> &data, int maxsize, int flag, POPMemspool *temp)
{
	if(flag&FLAG_MARSHAL) {
		int n=data.size();
		buf.Pack(&n,1);
		for(int i=0;i<n;i++) {
			marshal_METEO_DATASET(buf, data[i], maxsize, FLAG_MARSHAL, temp);
		}
	} else {
		int n=0;
		buf.UnPack(&n,1);
		data.clear();  
		for(int i=0;i<n;i++) {
			METEO_DATASET obj;
			marshal_METEO_DATASET(buf, obj, maxsize, !FLAG_MARSHAL, temp);
			data.push_back(obj);
		}
	}
}

void marshal_map_str_str(POPBuffer &buf, std::map<std::string, std::string> &data_map, int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if(flag&FLAG_MARSHAL) {
		int n=data_map.size();
		buf.Pack(&n,1);
		for(std::map<std::string, std::string>::const_iterator it = data_map.begin(); it != data_map.end(); ++it) {
			buf.Pack(&(it->first),1);
			buf.Pack(&(it->second),1);
		}

	} else {
		int n=0;
		std::string key;
		std::string value;
		buf.UnPack(&n,1);
		data_map.clear();
		for(int i=0;i<n;i++) {
			buf.UnPack(&key,1);
			buf.UnPack(&value,1);
			data_map[key] = value;
		}
	}
}

void marshal_vecstr(POPBuffer &buf, std::vector<std::string> &data, int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if(flag&FLAG_MARSHAL) {
		int n=data.size();
		buf.Pack(&n,1);
		for (int jj=0; jj<n; jj++) {
			buf.Pack(&(data[jj]),1);
		}
	} else {
		int n=0;
		std::string value;
		buf.UnPack(&n,1);
		data.clear();
		for(int jj=0;jj<n;jj++) {
			buf.UnPack(&value,1);
			data.push_back(value);
		}
	}
}

void marshal_map_str_vecstr(POPBuffer &buf, std::map<std::string, STR_VECTOR> &data_map, int maxsize, int flag, POPMemspool *temp)
{
	if(flag&FLAG_MARSHAL) {
		int n=data_map.size();
		buf.Pack(&n,1);
		for(std::map<std::string, STR_VECTOR>::const_iterator it = data_map.begin(); it != data_map.end(); ++it) {
			buf.Pack(&(it->first),1);
			STR_VECTOR tmp_strvec = it->second;
			marshal_vecstr(buf, tmp_strvec, maxsize, FLAG_MARSHAL, temp);
		}
	} else {
		int n=0;
		std::string key;
		std::vector<std::string> value;
		buf.UnPack(&n,1);
		data_map.clear();
		for(int i=0;i<n;i++) {
			buf.UnPack(&key,1);
			marshal_vecstr(buf, value, maxsize, !FLAG_MARSHAL, temp);
			data_map[key] = value;
		}
	}
}

void marshal_Coords(POPBuffer &buf, Coords &data, int maxsize, int flag, POPMemspool *temp) {
	(void)maxsize;
	(void)*temp;
	if(flag&FLAG_MARSHAL) {
		data.Serialize(buf,true);
	} else {
		data.Serialize(buf,false);
	}
}

void marshal_STATION_DATASET(POPBuffer &buf, STATION_DATASET &data, int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if(flag&FLAG_MARSHAL) {
		int n=data.size();
		buf.Pack(&n,1);
		for(int i=0;i<n;i++) {
			data[i].Serialize(buf,true);
		}
	} else {
		int n=0;
		buf.UnPack(&n,1);
		data.clear();
		for(int i=0;i<n;i++) {
			StationData obj;
			obj.Serialize(buf,false);
			data.push_back(obj);
		}
	}
}

void marshal_vector_STATION_DATASET(POPBuffer &buf, std::vector<STATION_DATASET> &data, int maxsize, int flag, POPMemspool *temp)
{
	if(flag&FLAG_MARSHAL) {
		int n=data.size();
		buf.Pack(&n,1);
		for(int i=0;i<n;i++) {
			marshal_STATION_DATASET(buf, data[i], maxsize, FLAG_MARSHAL, temp);
		}
	} else {
		int n=0;
		buf.UnPack(&n,1);
		data.clear();
		for(int i=0;i<n;i++) {
			STATION_DATASET obj;
			marshal_STATION_DATASET(buf, obj, maxsize, !FLAG_MARSHAL, temp);
			data.push_back(obj);
		}
	}
}

void marshal_vector_Grid2DObject(POPBuffer &buf, std::vector<Grid2DObject> &data, int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	assert(false); /* This line is here to check if the method is used*/
	if(flag&FLAG_MARSHAL) {
		int n=data.size();
		buf.Pack(&n,1);
		for(int i=0;i<n;i++) {
			data[i].Serialize(buf,true);
		}
	} else {
		int n=0;
		buf.UnPack(&n,1);
		data.clear();
		for(int i=0;i<n;i++) {
			Grid2DObject obj; 
      //buf.UnPack(&obj,1);
			obj.Serialize(buf,false);
      //marshal_Grid2DObject(buf, *obj, 0, flag, NULL);
			data.push_back(obj);
		}
	}
}

void marshal_TYPE_DOUBLE2D(POPBuffer &buf, TYPE_DOUBLE2D &data,int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if (flag & FLAG_MARSHAL) {
		unsigned int nx,ny;
		data.size(nx,ny);
		buf.Pack(&nx,1);
		buf.Pack(&ny,1);
		if (nx>0 && ny>0) {
			for (unsigned int i=0;i<nx;i++) buf.Pack(&data(i,0),ny);
		}
	} else {
		unsigned int nx,ny;
		buf.UnPack(&nx,1);
		buf.UnPack(&ny,1);
		if (nx>0 && ny>0) {
			data.resize(nx,ny);
			for (unsigned int i=0;i<nx;i++) buf.UnPack(&data(i,0),ny);
		} else
			data.clear();
	}
}

void marshal_TYPE_INT2D(POPBuffer &buf, TYPE_INT2D &data,int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if (flag & FLAG_MARSHAL) {
		unsigned int dim[2];
		data.size(dim[0],dim[1]);
		buf.Pack(dim,2);
		unsigned int nx=dim[0];
		unsigned int ny=dim[1];
		if (nx>0 && ny>0) {
			for (unsigned int i=0;i<nx;i++) buf.Pack(&data(i,0),ny);
		}
	} else {
		unsigned int dim[2];
		buf.UnPack(dim,2);
		unsigned int nx=dim[0];
		unsigned int ny=dim[1];
		if (nx>0 && ny>0) {
			data.resize(nx,ny);
			for (unsigned int i=0;i<nx;i++) buf.UnPack(&data(i,0),ny);
		} else
			data.clear();
	}
}

void marshal_CDoubleArray(POPBuffer &buf, CDoubleArray &data,int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	assert(false); /* This line is here to check if the method is used*/
	if (flag & FLAG_MARSHAL) {
		unsigned int n=data.size();
		buf.Pack(&n,1);
		if (n) buf.Pack((double *)&data[0],n);
	} else {
		unsigned int n;
		buf.UnPack(&n,1);
		if (n) {
			buf.UnPack((double *)&data[0],n);
			data.resize(n);
		} else
			data.clear();
	}

}

void marshal_CNodeArray(POPBuffer &buf,CNodeArray &data,int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	assert(false); /* This line is here to check if the method is used*/
	if (flag & FLAG_MARSHAL) {
		int n=data.size();
		buf.Pack(&n,1);
		if (n) {
			buf.Pack((double *)((NODE *)&data[0]),n*(sizeof(NODE)/sizeof(double)));
		}
	} else {
		int n;
		buf.UnPack(&n,1);
		if (n) {
			data.resize(n);
			buf.UnPack((double *)((NODE *)&data[0]),n*(sizeof(NODE)/sizeof(double)));
		} else
			data.clear();
	}
}

void marshal_update_CNodeArray(POPBuffer &buf,CNodeArray &data,int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	assert(false); /* This line is here to check if the method is used*/
	if (flag & FLAG_MARSHAL) {
		int n=data.size();
		buf.Pack(&n,1);
		NODE *tmp=&data[0];
		for (int i=0;i<n;i++,tmp++) {
			buf.Pack(&(tmp->u),1);
			buf.Pack(&(tmp->v),1);
			buf.Pack(&(tmp->w),1);
			buf.Pack(&(tmp->slope),1);
			buf.Pack(&(tmp->sl),1);

			buf.Pack(&(tmp->tet),1);
			buf.Pack(&(tmp->p),1);
			buf.Pack(&(tmp->Km),1);
			buf.Pack(&(tmp->lm),1);

			buf.Pack(&(tmp->wstar),1);
			buf.Pack(&(tmp->e),1);
		}
	} else {
		int n;
		buf.UnPack(&n,1);
		if(n>0)
			data.resize(n);
		else
			data.clear();
		NODE *tmp=&data[0];
		for (int i=0;i<n;i++,tmp++) {
			buf.UnPack(&(tmp->u),1);
			buf.UnPack(&(tmp->v),1);
			buf.UnPack(&(tmp->w),1);
			buf.UnPack(&(tmp->slope),1);
			buf.UnPack(&(tmp->sl),1);

			buf.UnPack(&(tmp->tet),1);
			buf.UnPack(&(tmp->p),1);
			buf.UnPack(&(tmp->Km),1);
			buf.UnPack(&(tmp->lm),1);

			buf.UnPack(&(tmp->wstar),1);
			buf.UnPack(&(tmp->e),1);
		}
	}
}

void marshal_input_CNodeArray(POPBuffer &buf,CNodeArray &data,int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	assert(false); /* This line is here to check if the method is used*/
	if (flag & FLAG_MARSHAL) {
		int n=data.size();
		buf.Pack(&n,1);
		NODE *tmp=&data[0];
		for (int i=0;i<n;i++,tmp++) {
			buf.Pack(&(tmp->u),1);
			buf.Pack(&(tmp->v),1);
			buf.Pack(&(tmp->w),1);
			buf.Pack(&(tmp->tet),1);
			buf.Pack(&(tmp->p),1);
			buf.Pack(&(tmp->Km),1);
			buf.Pack(&(tmp->lm),1);
			buf.Pack(&(tmp->e),1);
		}
	} else {
		int n;
		buf.UnPack(&n,1);
		if(n>0)
			data.resize(n);
		else
			data.clear();
		NODE *tmp=&data[0];
		for (int i=0;i<n;i++,tmp++) {
			buf.UnPack(&(tmp->u),1);
			buf.UnPack(&(tmp->v),1);
			buf.UnPack(&(tmp->w),1);
			buf.UnPack(&(tmp->tet),1);
			buf.UnPack(&(tmp->p),1);
			buf.UnPack(&(tmp->Km),1);
			buf.UnPack(&(tmp->lm),1);
			buf.UnPack(&(tmp->e),1);
		}
	}
}

#endif
