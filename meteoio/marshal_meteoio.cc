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

#include <meteoio/marshal_meteoio.h>

namespace mio {

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

void marshal_MeteoParameters(POPBuffer &buf, MeteoData::Parameters &data, int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if (flag & FLAG_MARSHAL) {
		int n=(int)data;
		buf.Pack(&n,1);
	} else {
		int n;
		buf.UnPack(&n,1);
		data=(MeteoData::Parameters)n;
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

void marshal_vec_coords(POPBuffer &buf,std::vector<Coords> &data, int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if (flag & FLAG_MARSHAL) {
		int n=data.size();
		buf.Pack(&n,1);
		for(int i=0;i<n;i++) {
			marshal_Coords(buf, data[i], maxsize, FLAG_MARSHAL, temp);
		}
	} else {
		int n=0;
		buf.UnPack(&n,1);
		data.clear();
		for(int i=0;i<n;i++) {
			Coords obj;
			marshal_Coords(buf, obj, maxsize, !FLAG_MARSHAL, temp);
			data.push_back(obj);
		}
	}
}

void marshal_METEO_TIMESERIE(POPBuffer &buf, METEO_TIMESERIE &data, int maxsize, int flag, POPMemspool *temp)
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

void marshal_vector_METEO_TIMESERIE(POPBuffer &buf, std::vector<METEO_TIMESERIE> &data, int maxsize, int flag, POPMemspool *temp)
{
	if(flag&FLAG_MARSHAL) {
		int n=data.size();
		buf.Pack(&n,1);
		for(int i=0;i<n;i++) {
			marshal_METEO_TIMESERIE(buf, data[i], maxsize, FLAG_MARSHAL, temp);
		}
	} else {
		int n=0;
		buf.UnPack(&n,1);
		data.clear();
		for(int i=0;i<n;i++) {
			METEO_TIMESERIE obj;
			marshal_METEO_TIMESERIE(buf, obj, maxsize, !FLAG_MARSHAL, temp);
			data.push_back(obj);
		}
	}
}

void marshal_map_str_dbl(POPBuffer &buf, std::map<std::string, double> &data_map, int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if(flag&FLAG_MARSHAL) {
		int n=data_map.size();
		buf.Pack(&n,1);
		for(std::map<std::string, double>::const_iterator it = data_map.begin(); it != data_map.end(); ++it) {
			buf.Pack(&(it->first),1);
			buf.Pack(&(it->second),1);
		}

	} else {
		int n=0;
		std::string key;
		double value;
		buf.UnPack(&n,1);
		data_map.clear();
		for(int i=0;i<n;i++) {
			buf.UnPack(&key,1);
			buf.UnPack(&value,1);
			data_map[key] = value;
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

//HACK is this still needed?
void marshal_Coords(POPBuffer &buf, Coords &data, int maxsize, int flag, POPMemspool *temp) {
	(void)maxsize;
	(void)*temp;
	if(flag&FLAG_MARSHAL) {
		data.Serialize(buf,true);
	} else {
		data.Serialize(buf,false);
	}
}

void marshal_STATION_TIMESERIE(POPBuffer &buf, STATION_TIMESERIE &data, int maxsize, int flag, POPMemspool *temp)
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

void marshal_vector_STATION_TIMESERIE(POPBuffer &buf, std::vector<STATION_TIMESERIE> &data, int maxsize, int flag, POPMemspool *temp)
{
	if(flag&FLAG_MARSHAL) {
		int n=data.size();
		buf.Pack(&n,1);
		for(int i=0;i<n;i++) {
			marshal_STATION_TIMESERIE(buf, data[i], maxsize, FLAG_MARSHAL, temp);
		}
	} else {
		int n=0;
		buf.UnPack(&n,1);
		data.clear();
		for(int i=0;i<n;i++) {
			STATION_TIMESERIE obj;
			marshal_STATION_TIMESERIE(buf, obj, maxsize, !FLAG_MARSHAL, temp);
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

void marshal_DOUBLE2D(POPBuffer &buf, DOUBLE2D &data, int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if (flag & FLAG_MARSHAL) {
		unsigned int nx,ny;
		data.size(nx,ny);
		buf.Pack(&nx,1);
		buf.Pack(&ny,1);
		bool keep_nodata = data.getKeepNodata();
		buf.Pack(&keep_nodata,1);
		if (nx>0 && ny>0) {
			const unsigned int nxy=nx*ny;
			buf.Pack(&data(0,0), nxy);
		}
	} else {
		unsigned int nx,ny;
		buf.UnPack(&nx,1);
		buf.UnPack(&ny,1);
		bool keep_nodata;
		buf.UnPack(&keep_nodata,1);
		data.setKeepNodata(keep_nodata);
		if (nx>0 && ny>0) {
			data.resize(nx,ny);
			const unsigned int nxy=nx*ny;
			buf.UnPack(&data(0,0),nxy);
		} else
			data.clear();
	}
}

void marshal_DOUBLE3D(POPBuffer &buf, DOUBLE3D &data, int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if (flag & FLAG_MARSHAL) {
		unsigned int nx,ny,nz;
		data.size(nx,ny,nz);
		buf.Pack(&nx,1);
		buf.Pack(&ny,1);
		buf.Pack(&nz,1);
		bool keep_nodata = data.getKeepNodata();
		buf.Pack(&keep_nodata,1);
		if (nx>0 && ny>0 && nz>0) {
			const unsigned int nxyz=nx*ny*nz;
			buf.Pack(&data(0,0,0),nxyz);
		}
	} else {
		unsigned int nx,ny,nz;
		buf.UnPack(&nx,1);
		buf.UnPack(&ny,1);
		buf.UnPack(&nz,1);
		bool keep_nodata;
		buf.UnPack(&keep_nodata,1);
		data.setKeepNodata(keep_nodata);
		if (nx>0 && ny>0 && nz>0) {
			data.resize(nx,ny,nz);
			const unsigned int nxyz=nx*ny*nz;
			buf.UnPack(&data(0,0,0),nxyz);
		} else
			data.clear();
	}
}

void marshal_INT2D(POPBuffer &buf, INT2D &data, int maxsize, int flag, POPMemspool *temp)
{
	(void)maxsize;
	(void)*temp;
	if (flag & FLAG_MARSHAL) {
		unsigned int nx, ny;
		data.size(nx,ny);
		buf.Pack(&nx,1);
		buf.Pack(&ny,1);
		bool keep_nodata = data.getKeepNodata();
		buf.Pack(&keep_nodata,1);
		if (nx>0 && ny>0) {
			const unsigned int nxy=nx*ny;
			buf.Pack(&data(0,0), nxy);
		}
	} else {
		unsigned int nx,ny;
		buf.UnPack(&nx,1);
		buf.UnPack(&ny,1);
		bool keep_nodata;
		buf.UnPack(&keep_nodata,1);
		data.setKeepNodata(keep_nodata);
		if (nx>0 && ny>0) {
			data.resize(nx,ny);
			const unsigned int nxy=nx*ny;
			buf.UnPack(&data(0,0), nxy);
		} else
			data.clear();
	}
}

} //namespace

#endif
