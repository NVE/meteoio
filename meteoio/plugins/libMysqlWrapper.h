// SPDX-License-Identifier: LGPL-3.0-or-later
/***********************************************************************************/
/*  Copyright 2022 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef LIBMYSQLWRAPPER_H
#define LIBMYSQLWRAPPER_H

#ifdef _WIN32
	#include <winsock.h>
#endif // _WIN32

#include <mysql.h>
#include <cstring>
#include <meteoio/dataClasses/Date.h>
#include <meteoio/IOUtils.h>

namespace mysql_wrp {
	static const int STRING_SIZE = 50;
	
	#if LIBMYSQL_VERSION_ID > 80001
		typedef bool BOOL_TYPE;
	#else
		typedef my_bool BOOL_TYPE;
	#endif
	
	enum MysqlOptions {
		COMPRESSION	= 1,
		ENCRYPTION	= 1 << 1
	};
	
	enum unitsConversions {
		NONE=0,
		C_TO_K,
		CM_TO_M,
		NORMALIZE_PC
	};
	
	typedef struct MYSQL_FIELD {
		MYSQL_FIELD() : str(""), dt(), str_len(0), buffer_len(0), val(mio::IOUtils::nodata), is_null(false), error(false), MysqlType(MYSQL_TYPE_NULL) {}
		MYSQL_FIELD(const enum_field_types &type) : str(""), dt(), str_len(0), buffer_len(0), val(mio::IOUtils::nodata), is_null(false), error(false), MysqlType(type) {}
		MYSQL_FIELD(const std::string& i_str) : str(""), dt(), str_len(0), buffer_len(0), val(mio::IOUtils::nodata), is_null(false), error(false), MysqlType(MYSQL_TYPE_STRING) {setString(i_str);}
		MYSQL_FIELD(const mio::Date& i_dt) : str(""), dt(), str_len(0), buffer_len(0), val(mio::IOUtils::nodata), is_null(false), error(false), MysqlType(MYSQL_TYPE_DATETIME) {setFromDate(i_dt, dt);}
		
		void resetDate() {dt.year=0; dt.month=0; dt.day=0; dt.hour=0; dt.minute=0; dt.second=0; dt.second_part=0;}
		void reset() {str[0]='\0'; resetDate(); str_len=0; val=mio::IOUtils::nodata; MysqlType=MYSQL_TYPE_NULL;}
		void setFromDate(const mio::Date& i_dt, MYSQL_TIME &ts) { int year, month, day, hour, minute; double second;
			i_dt.getDate(year, month, day, hour, minute, second);
			ts.year = static_cast<unsigned int>( year );
			ts.month = static_cast<unsigned int>( month );
			ts.day = static_cast<unsigned int>( day );
			ts.hour = static_cast<unsigned int>( hour );
			ts.minute = static_cast<unsigned int>( minute );
			ts.second = static_cast<unsigned int>( floor(second) );
			ts.second_part = static_cast<unsigned long int>( floor((second - floor(second))*1e6) );
		}
		void setString(const std::string& i_str) {reset(); strncpy(str, i_str.c_str(), std::min(static_cast<int>(i_str.size()), STRING_SIZE)); str_len=strlen(str); MysqlType=MYSQL_TYPE_STRING;}
		void setDate(const mio::Date& i_dt) {reset(); setFromDate(i_dt, dt); MysqlType=MYSQL_TYPE_DATETIME;}
		void setDouble(const double& i_val) {reset(); val=i_val; MysqlType=MYSQL_TYPE_DOUBLE;}
		mio::Date getDate(const double& TZ) const {mio::Date o_dt(static_cast<int>(dt.year), static_cast<int>(dt.month), static_cast<int>(dt.day), static_cast<int>(dt.hour), static_cast<int>(dt.minute), static_cast<double>(dt.second)+static_cast<double>(dt.second_part)*1e-6, TZ); return o_dt;}
		
		char str[STRING_SIZE];
		MYSQL_TIME dt;
		unsigned long int str_len;
		unsigned long int buffer_len; //reserved for mysql
		double val;
		BOOL_TYPE is_null, error;
		enum_field_types MysqlType; //see mysql/field_types.h
	} fType;
	
	MYSQL* initMysql(const std::string& mysqlhost, const std::string& mysqluser, const std::string& mysqlpass, const std::string& mysqldb, const unsigned int& options=0);
	MYSQL_STMT* initStmt(MYSQL **mysql, const std::string& query, const long unsigned int& ref_param_count);
	
	void bindParams(MYSQL_STMT **stmt, std::vector<fType> &params_fields);
	void bindResults(MYSQL_STMT **stmt, std::vector<fType> &result_fields);
	double retrieveData(const fType &field, const unsigned int& conversion=NONE);
}

#endif
