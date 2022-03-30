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
	
	typedef struct MYSQL_FIELD {
		MYSQL_FIELD() : str(""), dt(), str_len(0), val(mio::IOUtils::nodata), MysqlType(MYSQL_TYPE_NULL) {}
		MYSQL_FIELD(const enum_field_types &type) : str(""), dt(), str_len(0), val(mio::IOUtils::nodata), MysqlType(type) {}
		MYSQL_FIELD(const std::string& i_str) : str(""), dt(), str_len(0), val(mio::IOUtils::nodata), MysqlType(MYSQL_TYPE_STRING) {setString(i_str);}
		
		void reset() {str[0]='\0'; dt.setUndef(); str_len=0; val=mio::IOUtils::nodata; MysqlType=MYSQL_TYPE_NULL;}
		void setString(const std::string& i_str) {reset(); strncpy(str, i_str.c_str(), std::min(static_cast<int>(i_str.size()), STRING_SIZE)); str_len=strlen(str); MysqlType=MYSQL_TYPE_STRING;}
		void setDate(const mio::Date& i_dt) {reset(); dt=i_dt; MysqlType=MYSQL_TYPE_DATETIME;}
		void setDouble(const double& i_val) {reset(); val=i_val; MysqlType=MYSQL_TYPE_DOUBLE;}
		
		char str[STRING_SIZE];
		mio::Date dt;
		long unsigned int str_len;
		double val;
		enum_field_types MysqlType; //see mysql/field_types.h
	} fType;
	
	MYSQL* initMysql(const std::string& mysqlhost, const std::string& mysqluser, const std::string& mysqlpass, const std::string& mysqldb);
	MYSQL_STMT* initStmt(MYSQL **mysql, const std::string& query, const long unsigned int& ref_param_count);
	
	void bindParams(MYSQL_STMT **stmt, std::vector<fType> &params_fields);
	void bindResults(MYSQL_STMT **stmt, std::vector<fType> &result_fields);
}

#endif
