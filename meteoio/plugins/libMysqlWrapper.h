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

class SQL_FIELD {
	public:
		static const int STRING_SIZE = 50;
		
		#if LIBMYSQL_VERSION_ID > 80001
			typedef bool BOOL_TYPE;
		#else
			typedef my_bool BOOL_TYPE;
		#endif
		
		enum unitsConversions {
			NONE=0,
			C_TO_K,
			CM_TO_M,
			NORMALIZE_PC
		};
	
		SQL_FIELD();
		SQL_FIELD(const enum_field_types &type);
		SQL_FIELD(const std::string& i_str);
		SQL_FIELD(const mio::Date& i_dt);
		
		void resetDate();
		void reset();
		void setFromDate(const mio::Date& i_dt, MYSQL_TIME &ts);
		void setString(const std::string& i_str);
		void setDate(const mio::Date& i_dt);
		void setDouble(const double& i_val);
		mio::Date getDate(const double& TZ) const;
		
		//const std::string param;
		char str[STRING_SIZE];
		MYSQL_TIME dt;
		unsigned long int str_len;
		unsigned long int buffer_len; //reserved for mysql
		double val;
		//const unsigned int processing;
		BOOL_TYPE is_null, error;
		//const bool isDate;
		enum_field_types MysqlType; //see mysql/field_types.h
};

namespace mysql_wrp {
	enum MysqlOptions {
		COMPRESSION	= 1,
		ENCRYPTION	= 1 << 1
	};
	
	MYSQL* initMysql(const std::string& mysqlhost, const std::string& mysqluser, const std::string& mysqlpass, const std::string& mysqldb, const unsigned int& options=0);
	MYSQL_STMT* initStmt(MYSQL **mysql, const std::string& query, const long unsigned int& ref_param_count);
	
	void bindParams(MYSQL_STMT **stmt, std::vector<SQL_FIELD> &params_fields);
	void bindResults(MYSQL_STMT **stmt, std::vector<SQL_FIELD> &result_fields);
	double retrieveData(const SQL_FIELD &field, const unsigned int& conversion=SQL_FIELD::NONE);
}

#endif
