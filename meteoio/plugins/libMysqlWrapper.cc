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
#include <meteoio/plugins/libMysqlWrapper.h>
#include <meteoio/IOExceptions.h>

#include <stdio.h>

using namespace std;
using namespace mio;

namespace mysql_wrp {

MYSQL* initMysql(const std::string& mysqlhost, const std::string& mysqluser, const std::string& mysqlpass, const std::string& mysqldb)
{
	MYSQL *mysql = mysql_init(nullptr);
	//mysql_options(mysql, MYSQL_OPT_RECONNECT, &reconnect);
	if (!mysql_real_connect(mysql, mysqlhost.c_str(), mysqluser.c_str(), mysqlpass.c_str(), mysqldb.c_str(), 0, NULL, 0))
		throw IOException("Could not initiate connection to Mysql server "+mysqlhost+": "+std::string(mysql_error(mysql)), AT);

	return mysql;
}

MYSQL_STMT* initStmt(MYSQL **mysql, const std::string& query, const long unsigned int& ref_param_count)
{
	MYSQL_STMT* stmt = mysql_stmt_init(*mysql);
	if (!stmt) throw IOException("Could not allocate memory for mysql statement", AT);

	if (mysql_stmt_prepare(stmt, query.c_str(), query.size())) {
		throw IOException("Error preparing mysql statement", AT);
	} else {
		const long unsigned int param_count = mysql_stmt_param_count(stmt);
		if (param_count!=ref_param_count) throw IOException("Wrong number of parameters in mysql statement", AT);
	}
	
	return stmt;
}

void bindParams(MYSQL_STMT **stmt, std::vector<fType> &params_fields)
{
	const size_t params_count = params_fields.size();
	MYSQL_BIND stmtParams[params_count];
	memset(stmtParams, 0, sizeof(stmtParams));
	
	for(size_t ii=0; ii<params_count; ++ii) {
		stmtParams[ii].buffer_type = params_fields[ii].MysqlType;
		if (params_fields[ii].MysqlType==MYSQL_TYPE_STRING) {
			stmtParams[ii].buffer = (char *)params_fields[ii].str;
			stmtParams[ii].buffer_length = STRING_SIZE;
			stmtParams[ii].length = &params_fields[ii].str_len;
		} else if(params_fields[ii].MysqlType==MYSQL_TYPE_DOUBLE) {
			stmtParams[ii].buffer_type= MYSQL_TYPE_DOUBLE;
			stmtParams[ii].buffer= (char *)&params_fields[ii].val;
		}
		
		stmtParams[ii].is_null = 0;
	}
	
	if (mysql_stmt_bind_param(*stmt, stmtParams)) {
		throw IOException("Error binding parameters", AT);
	}
}

void bindResults(MYSQL_STMT **stmt, std::vector<fType> &result_fields)
{
	MYSQL_RES *prepare_meta_result = mysql_stmt_result_metadata(*stmt);
	if (!prepare_meta_result) throw IOException("Error executing meta statement", AT);
	const size_t column_count = static_cast<size_t>( mysql_num_fields(prepare_meta_result) );
	if (column_count!=result_fields.size()) throw IOException("Wrong number of columns returned", AT);
	mysql_free_result(prepare_meta_result);
	
	MYSQL_BIND result[column_count];
	memset(result, 0, sizeof(result));
	unsigned long length[column_count];
#if LIBMYSQL_VERSION_ID > 80001 
	bool is_null[column_count];
	bool error[column_count];
#else
	my_bool is_null[column_count];
	my_bool error[column_count];
#endif
	
	for(size_t ii=0; ii<column_count; ++ii) {
		result[ii].buffer_type = result_fields[ii].MysqlType;
		if (result_fields[ii].MysqlType==MYSQL_TYPE_STRING) {
			result[ii].buffer = (char *)result_fields[ii].str;
			result[ii].buffer_length = STRING_SIZE;
		} else if(result_fields[ii].MysqlType==MYSQL_TYPE_DOUBLE) {
			result[ii].buffer_type= MYSQL_TYPE_DOUBLE;
			result[ii].buffer= (char *)&result_fields[ii].val;
		}
		
		result[ii].is_null = &is_null[ii];
		result[ii].length = &length[ii];
		result[ii].error = &error[ii];
	}
	
	if (mysql_stmt_bind_result(*stmt, result)) throw IOException("Error binding results", AT);
	if (mysql_stmt_store_result(*stmt)) throw IOException("mysql_stmt_store_result failed", AT);
}

}
