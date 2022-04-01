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

MYSQL* initMysql(const std::string& mysqlhost, const std::string& mysqluser, const std::string& mysqlpass, const std::string& mysqldb, const unsigned int& options)
{
	MYSQL *mysql = mysql_init(nullptr);
	
	//set some options
	unsigned int timeout = 2; // in seconds
	mysql_options(mysql, MYSQL_OPT_CONNECT_TIMEOUT, &timeout);
	if ((COMPRESSION & options) == COMPRESSION)
		mysql_options(mysql, MYSQL_OPT_COMPRESS, 0);
	if ((ENCRYPTION & options) == ENCRYPTION) {
		unsigned int enforce_ssl = SSL_MODE_REQUIRED;
		mysql_options(mysql, MYSQL_OPT_SSL_MODE, &enforce_ssl);
	}
	
	if (!mysql_real_connect(mysql, mysqlhost.c_str(), mysqluser.c_str(), mysqlpass.c_str(), mysqldb.c_str(), 0, NULL, 0))
		throw AccessException("Could not initiate connection to Mysql server "+mysqlhost+": "+std::string(mysql_error(mysql)), AT);

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
		if (param_count!=ref_param_count) throw InvalidArgumentException("Wrong number of parameters in mysql statement", AT);
	}
	
	return stmt;
}

void bindParams(MYSQL_STMT **stmt, std::vector<fType> &params_fields)
{
	const size_t params_count = params_fields.size();
	MYSQL_BIND *stmtParams = (MYSQL_BIND*)calloc(params_count, sizeof(MYSQL_BIND));
	if (stmtParams==nullptr) throw IOException("Could not allocate memory for parameter binding to Mysql query", AT);
	
	for(size_t ii=0; ii<params_count; ++ii) {
		stmtParams[ii].buffer_type = params_fields[ii].MysqlType;
		stmtParams[ii].is_null = nullptr;
		
		if (params_fields[ii].MysqlType==MYSQL_TYPE_STRING) {
			stmtParams[ii].buffer = (char *)params_fields[ii].str;
			stmtParams[ii].buffer_length = STRING_SIZE;
			stmtParams[ii].length = &params_fields[ii].str_len;
		} else if(params_fields[ii].MysqlType==MYSQL_TYPE_DOUBLE) {
			stmtParams[ii].buffer = (char *)&params_fields[ii].val;
		} else if(params_fields[ii].MysqlType==MYSQL_TYPE_DATETIME) {
			stmtParams[ii].buffer = (char *)&params_fields[ii].dt;
		}
	}
	
	if (mysql_stmt_bind_param(*stmt, stmtParams)) {
		free( stmtParams );
		throw IOException("Error binding parameters", AT);
	}
	free( stmtParams );
}

void bindResults(MYSQL_STMT **stmt, std::vector<fType> &result_fields)
{
	MYSQL_RES *prepare_meta_result = mysql_stmt_result_metadata(*stmt);
	if (!prepare_meta_result) throw IOException("Error executing meta statement", AT);
	const size_t column_count = static_cast<size_t>( mysql_num_fields(prepare_meta_result) );
	if (column_count!=result_fields.size()) throw InvalidArgumentException("Wrong number of columns returned", AT);
	mysql_free_result(prepare_meta_result);
	
	MYSQL_BIND *result = (MYSQL_BIND*)calloc(column_count, sizeof(MYSQL_BIND));
	if (result==nullptr) throw IOException("Could not allocate memory for results binding to Mysql query", AT);
	
	for(size_t ii=0; ii<column_count; ++ii) {
		result[ii].buffer_type = result_fields[ii].MysqlType;
		if (result_fields[ii].MysqlType==MYSQL_TYPE_STRING) {
			result[ii].buffer = (char *)result_fields[ii].str;
			result[ii].buffer_length = STRING_SIZE;
		} else if(result_fields[ii].MysqlType==MYSQL_TYPE_DOUBLE) {
			result[ii].buffer_type = MYSQL_TYPE_DOUBLE;
			result[ii].buffer = (char *)&result_fields[ii].val;
		} else if(result_fields[ii].MysqlType==MYSQL_TYPE_DATETIME) {
			result[ii].buffer_type = MYSQL_TYPE_DATETIME;
			result[ii].buffer = (char *)&result_fields[ii].dt;
		}
		
		result[ii].is_null = &result_fields[ii].is_null;
		result[ii].length = &result_fields[ii].buffer_len;
		result[ii].error = &result_fields[ii].error;
	}
	
	if (mysql_stmt_bind_result(*stmt, result)) throw IOException("Error binding results", AT);
	if (mysql_stmt_store_result(*stmt)) throw IOException("mysql_stmt_store_result failed", AT);
	free( result );
}

double retrieveData(const fType &field, const unsigned int& conversion)
{
	const double val = field.val;
	if (field.is_null==1) return IOUtils::nodata;
	
	if (conversion==C_TO_K) return IOUtils::C_TO_K( val );
	if (conversion==NORMALIZE_PC || conversion==CM_TO_M) return val / 100.;
	
	return val;
}

}
