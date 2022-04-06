#SPDX-License-Identifier: LGPL-3.0-or-later
# - Find mysqlclient
# Find the native MySQL includes and library
#
#  MYSQL_INCLUDE_DIR - where to find mysql.h, etc.
#  MYSQL_LIBRARY   - List of libraries when using MySQL.
#  MYSQL_FOUND       - True if MySQL found.

INCLUDE(LibFindMacros)

IF(WIN32)
	SET(SEARCH_PATH
		ENV LIB
		"$ENV{PROGRAMFILES}/MySQL/*/lib"
		"$ENV{PROGRAMFILES\(x86\)}/MySQL/*/lib"
		"$ENV{SYSTEMDRIVE}/MySQL/*/lib"
	)
	FIND_LIBRARY(MySQL_LIBRARY
		NAMES libmysql
		HINTS ${SEARCH_PATH}
		DOC "Location of the libmysql dynamic library (dll)"
		)
ELSE(WIN32)
	IF(APPLE)
		FIND_LIBRARY(MySQL_LIBRARY
		NAMES libmysqlclient.dylib
		HINTS
			ENV LD_LIBRARY_PATH
			ENV DYLD_FALLBACK_LIBRARY_PATH
			"~/usr/lib"
			"/usr/local/lib"
			"/usr/lib"
			"/opt/lib"
		DOC "Location of the libmysql dynamic library"
		)
	ELSE(APPLE)
		FIND_LIBRARY(MySQL_LIBRARY
		NAMES libmysqlclient.so
		HINTS
			ENV LD_LIBRARY_PATH
			"~/usr/lib"
			"/usr/local/lib"
			"/usr/local/lib64"
			"/usr/lib"
			"/usr/lib64"
			"/usr/lib/x86_64-linux-gnu"
			"/opt/lib"
			"/opt/lib64"
		DOC "Location of the libmysql dynamic library"
		)
	ENDIF(APPLE)
ENDIF(WIN32)

IF(MySQL_LIBRARY)
	#build MySQL_ROOT so we can provide a hint for searching for the matching header file
	GET_FILENAME_COMPONENT(MySQL_ROOT ${MySQL_LIBRARY} DIRECTORY) #get PATH
	GET_FILENAME_COMPONENT(MySQL_ROOT ${MySQL_ROOT} DIRECTORY) #go up one level

	# locate main header file
	FIND_PATH(MySQL_INCLUDE_DIR
	NAMES mysql.h
	HINTS
		"${MySQL_ROOT}/include"
		"~/usr/include/mysql"
		"/usr/local/include/mysql"
		"/usr/include/mysql"
		"/opt/include/mysql"
	DOC "Location of the MySQL headers, like /usr/include/mysql/include"
	)
ENDIF(MySQL_LIBRARY)

if( MySQL_INCLUDE_DIR AND EXISTS "${MySQL_INCLUDE_DIR}/mysql_version.h" )
	file( STRINGS "${MySQL_INCLUDE_DIR}/mysql_version.h"
		MYSQL_VERSION_H REGEX "^#define[ \t]+MYSQL_SERVER_VERSION[ \t]+\"[^\"]+\".*$" )
	string( REGEX REPLACE
		"^.*MYSQL_SERVER_VERSION[ \t]+\"([^\"]+)\".*$" "\\1" MYSQL_VERSION_STRING
		"${MYSQL_VERSION_H}" )
endif()

# handle the QUIETLY and REQUIRED arguments and set MYSQL_FOUND to TRUE if
# all listed variables are TRUE
include( FindPackageHandleStandardArgs )
find_package_handle_standard_args(MySQL DEFAULT_MSG
                                  MySQL_LIBRARY MySQL_INCLUDE_DIR)

set( MYSQL_INCLUDE_DIR ${MySQL_INCLUDE_DIR} )
set( MYSQL_LIBRARY ${MySQL_LIBRARY} )

mark_as_advanced( MySQL_INCLUDE_DIR MySQL_LIBRARY )
mark_as_advanced( MYSQL_INCLUDE_DIR MYSQL_LIBRARY )
