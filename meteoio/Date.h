/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef __DATE_H__
#define __DATE_H__

#include <meteoio/IOExceptions.h>

#include <string>
#include <sstream>
#include <utility>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <ctime>

///Using the following namespace for the comparison operator overloading
//using namespace rel_ops;

namespace mio {

/**
 * @class Date
 * @brief  A class to handle timestamps.
 * This class handles conversion between different time display formats (ISO, numeric) as well as different
 * time representation (julian date, modified julian date, etc). It also handles time zones as well as
 * very basic Daylight Saving Time (DST). Since the activation dates of DST are political and not technical,
 * it can not be automatically calculated. Therefore, it has to be provided by the caller: when the dst flag 
 * is set, the dst time shift is automatically applied. When the dst flag ceases to be set, the dst time shift
 * is no longer applied. This is very crude, but please keep in mind that using DST for monitoring data is
 * usually a bad idea...
 * 
 * Internally, the date is stored as true julian date in GMT.
 * The maximal precision is 1 minute (that can be easily brought to 1 seconds if
 * it would appear necessary/useful, with the limitation that leap seconds are currently not handled).
 *
 * Please see Date::FORMATS for supported display formats and http://en.wikipedia.org/wiki/Julian_day for
 * the various date representation definitions. The following data representation are currently supported:
 * - julian date, see Date::getJulianDate
 * - modified julian date, see Date::getModifiedJulianDate
 * - truncated julian date, see Date::getTruncatedJulianDate
 * - Unix date, see Date::getUnixDate
 * - Excel date, see Date::getExcelDate
 * 
 * @author Mathias Bavay
 * @date 2010-04-15
 */

#ifdef _POPC_
class DateDummy {}; //HACK for POPC

class Date : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class Date {
#endif  
	public:
		///Keywords for selecting the date formats
		typedef enum {
			ISO, ///< ISO 8601 extended format combined date: YYYY-MM-DDTHH:mm:SS (fields might be dropped, in the least to the most significant order)
			FULL, ///< ISO 8601 followed by the julian date (in parenthesis)
			NUM, ///< ISO 8601 basic format date: YYYYMMDDHHmmSS (fields might be dropped, in the least to the most significant order)
			DIN ///<DIN5008 format: DD.MM.YYYY HH:MM
		} FORMATS;
		static const int daysLeapYear[];
		static const int daysNonLeapYear[];
		static const double DST_shift;
		static const float MJD_offset;
		static const float Unix_offset;
		static const float Excel_offset;

		Date();
		Date(const double& julian_in, const double& _timezone=undefined, const bool& _dst=false);
		Date(const int& year, const int& month, const int& day, const int& hour, const int& minute, const double& _timezone=undefined, const bool& _dst=false);
		Date(const time_t&, const double& _timezone=undefined, const bool& _dst=false);
		Date(const Date& _date_in);

		void setTimeZone(const double& _timezone, const bool& _dst=false);
		void setDate(const double& julian_in, const double& _timezone=undefined, const bool& _dst=false);
		void setDate(const int& year, const int& month, const int& day, const int& hour, const int& minute, const double& _timezone=undefined, const bool& _dst=false);
		void setDate(const time_t& _time, const double& _timezone=undefined, const bool& _dst=false);
		void setModifiedJulianDate(const double& julian_in, const double& _timezone=undefined, const bool& _dst=false);
		void setUnixDate(const time_t& _time, const double& _timezone=undefined, const bool& _dst=false);
		void setExcelDate(const double excel_in, const double& _timezone=undefined, const bool& _dst=false);

		double getTimeZone() const;
		bool getDST() const;
		double getJulianDate(const bool& gmt=false) const;
		double getModifiedJulianDate(const bool& gmt=false) const;
		double getTruncatedJulianDate(const bool& gmt=false) const;
		time_t getUnixDate(const bool& gmt=false) const;
		double getExcelDate(const bool& gmt=false) const;

		void getDate(double& julian_out, const bool& gmt=false) const;
		void getDate(int& year, int& month, int& day, const bool& gmt=false) const;
		void getDate(int& year, int& month, int& day, int& hour, const bool& gmt=false) const;
		void getDate(int& year, int& month, int& day, int& hour, int& minute, const bool& gmt=false) const;
		int getYear(const bool& gmt=false) const;

		int getJulianDayNumber() const;
		bool isLeapYear() const;

		const std::string toString(FORMATS type, const bool& gmt=false) const;

		friend std::ostream& operator<<(std::ostream& os, const Date& date);

		//Operator Prototypes
		///Can be used to add an interval to an existing Date object.
		///Construct a Date object representing the interval e.g. Date(1.0) for 1 day and add that to another Date object.
		Date& operator+=(const Date&);
		///Can be used to subtract an interval from an existing Date object
		Date& operator-=(const Date&);
		bool operator==(const Date&) const;
		bool operator!=(const Date&) const;
		bool operator<(const Date&) const;
		bool operator<=(const Date&) const;
		bool operator>(const Date&) const;
		bool operator>=(const Date&) const;

		const Date operator+(const Date&) const;
		const Date operator-(const Date&) const;

	private:
		double localToGMT(const double& _julian)const;
		double GMTToLocal(const double& _gmt_julian) const;
		double calculateJulianDate(const int& _year, const int& _month, const int& _day, const int& _hour, const int& _minute) const;
		void calculateValues(const double& julian, int& _year, int& _month, int& _day, int& _hour, int& _minute) const;
		long getJulianDayNumber(const int&, const int&, const int&) const;
		bool isLeapYear(const int&) const;
		void plausibilityCheck(const int& in_year, const int& in_month, const int& in_day, const int& in_hour, const int& in_minute) const;

		double timezone;
		bool dst;
		double gmt_julian;
		int gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute;
		static const double undefined;
};
} //end namespace

#endif
