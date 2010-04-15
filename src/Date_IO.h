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
#ifndef __DATE_IO_H__
#define __DATE_IO_H__

#include <string>
#include <sstream>
#include <utility>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <ctime>
#include "IOExceptions.h"

///Using the following namespace for the comparison operator overloading
//using namespace rel_ops; 

/**
 * @class Date_IO
 * @brief  A relatively basic class designed to represent julian dates and to provide basic conversions into yyyy/mm/dd/hh/mm representations
 *
 * The maximal precision as implemented is 1 minute.
 *
 * @author Thomas Egger
 */
#ifdef _POPC_
class Date_IO : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class Date_IO {
#endif  
	public:
		///Keywords for selecting the date formats
		typedef enum {
			ISO, ///< ISO 8601 extended format combined date: YYYY-MM-DDTHH:mm:SS (fields might be dropped, in the least to the most significant order)
			FULL, ///< ISO 8601 followed by the julian date (in parenthesis)
			NUM ///< ISO 8601 basic format date: YYYYMMDDHHmmSS (fields might be dropped, in the least to the most significant order)
		} FORMATS;
		static const int daysLeapYear[];
		static const int daysNonLeapYear[];
		static const double DST_shift;
		static const float MJD_offset;
		static const float Unix_offset;
		static const float Excel_offset;

		///Note that constructing an Date_IO object without any parameter results 
		///in constructing Date_IO(0.0) which in its current expression results in a date of -4713-01-01T12:00
		Date_IO();
		Date_IO(const double& julian_in, const double& _timezone=0.0, const bool& _dst=false);
		///All values passed will be checked for plausibility 
		Date_IO(const int& year, const int& month, const int& day, const int& hour, const int& minute, const double& _timezone=0.0, const bool& _dst=false);
		Date_IO(const time_t&, const double& _timezone=0.0, const bool& _dst=false);
		Date_IO(const Date_IO& _date_in);

		void setTimeZone(const double& _timezone, const bool& _dst);
		void setDate(const double& julian_in, const double& _timezone=0.0, const bool& _dst=false);
		void setDate(const int& year, const int& month, const int& day, const int& hour, const int& minute, const double& _timezone=0.0, const bool& _dst=false);
		void setDate(const time_t& _time, const double& _timezone=0.0, const bool& _dst=false);
		void setModifiedJulianDate(const double& julian_in, const double& _timezone=0.0, const bool& _dst=false);
		void setUnixDate(const time_t& _time, const double& _timezone=0.0, const bool& _dst=false);
		void setExcelDate(const double julian_in, const double& _timezone=0.0, const bool& _dst=false);

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

		const std::string toString(FORMATS type, const bool& gmt=false) const;

		//Operator Prototypes
		///Can be used to add an interval to an existing Date_IO object.
		///Construct a Date_IO object representing the interval e.g. Date_IO(1.0) for 1 day and add that to another Date_IO object.
		Date_IO& operator+=(const Date_IO&);
		///Can be used to subtract an interval from an existing Date_IO object
		Date_IO& operator-=(const Date_IO&);
		bool operator==(const Date_IO&) const;
		bool operator!=(const Date_IO&) const;
		bool operator<(const Date_IO&) const;
		bool operator<=(const Date_IO&) const;
		bool operator>(const Date_IO&) const;
		bool operator>=(const Date_IO&) const;

		const Date_IO operator+(const Date_IO&) const;
		const Date_IO operator-(const Date_IO&) const;

		friend std::ostream& operator<<(std::ostream& os, const Date_IO& date);

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
};


#endif
