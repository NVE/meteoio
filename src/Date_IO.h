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

using namespace std; 
///Using the following namespace for the comparison operator overloading
using namespace rel_ops; 

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
		typedef enum {
			ISO, ///< ISO 8601 extended format combined date: YYYY-MM-DDTHH:mm:SS (fields might be dropped, in the least to the most significant order)
			NUM ///< ISO 8601 basic format date: YYYYMMDDHHmmSS (fields might be dropped, in the least to the most significant order)
		} FORMATS;
		static const int daysLeapYear[];
		static const int daysNonLeapYear[];
		static const long offset;

		///Note that constructing an Date_IO object without any parameter results 
		///in constructing Date_IO(0.0) which in its current expression results in a date of 1900/1/1 00:00:00
		//TODO: use modified julian date instead (ie: 0=1858-11-17T00:00:00)
		Date_IO(const double& julian_in=0.0);
		///All values passed will be checked for plausibility 
		Date_IO(const int& year, const int& month, const int& day=1, const int& hour=0, const int& minute=0);
		Date_IO(const time_t&);
		Date_IO(const Date_IO& _date_in);

		void setDate(const double& julian_in);
		void setDate(const int& year, const int& month, const int& day=1, const int& hour=0, const int& minute=0);
		void setDate(const time_t& _time);

		double getJulian() const;

		void getDate(double& julian_out) const;
		void getDate(int& year, int& month, int& day) const;
		void getDate(int& year, int& month, int& day, int& hour) const;
		void getDate(int& year, int& month, int& day, int& hour, int& minute) const;
		int getYear() const;

		//Dealing with epoch time (seconds since 00:00:00 on January 1 (UTC))
		time_t getEpochTime() const;

		///Since at SLF julian dates are always treated with an offset, this function provides a way to deal with real julian dates
		void setRealJulianDate(const double& julian_in);
		///Since at SLF julian dates are always treated with an offset, this function provides a way to deal with real julian dates
		void getRealJulianDate(double& julian_out) const;

		///The toString representation outputs the date in the form of "[double]julian_date yyyy/mm/dd hh:mm"
		const std::string toString(void) const;
		friend ostream& operator<<(ostream& os, const Date_IO& date);

		const std::string toString(FORMATS type) const;

		//Operator Prototypes
		///Can be used to add an interval to an existing Date_IO object. 
		///Construct a Date_IO object representing the interval e.g. Date_IO(1.0) for 1 day and add that to another Date_IO object.
		Date_IO& operator+=(const Date_IO&);
		///Can be used to subtract an interval from an existing Date_IO object
		Date_IO& operator-=(const Date_IO&);
		bool operator==(const Date_IO&) const;
		bool operator<(const Date_IO&) const;

		const Date_IO operator+(const Date_IO&) const;
		const Date_IO operator-(const Date_IO&) const;

	private:
		void calculateJulianDate(void);
		void calculateValues(void);
		long getJulianDay(const int&, const int&, const int&) const;
		bool isLeapYear(const int&) const;
		void plausibilityCheck(const int& in_year, const int& in_month, const int& in_day, const int& in_hour, const int& in_minute) const;

		double julian, realjulian;
		int year, month, day, hour, minute;
};


#endif
