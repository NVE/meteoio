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
#include <meteoio/Date.h>

//const long Date::offset = 2415021;	///snowpack offset

using namespace std;

namespace mio {

const int Date::daysLeapYear[12] = {31,29,31,30,31,30,31,31,30,31,30,31};
const int Date::daysNonLeapYear[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
const double Date::DST_shift = 1.0; //in hours
const float Date::MJD_offset = 2400000.5; ///<offset between julian date and modified julian date
const float Date::Unix_offset = 2440587.5; ///<offset between julian date and Unix Epoch time
const float Date::Excel_offset = 2415018.5;  ///<offset between julian date and Excel dates (note that excel invented some days...)

// CONSTUCTORS
/**
* @brief Default constructor: timezone is set to GMT without DST, julian date is set to 0 (meaning -4713-01-01T12:00)
*/
Date::Date() {
	undef = true;
	timezone = 0;
	dst = false;
	setDate(0., 0., false);
}

/**
* @brief Julian date constructor.
* @param julian_in julian date to set
* @param in_timezone timezone as an offset to GMT (in hours, optional)
* @param in_dst is it DST? (default: no)
*/
Date::Date(const double& julian_in, const double& in_timezone, const bool& in_dst) {
	timezone = 0;
	dst = false;
	setDate(julian_in, in_timezone, in_dst);
}

/**
* @brief Unix date constructor.
* @param in_time unix time (ie: as number of seconds since Unix Epoch)
* @param in_timezone timezone as an offset to GMT (in hours, optional)
* @param in_dst is it DST? (default: no)
*/
Date::Date(const time_t& in_time, const double& in_timezone, const bool& in_dst) {
	timezone = 0;
	dst = false;
	setDate(in_time, in_timezone, in_dst);
}
//HACK: is it needed? Why paroc_base instead of POPC??
/**
* @brief Copy constructor.
* @param in_date Date object to copy
*/
#ifdef _POPC_
Date::Date(const Date& in_date) : paroc_base()
#else
Date::Date(const Date& in_date)
#endif
{
	setDate(in_date.getJulianDate(), in_date.getTimeZone(), in_date.getDST());
}

/**
* @brief Date constructor by elements.
* All values are checked for plausibility.
* @param in_year in 4 digits
* @param in_month please keep in mind that first month of the year is 1 (ie: not 0!)
* @param in_day please keep in mind that first day of the month is 1 (ie: not 0!)
* @param in_hour
* @param in_minute
* @param in_timezone timezone as an offset to GMT (in hours, optional)
* @param in_dst is it DST? (default: no)
*/
Date::Date(const int& in_year, const int& in_month, const int& in_day, const int& in_hour, const int& in_minute, const double& in_timezone, const bool& in_dst)
{
	timezone = 0;
	dst = false;
	setDate(in_year, in_month, in_day, in_hour, in_minute, in_timezone, in_dst);
}

// SETTERS
void Date::setUndef(const bool& flag) {
	undef = flag;
}
/**
* @brief Set internal gmt time from system time
*/
void Date::setFromSys(const double& in_timezone) {
	setDate( time(NULL), in_timezone ); //Unix time_t setter
}

/**
* @brief Set timezone and Daylight Saving Time flag.
* @param in_timezone timezone as an offset to GMT (in hours)
* @param in_dst is it DST?
*/
void Date::setTimeZone(const double& in_timezone, const bool& in_dst) {
//please keep in mind that timezone might be fractional (ie: 15 minutes, etc)
	if(abs(in_timezone) > 12) {
		throw InvalidArgumentException("[E] Time zone can NOT be greater than +/-12!!", AT);
	}

	timezone = in_timezone;
	dst = in_dst;
}

/**
* @brief Set date by elements.
* All values are checked for plausibility.
* @param _year in 4 digits
* @param _month please keep in mind that first month of the year is 1 (ie: not 0!)
* @param _day please keep in mind that first day of the month is 1 (ie: not 0!)
* @param _hour
* @param _minute
* @param _timezone timezone as an offset to GMT (in hours, optional)
* @param _dst is it DST? (default: no)
*/
void Date::setDate(const int& _year, const int& _month, const int& _day, const int& _hour, const int& _minute, const double& _timezone, const bool& _dst)
{
	plausibilityCheck(_year, _month, _day, _hour, _minute);
	undef = false;
	setTimeZone(_timezone, _dst);

	if(timezone==0 && dst==false) {
		//data is GMT and no DST
		//setting values and computing GMT julian date
		gmt_year = _year;
		gmt_month = _month;
		gmt_day = _day;
		gmt_hour = _hour;
		gmt_minute = _minute;
		gmt_julian = calculateJulianDate(gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute);
	} else {
		//computing local julian date
		const double local_julian = calculateJulianDate(_year, _month, _day, _hour, _minute);
		//converting local julian date to GMT julian date
		gmt_julian = localToGMT(local_julian);
		//updating values to GMT
		calculateValues(gmt_julian, gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute);
	}

}

/**
* @brief Set date from a julian date (JD).
* @param julian_in julian date to set
* @param in_timezone timezone as an offset to GMT (in hours, optional)
* @param in_dst is it DST? (default: no)
*/
void Date::setDate(const double& julian_in, const double& in_timezone, const bool& in_dst) {
	setTimeZone(in_timezone, in_dst);
	gmt_julian = localToGMT(julian_in);
	undef = false;
	calculateValues(gmt_julian, gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute);
}

/**
* @brief Set date from a Unix date.
* @param _time unix time (ie: as number of seconds since Unix Epoch)
* @param _timezone timezone as an offset to GMT (in hours, optional)
* @param _dst is it DST? (default: no)
*/
void Date::setDate(const time_t& _time, const double& _timezone, const bool& _dst) {
	setUnixDate(_time, _timezone, _dst);
}

/**
* @brief Set date from a modified julian date (MJD).
* @param julian_in julian date to set
* @param _timezone timezone as an offset to GMT (in hours, optional)
* @param _dst is it DST? (default: no)
*/
void Date::setModifiedJulianDate(const double& julian_in, const double& _timezone, const bool& _dst) {
	const double _julian = julian_in + MJD_offset;
	setDate(_julian, _timezone, _dst);
}

/**
* @brief Set date from a Unix date.
* @param in_time unix time (ie: as number of seconds since Unix Epoch)
* @param in_timezone timezone as an offset to GMT (in hours, optional)
* @param in_dst is it DST? (default: no)
*/
void Date::setUnixDate(const time_t& in_time, const double& in_timezone, const bool& in_dst) {
	const double in_julian = (double)(in_time)/(24.*60.*60.) + Unix_offset;
	setDate(in_julian, in_timezone, in_dst);
}

/**
* @brief Set date from an Excel date.
* @param excel_in Excel date to set
* @param _timezone timezone as an offset to GMT (in hours, optional)
* @param _dst is it DST? (default: no)
*/
void Date::setExcelDate(const double excel_in, const double& _timezone, const bool& _dst) {
	//TODO: handle date < 1900-01-00 and date before 1900-03-01
	const double _julian = excel_in + Excel_offset;
	setDate(_julian, _timezone, _dst);
}

// GETTERS
bool Date::isUndef() const {
	return undef;
}

/**
* @brief Returns timezone.
* @return timezone as an offset to GMT
*/
double Date::getTimeZone() const {
	if(undef==true)
		throw UnknownValueException("Date object is undefined!", AT);

	return timezone;
}

/**
* @brief Returns Daylight Saving Time flag.
* @return dst enabled?
*/
bool Date::getDST() const {
	if(undef==true)
		throw UnknownValueException("Date object is undefined!", AT);

	return dst;
}

/**
* @brief Return julian date (JD).
* The julian date is defined as the fractional number of days since -4713-01-01T12:00 UTC.
* @param gmt convert returned value to GMT? (default: false)
* @return julian date in the current timezone / in GMT depending on the gmt parameter
*/
double Date::getJulianDate(const bool& gmt) const {
	if(undef==true)
		throw UnknownValueException("Date object is undefined!", AT);

	if(gmt) {
		return gmt_julian;
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		return local_julian;
	}
}

/**
* @brief Return modified julian date (MJD).
* The modified julian date is defined as the fractional number of days since 1858-11-17T00:00 UTC
* (definition by the Smithsonian Astrophysical Observatory, MA).
* @param gmt convert returned value to GMT? (default: false)
* @return modified julian date in the current timezone / in GMT depending on the gmt parameter
*/
double Date::getModifiedJulianDate(const bool& gmt) const {
	if(undef==true)
		throw UnknownValueException("Date object is undefined!", AT);

	if(gmt) {
		return (gmt_julian - MJD_offset);
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		return (local_julian - MJD_offset);
	}
}

/**
* @brief Return truncated julian date (TJD).
* The truncated julian date is defined as the julian day shifted to start at 00:00 and modulo 10000 days.
* The last origin (ie: 0) was 1995-10-10T00:00 
* (definition by National Institute of Standards and Technology).
* @param gmt convert returned value to GMT? (default: false)
* @return truncated julian date in the current timezone / in GMT depending on the gmt parameter
*/
double Date::getTruncatedJulianDate(const bool& gmt) const {
	if(undef==true)
		throw UnknownValueException("Date object is undefined!", AT);

	if(gmt) {
		return (fmod( (gmt_julian - 0.5), 10000. ));
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		return (fmod( (local_julian - 0.5), 10000. ));
	}
}

/**
* @brief Return Unix time (or POSIX time).
* The Unix time is defined as the number of seconds since 1970-01-01T00:00 UTC (Unix Epoch).
* (defined as IEEE P1003.1 POSIX. See http://www.mail-archive.com/leapsecs@rom.usno.navy.mil/msg00109.html
* for some technical, historical and funny insight into the standardization process)
* @param gmt convert returned value to GMT? (default: false)
* @return Unix time in the current timezone / in GMT depending on the gmt parameter
*/
time_t Date::getUnixDate(const bool& gmt) const { //HACK: should Unix date always be GMT?
	if(undef==true)
		throw UnknownValueException("Date object is undefined!", AT);

	if (gmt_julian < Unix_offset)
			throw IOException("Dates before 1970 cannot be displayed in Unix epoch time", AT);

	if(gmt) {
		return ( (time_t)floor( (gmt_julian - Unix_offset) * (24*60*60) ));
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		return ( (time_t)floor( (local_julian - Unix_offset) * (24*60*60) ));
	}
}

/**
* @brief Return Excel date.
* The (sick) Excel date is defined as the number of days since 1900-01-00T00:00 (no, this is NOT a typo).
* Moreover, it (wrongly) considers that 1900 was a leap year (in order to remain compatible with an old Lotus123 bug).
* This practically means that for dates after 1900-03-01, an Excel date really represents the number of days since 1900-01-01T00:00 PLUS 2.
* @param gmt convert returned value to GMT? (default: false)
* @return Excel date in the current timezone / in GMT depending on the gmt parameter
*/
double Date::getExcelDate(const bool& gmt) const {
	if(undef==true)
		throw UnknownValueException("Date object is undefined!", AT);

	if (gmt_julian < Excel_offset)
		throw IOException("Dates before 1900 cannot be converted to Excel date", AT);

	if(gmt) {
		return ( gmt_julian - Excel_offset);
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		return ( local_julian - Excel_offset);
	}
}

/**
* @brief Retrieve julian date.
* This method is a candidate for deletion: it should now be obsolete.
* @param julian_out julian date (in local time zone or GMT depending on the gmt flag)
* @param gmt convert returned value to GMT? (default: false)
*/
void Date::getDate(double& julian_out, const bool& gmt) const {
	if(undef==true)
		throw UnknownValueException("Date object is undefined!", AT);

	if(gmt) {
		julian_out = gmt_julian;
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		julian_out = local_julian;
	}
}

/**
* @brief Return year.
* @param gmt convert returned value to GMT? (default: false)
* @return year
*/
int Date::getYear(const bool& gmt) const {
	if(undef==true)
		throw UnknownValueException("Date object is undefined!", AT);

	if(gmt) {
		return gmt_year;
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		int local_year, local_month, local_day, local_hour, local_minute;
		calculateValues(local_julian, local_year, local_month, local_day, local_hour, local_minute);
		return local_year;
	}
}

/**
* @brief Return year, month, day.
* @param year_out
* @param month_out
* @param day_out
* @param gmt convert returned value to GMT? (default: false)
*/
void Date::getDate(int& year_out, int& month_out, int& day_out, const bool& gmt) const {
	if(undef==true)
		throw UnknownValueException("Date object is undefined!", AT);

	if(gmt) {
		year_out = gmt_year;
		month_out = gmt_month;
		day_out = gmt_day;
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		int local_hour, local_minute;
		calculateValues(local_julian, year_out, month_out, day_out, local_hour, local_minute);
	}
}

/**
* @brief Return year, month, day.
* @param year_out
* @param month_out
* @param day_out
* @param hour_out
* @param gmt convert returned value to GMT? (default: false)
*/
void Date::getDate(int& year_out, int& month_out, int& day_out, int& hour_out, const bool& gmt) const {
	if(undef==true)
		throw UnknownValueException("Date object is undefined!", AT);

	if(gmt) {
		year_out = gmt_year;
		month_out = gmt_month;
		day_out = gmt_day;
		hour_out = gmt_hour;
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		int local_minute;
		calculateValues(local_julian, year_out, month_out, day_out, hour_out, local_minute);
	}
}

/**
* @brief Return year, month, day.
* @param year_out
* @param month_out
* @param day_out
* @param hour_out
* @param minute_out
* @param gmt convert returned value to GMT? (default: false)
*/
void Date::getDate(int& year_out, int& month_out, int& day_out, int& hour_out, int& minute_out, const bool& gmt) const {
	if(undef==true)
		throw UnknownValueException("Date object is undefined!", AT);

	if(gmt) {
		year_out = gmt_year;
		month_out = gmt_month;
		day_out = gmt_day;
		hour_out = gmt_hour;
		minute_out = gmt_minute;
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		calculateValues(local_julian, year_out, month_out, day_out, hour_out, minute_out);
	}
}

/**
* @brief Return year, month, day.
* Return the day of the year index for the current Date object
* @return julian day number
*/
int Date::getJulianDayNumber() const {
	if(undef==true)
		throw UnknownValueException("Date object is undefined!", AT);

	//this is quite inefficient... we might want to deal with leap years with their rule + days arrays instead
	int local_year, local_month, local_day, local_hour, local_minute;
	getDate(local_year, local_month, local_day, local_hour, local_minute);

	return (getJulianDayNumber(local_year, local_month, local_day));
}

/**
* @brief Return true if the current year is a leap year
* @return true if the current year is a leap year
*/
bool Date::isLeapYear() const {
	if(undef==true)
		throw UnknownValueException("Date object is undefined!", AT);

	//this is quite inefficient... we might want to deal with leap years with their rule instead
	int local_year, local_month, local_day, local_hour, local_minute;
	getDate(local_year, local_month, local_day, local_hour, local_minute);

	return (isLeapYear(local_year));
	/*
	if( ((local_year%4 == 0) && (local_year%100 != 0)) || (local_year%400 == 0) ) {
	return true;
} else {
	return false;
}
	*/
}

/**
 * @brief Rounds julian date
 */
double Date::rndJulianDate() const {
	return floor(getJulianDate()*24.*2.) / (24.*2.);
}

// OPERATORS //HACK this will have to handle Durations
//HACK: handle undef dates!!
Date& Date::operator+=(const Date& indate) {
	gmt_julian += indate.gmt_julian;
	calculateValues(gmt_julian, gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute);
	return *this;
}

Date& Date::operator-=(const Date& indate) {
	gmt_julian -= indate.gmt_julian;
	calculateValues(gmt_julian, gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute);
	return *this;
}

Date& Date::operator+=(const double& indate) {
	gmt_julian += indate;
	calculateValues(gmt_julian, gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute);
	return *this;
}

Date& Date::operator-=(const double& indate) {
	gmt_julian -= indate;
	calculateValues(gmt_julian, gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute);
	return *this;
}

Date& Date::operator*=(const double& value) {
	gmt_julian *= value;
	calculateValues(gmt_julian, gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute);
	return *this;
}

Date& Date::operator/=(const double& value) {
	gmt_julian /= value;
	calculateValues(gmt_julian, gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute);
	return *this;
}

bool Date::operator==(const Date& indate) const { 
	const double epsilon=1./(24.*3600.); //that is, 1 second in units of days
	return( fabs(indate.gmt_julian - gmt_julian) < epsilon );
}

bool Date::operator!=(const Date& indate) const {
	return !(*this==indate);
}

bool Date::operator<(const Date& indate) const {
	if (*this == indate) {
		return false;
	}

	return (gmt_julian < indate.gmt_julian);
}

bool Date::operator<=(const Date& indate) const {
	if (*this == indate) {
		return true;
	}
	return (gmt_julian <= indate.gmt_julian);
}

bool Date::operator>(const Date& indate) const {
	if (*this == indate) {
		return false;
	}

	return (gmt_julian > indate.gmt_julian);
}

bool Date::operator>=(const Date& indate) const {
	if (*this == indate) {
		return true;
	}
	return (gmt_julian >= indate.gmt_julian);
}

const Date Date::operator+(const Date& indate) const {
	Date tmp(gmt_julian + indate.gmt_julian, timezone);
	return tmp;
}

const Date Date::operator-(const Date& indate) const {
	Date tmp(gmt_julian - indate.gmt_julian, timezone);
	return tmp;
}

const Date Date::operator+(const double& indate) const {
	Date tmp(gmt_julian + indate, timezone);
	return tmp;
}

const Date Date::operator-(const double& indate) const {
	Date tmp(gmt_julian - indate, timezone);
	return tmp;
}

const Date Date::operator*(const double& value) const {
	Date tmp(gmt_julian * value, timezone);
	return tmp;
}

const Date Date::operator/(const double& value) const {
	Date tmp(gmt_julian / value, timezone);
	return tmp;
}

std::ostream& operator<<(std::ostream &os, const Date &date) {
	os << "<date>\n";
	if(date.undef==true)
		os << "Date is undefined\n";
	else {
		os << date.toString(Date::ISO) << "\n";
		os << "TZ=GMT" << showpos << date.timezone << noshowpos << "\n";
		os << "DST=" << date.dst << "\n";
		os << "julian:\t\t\t" << setprecision(10) << date.getJulianDate() << "\t(GMT=" << date.getJulianDate(true) << ")\n";
		os << "ModifiedJulian:\t\t" << date.getModifiedJulianDate() << "\n";
		os << "TruncatedJulian:\t" << date.getTruncatedJulianDate() << "\n";
		try {
			os << "Unix:\t\t\t" << date.getUnixDate() << "\n";
		} catch (...) {}
		try {
			os << "Excel:\t\t\t" << date.getExcelDate() << "\n";
		} catch (...) {}
	}
	os << "</date>\n";
	return os;
}

/**
* @brief Return a nicely formated string.
* @param type select the formating to apply (see the definition of Date::FORMATS)
* @param gmt convert returned value to GMT? (default: false)
* @return formatted time in a string
*/
const string Date::toString(FORMATS type, const bool& gmt) const
{//the date are displayed in LOCAL timezone (more user friendly)
	int year_out, month_out, day_out, hour_out, minute_out;
	double julian_out;

	if(undef==true)
		throw UnknownValueException("Date object is undefined!", AT);

	if(gmt) {
		julian_out = gmt_julian;
		year_out = gmt_year;
		month_out = gmt_month;
		day_out = gmt_day;
		hour_out = gmt_hour;
		minute_out = gmt_minute;
	} else {
		julian_out = GMTToLocal(gmt_julian);
		calculateValues(julian_out, year_out, month_out, day_out, hour_out, minute_out);
	}

	stringstream tmpstr;
	if(type==ISO) {
			tmpstr 
			<< setw(4) << setfill('0') << year_out << "-"
			<< setw(2) << setfill('0') << month_out << "-"
			<< setw(2) << setfill('0') << day_out << "T"
			<< setw(2) << setfill('0') << hour_out << ":"
			<< setw(2) << setfill('0') << minute_out;
	} else if(type==NUM) {
			tmpstr 
			<< setw(4) << setfill('0') << year_out
			<< setw(2) << setfill('0') << month_out
			<< setw(2) << setfill('0') << day_out
			<< setw(2) << setfill('0') << hour_out
			<< setw(2) << setfill('0') << minute_out ;
	} else if(type==FULL) {
			tmpstr 
			<< setw(4) << setfill('0') << year_out << "-"
			<< setw(2) << setfill('0') << month_out << "-"
			<< setw(2) << setfill('0') << day_out << "T"
			<< setw(2) << setfill('0') << hour_out << ":"
			<< setw(2) << setfill('0') << minute_out << " ("
			<< setprecision(10) << julian_out << ") GMT"
			<< setw(2) << setfill('0') << showpos << timezone << noshowpos;
	} else if(type==DIN) {
			tmpstr 
			<< setw(2) << setfill('0') << day_out << "."
			<< setw(2) << setfill('0') << month_out << "."
			<< setw(4) << setfill('0') << year_out << " "
			<< setw(2) << setfill('0') << hour_out << ":"
			<< setw(2) << setfill('0') << minute_out;
	} else {
		throw InvalidArgumentException("Wrong date conversion format requested", AT);
	}

	return tmpstr.str();
}

// PRIVATE METHODS
double Date::calculateJulianDate(const int& _year, const int& _month, const int& _day, const int& _hour, const int& _minute) const
{
	const long julday = getJulianDayNumber(_year, _month, _day);
	const double frac = (_hour-12.)/24. + _minute/(24.*60.); //the julian date reference is at 12:00

	return (((double)julday) + frac);
}

void Date::calculateValues(const double& _julian, int& _year, int& _month, int& _day, int& _hour, int& _minute) const
{ //given a julian day, calculate the year, month, day, hours and minutes
 //see Fliegel, H. F. and van Flandern, T. C. 1968. Letters to the editor: a machine algorithm for processing calendar dates. Commun. ACM 11, 10 (Oct. 1968), 657. DOI= http://doi.acm.org/10.1145/364096.364097 
	long t1;
	const long julday = (long) floor(_julian+0.5);

	t1 = julday + 68569L;
	const long t2 = 4L * t1 / 146097L;
	t1 = t1 - ( 146097L * t2 + 3L ) / 4L;
	const long yr = 4000L * ( t1 + 1L ) / 1461001L;
	t1 = t1 - 1461L * yr / 4L + 31L;
	const long mo = 80L * t1 / 2447L;

	_day = (int) ( t1 - 2447L * mo / 80L );
	t1 = mo / 11L;
	_month = (int) ( mo + 2L - 12L * t1 );
	_year = (int) ( 100L * ( t2 - 49L ) + yr + t1 );

	// Correct for BC years -> astronomical year, that is from year -1 to year 0
	if ( _year <= 0 ) {
		_year--;
	}

	const double frac = (_julian + 0.5) - floor(_julian+0.5); //the julian date reference is at 12:00
	_minute = ((int)round(frac*((double)24.0*60.0))) % 60;
	_hour = (int) round((((double)1440.0)*frac-(double)_minute)/(double)60.0);
}

bool Date::isLeapYear(const int& _year) const
{
	long jd1, jd2;
	jd1 = getJulianDayNumber( _year, 2, 28 );
	jd2 = getJulianDayNumber( _year, 3, 1 );
	return ( (jd2-jd1) > 1 );
}

long Date::getJulianDayNumber(const int& _year, const int& _month, const int& _day) const
{ //given year, month, day, calculate the matching julian day
 //see Fliegel, H. F. and van Flandern, T. C. 1968. Letters to the editor: a machine algorithm for processing calendar dates. Commun. ACM 11, 10 (Oct. 1968), 657. DOI= http://doi.acm.org/10.1145/364096.364097 
	const long lmonth = (long) _month, lday = (long) _day;
	long lyear = (long) _year;

	// Correct for BC years -> astronomical year, that is from year -1 to year 0
	if ( lyear < 0 ) {
		lyear++;
	}

	const long jdn = lday - 32075L +
		1461L * ( lyear + 4800L + ( lmonth - 14L ) / 12L ) / 4L +
		367L * ( lmonth - 2L - ( lmonth - 14L ) / 12L * 12L ) / 12L -
		3L * ( ( lyear + 4900L + ( lmonth - 14L ) / 12L ) / 100L ) / 4L;

	return jdn;
}

void Date::plausibilityCheck(const int& in_year, const int& in_month, const int& in_day, const int& in_hour, const int& in_minute) const {
	if ((in_year < -4713) || (in_year >3000)
	    || (in_month < 1) || (in_month > 12) 
	    || (in_day < 1) || ((in_day > daysNonLeapYear[in_month-1]) && !isLeapYear(in_year)) 
	    || ((in_day > daysLeapYear[in_month-1]) && isLeapYear(in_year)) 
	    || (in_hour < 0) || (in_hour > 24) 
	    || (in_minute < 0) || (in_minute > 59)) {
		throw IOException("InvalidDateFormat", AT);
	}

	if ((in_hour == 24) && (in_minute != 0)) {
		throw IOException("InvalidDateFormat", AT);
	}
}

double Date::localToGMT(const double& _julian) const {
	if(dst) {
		return (_julian - timezone/24. - DST_shift/24.);
	} else {
		return (_julian - timezone/24.);
	}
}

double Date::GMTToLocal(const double& _gmt_julian) const {
	if(dst) {
		return (_gmt_julian + timezone/24. + DST_shift/24.);
	} else {
		return (_gmt_julian + timezone/24.);
	}
}

#ifdef _POPC_
void Date::Serialize(POPBuffer &buf, bool pack)
{
	if (pack) {
		buf.Pack(&timezone,1);
		buf.Pack(&dst,1);
		buf.Pack(&gmt_julian,1);
		buf.Pack(&gmt_year,1);
		buf.Pack(&gmt_month,1);
		buf.Pack(&gmt_day,1);
		buf.Pack(&gmt_hour,1);
		buf.Pack(&gmt_minute,1);
		buf.Pack(&undef,1);
	} else {
		buf.UnPack(&timezone,1);
		buf.UnPack(&dst,1);
		buf.UnPack(&gmt_julian,1);
		buf.UnPack(&gmt_year,1);
		buf.UnPack(&gmt_month,1);
		buf.UnPack(&gmt_day,1);
		buf.UnPack(&gmt_hour,1);
		buf.UnPack(&gmt_minute,1);
		buf.UnPack(&undef,1);
	}
}
#endif

} //namespace
