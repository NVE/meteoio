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
#include "Date_IO.h"

//const long Date_IO::offset = 2415021;	///we define julian date as days since 1900/01/01, so we have an offset compared to std julian dates


using namespace std;

// see http://en.wikipedia.org/wiki/Julian_day for the various definitions of Julian dates

const int Date_IO::daysLeapYear[12] = {31,29,31,30,31,30,31,31,30,31,30,31};
const int Date_IO::daysNonLeapYear[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
const double Date_IO::DST_shift = 1.0;
const float Date_IO::MJD_offset = 2400000.5; ///offset between julian date and modified julian date
const float Date_IO::Unix_offset = 2440587.5; ///offset between julian date and Unix Epoch time
const float Date_IO::Excel_offset = 2415018.5;  ///offset between julian date and Excel dates (note that excel inveted some days...)

// CONSTUCTORS
Date_IO::Date_IO() {
	setDate(0., 0., false);
}

Date_IO::Date_IO(const double& julian_in, const double& _timezone, const bool& _dst) {
	setDate(julian_in, _timezone, _dst);
}

Date_IO::Date_IO(const time_t& _time, const double& _timezone, const bool& _dst) {
	setDate(_time, _timezone, _dst);
}
//HACK: is it needed? Why paroc_base instead of POPC??
#ifdef _POPC_
Date_IO::Date_IO(const Date_IO& _date_in) : paroc_base()
#else
Date_IO::Date_IO(const Date_IO& _date_in)
#endif
{
	setDate(_date_in.getJulianDate(), _date_in.getTimeZone(), _date_in.getDST());
}

Date_IO::Date_IO(const int& _year, const int& _month, const int& _day, const int& _hour, const int& _minute, const double& _timezone, const bool& _dst)
{
	setDate(_year, _month, _day, _hour, _minute, _timezone, _dst);
}

// SETTERS
void Date_IO::setTimeZone(const double& _timezone, const bool& _dst) {
//please keep in mind that timezone might be fractional (ie: 15 minutes, etc)
	if(abs(_timezone) > 12) {
		throw InvalidArgumentException("[E] Time zone can NOT be greater than +/-12!!", AT);
	}

	timezone = _timezone;
	dst = _dst;
}

void Date_IO::setDate(const int& _year, const int& _month, const int& _day, const int& _hour, const int& _minute, const double& _timezone, const bool& _dst)
{
	plausibilityCheck(_year, _month, _day, _hour, _minute);
	setTimeZone(_timezone, _dst);

	if(timezone!=0 || dst!=false) {
		//computing local julian date
		const double local_julian = calculateJulianDate(_year, _month, _day, _hour, _minute);
		//converting local julian date to GMT julian date
		gmt_julian = localToGMT(local_julian);
		//updating values to GMT
		calculateValues(gmt_julian, gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute);
	} else {
		//setting values and computing GMT julian date
		gmt_year = _year;
		gmt_month = _month;
		gmt_day = _day;
		gmt_hour = _hour;
		gmt_minute = _minute;
		gmt_julian = calculateJulianDate(gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute);
	}

}

void Date_IO::setDate(const double& julian_in, const double& _timezone, const bool& _dst) {
	setTimeZone(_timezone, _dst);
	gmt_julian = localToGMT(julian_in);
	calculateValues(gmt_julian, gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute);
}

void Date_IO::setDate(const time_t& _time, const double& _timezone, const bool& _dst) {
	setUnixDate(_time, _timezone, _dst);
}

void Date_IO::setModifiedJulianDate(const double& julian_in, const double& _timezone, const bool& _dst) {
	const double _julian = julian_in + MJD_offset;
	setDate(_julian, _timezone, _dst);
}

void Date_IO::setUnixDate(const time_t& _time, const double& _timezone, const bool& _dst) {
	const double _julian = (double)(_time)/(24.*60.*60.) + Unix_offset;
	setDate(_julian, _timezone, _dst);
}

void Date_IO::setExcelDate(const double julian_in, const double& _timezone, const bool& _dst) {
	const double _julian = julian_in + Excel_offset;
	setDate(_julian, _timezone, _dst);
}

// GETTERS
double Date_IO::getTimeZone() const {
	return timezone;
}

bool Date_IO::getDST() const {
	return dst;
}

double Date_IO::getJulianDate(const bool& gmt) const {
	if(gmt) {
		return gmt_julian;
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		return local_julian;
	}
}

double Date_IO::getModifiedJulianDate(const bool& gmt) const {
	if(gmt) {
		return (gmt_julian - MJD_offset);
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		return (local_julian - MJD_offset);
	}
}

//following NIST definition
double Date_IO::getTruncatedJulianDate(const bool& gmt) const {
	if(gmt) {
		return (fmod( (gmt_julian - 0.5), 10000. ));
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		return (fmod( (local_julian - 0.5), 10000. ));
	}
}

time_t Date_IO::getUnixDate(const bool& gmt) const {
	if (gmt_julian < Unix_offset)
			throw IOException("Dates before 1970 cannot be displayed in Unix epoch time", AT);

	if(gmt) {
		return ( (time_t)floor( (gmt_julian - Unix_offset) * (24*60*60) ));
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		return ( (time_t)floor( (local_julian - Unix_offset) * (24*60*60) ));
	}
}

double Date_IO::getExcelDate(const bool& gmt) const {
	if (gmt_julian < Excel_offset)
		throw IOException("Dates before 1900 cannot be converted to Excel date", AT);

	if(gmt) {
		return ( gmt_julian - Excel_offset);
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		return ( local_julian - Excel_offset);
	}
}

void Date_IO::getDate(double& julian_out, const bool& gmt) const {
	if(gmt) {
		julian_out = gmt_julian;
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		julian_out = local_julian;
	}
}

int Date_IO::getYear(const bool& gmt) const {
	if(gmt) {
		return gmt_year;
	} else {
		const double local_julian = GMTToLocal(gmt_julian);
		int local_year, local_month, local_day, local_hour, local_minute;
		calculateValues(local_julian, local_year, local_month, local_day, local_hour, local_minute);
		return local_year;
	}
}

void Date_IO::getDate(int& year_out, int& month_out, int& day_out, const bool& gmt) const {
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

void Date_IO::getDate(int& year_out, int& month_out, int& day_out, int& hour_out, const bool& gmt) const {
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
 
void Date_IO::getDate(int& year_out, int& month_out, int& day_out, int& hour_out, int& minute_out, const bool& gmt) const {
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


// OPERATORS
Date_IO& Date_IO::operator+=(const Date_IO& indate) {
	gmt_julian += indate.gmt_julian;
	calculateValues(gmt_julian, gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute);
	return *this;
}

Date_IO& Date_IO::operator-=(const Date_IO& indate) {
	gmt_julian -= indate.gmt_julian;
	calculateValues(gmt_julian, gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute);
	return *this;
}

bool Date_IO::operator==(const Date_IO& indate) const { 
	const double epsilon=1./(24.*3600.); //that is, 1 second in units of days
	return( fabs(indate.gmt_julian - gmt_julian) < epsilon );
}

bool Date_IO::operator!=(const Date_IO& indate) const {
	return !(*this==indate);
}

bool Date_IO::operator<(const Date_IO& indate) const {
	if (*this == indate) {
		return false;
	}

	return (gmt_julian < indate.gmt_julian);
}

bool Date_IO::operator<=(const Date_IO& indate) const {
	if (*this == indate) {
		return true;
	}
	return (gmt_julian <= indate.gmt_julian);
}

bool Date_IO::operator>(const Date_IO& indate) const {
	if (*this == indate) {
		return false;
	}

	return (gmt_julian > indate.gmt_julian);
}

bool Date_IO::operator>=(const Date_IO& indate) const {
	if (*this == indate) {
		return true;
	}
	return (gmt_julian >= indate.gmt_julian);
}

const Date_IO Date_IO::operator+(const Date_IO& indate) const {
	Date_IO tmp(gmt_julian + indate.gmt_julian);
	return tmp;
}

const Date_IO Date_IO::operator-(const Date_IO& indate) const {
	Date_IO tmp(gmt_julian - indate.gmt_julian);
	return tmp;
}

std::ostream& operator<<(std::ostream &os, const Date_IO &date) {
	os << "<date>\n";
	os << date.toString(Date_IO::ISO) << "\n";
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
	os << "</date>\n";
	return os;
}


// PRIVATE METHODS
double Date_IO::calculateJulianDate(const int& _year, const int& _month, const int& _day, const int& _hour, const int& _minute) const
{
	const long julday = getJulianDayNumber(_year, _month, _day);
	const double frac = (_hour-12.)/24. + _minute/(24.*60.); //the julian date reference is at 12:00

	return (((double)julday) + frac);
}

void Date_IO::calculateValues(const double& _julian, int& _year, int& _month, int& _day, int& _hour, int& _minute) const
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

bool Date_IO::isLeapYear(const int& iYear) const
{
	long jd1, jd2;
	jd1 = getJulianDayNumber( iYear, 2, 28 );
	jd2 = getJulianDayNumber( iYear, 3, 1 );
	return ( (jd2-jd1) > 1 );
}

long Date_IO::getJulianDayNumber(const int& _year, const int& _month, const int& _day) const
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

void Date_IO::plausibilityCheck(const int& in_year, const int& in_month, const int& in_day, const int& in_hour, const int& in_minute) const{
	if ((in_year < -4713) || (in_year >3000)
	    || (in_month < 1) || (in_month > 12) 
	    || (in_day < 1) || ((in_day > daysNonLeapYear[in_month-1]) && !isLeapYear(in_year)) 
	    || ((in_day > daysLeapYear[in_month-1]) && isLeapYear(in_year)) 
	    || (in_hour < 0) || (in_hour > 24) 
	    || (in_minute < 0) || (in_minute > 59)) {
		throw IOException("InvalidDate_IOFormat", AT);
	}

	if ((in_hour == 24) && (in_minute != 0)) {
		throw IOException("InvalidDate_IOFormat", AT);
	}
}

const string Date_IO::toString(FORMATS type, const bool& gmt) const
{//the date are displayed in LOCAL timezone (more user friendly)
	int year_out, month_out, day_out, hour_out, minute_out;
	double julian_out;

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
			<< setprecision(10) << julian_out << ")" ;
	} else {
		throw InvalidArgumentException("Wrong date conversion format requested", AT);
	}

	return tmpstr.str();
}

double Date_IO::localToGMT(const double& _julian) const {
	if(dst) {
		return (_julian - timezone/24. - DST_shift/24.);
	} else {
		return (_julian - timezone/24.);
	}
}

double Date_IO::GMTToLocal(const double& _gmt_julian) const {
	if(dst) {
		return (_gmt_julian + timezone/24. + DST_shift/24.);
	} else {
		return (_gmt_julian + timezone/24.);
	}
}

#ifdef _POPC_
void Date_IO::Serialize(POPBuffer &buf, bool pack)
{
	if (pack){
		TODO
		buf.Pack(&julian,1);
		buf.Pack(&year,1);
		buf.Pack(&month,1);
		buf.Pack(&day,1);
		buf.Pack(&hour,1);
		buf.Pack(&minute,1);
	}else{
		buf.UnPack(&julian,1);
		buf.UnPack(&year,1);
		buf.UnPack(&month,1);
		buf.UnPack(&day,1);
		buf.UnPack(&hour,1);
		buf.UnPack(&minute,1);
	}
}
#endif
