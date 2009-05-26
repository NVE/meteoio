#ifndef __DATE_H__
#define __DATE_H__

#include <string>
#include <sstream>
#include <utility>
#include <cmath>
#include <iomanip>
#include <iostream>

//#include <ctime>
//#include "slfutils.h"
#include "slfexceptions.h"

using namespace std; 
///Using the following namespace for the comparison operator overloading
using namespace rel_ops; 

/**
 * @class Date
 * @brief  A relatively basic class designed to represent julian dates and to provide basic conversions into yyyy/mm/dd/hh/mm representations
 *
 * The maximal precision as implemented is 1 minute.
 *
 * @author Thomas Egger
 */
#ifdef _PAROC_
class Date : POPBase {
  public:
    void Serialize(paroc_buffer &buf, bool pack);
#else
class Date {
#endif  
  public:
  static const int daysLeapYear[];
  static const int daysNonLeapYear[];
  static const long offset;

  ///Explicitly implemented Copy Constructor
  Date(const Date&); //Copy Constructor
  ///Note that constructing an Date object without any parameter results 
  ///in constructing Date(0.0) which in its current expression results in a date of 1900/1/1 00:00:00
  Date(const double& julian_in=0.0);
  ///All values passed will be checked for plausibility 
  Date(const int& year, const int& month, const int& day=1, const int& hour=0, const int& minute=0);
  //Date(const time_t&);

  void setDate(const double& julian_in);
  void setDate(const int& year, const int& month, const int& day=1, const int& hour=0, const int& minute=0);

  double getJulian() const;

  void getDate(double& julian_out) const;
  void getDate(int& year, int& month, int& day) const;
  void getDate(int& year, int& month, int& day, int& hour) const;
  void getDate(int& year, int& month, int& day, int& hour, int& minute) const;

  ///Since at SLF julian dates are always treated with an offset, this function provides a way to deal with real julian dates
  void setRealJulianDate(const double& julian_in);
  ///Since at SLF julian dates are always treated with an offset, this function provides a way to deal with real julian dates
  void getRealJulianDate(double& julian_out) const;

  ///The toString representation outputs the date in the form of "[double]julian_date yyyy/mm/dd hh:mm"
  const std::string toString(void) const;
  friend ostream& operator<<(ostream& os, const Date& date);

  //Operator Prototypes
  Date& operator=(const Date&);
  ///Can be used to add an interval to an existing Date object. 
  ///Construct a Date object representing the interval e.g. Date(1.0) for 1 day and add that to another Date object.
  Date& operator+=(const Date&);
  ///Can be used to subtract an interval from an existing Date object
  Date& operator-=(const Date&);
  bool operator==(const Date&) const;
  bool operator<(const Date&) const;

  const Date operator+(const Date&) const;
  const Date operator-(const Date&) const;
  
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
