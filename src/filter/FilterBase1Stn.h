#ifndef FILTERBASE1STN_H_INCLUDED
#define FILTERBASE1STN_H_INCLUDED

#include "FilterBase.h"

using namespace std;

/**
 * @class FilterBase1Stn
 * @brief Base class (interface) for data filters working on a single station.
 * @author Florian Hof
 * @date   2009-03-13
 */
class FilterBase1Stn : public FilterBase {

 public:

  // check handling

  /**
   * Check the meteo data at a single time and for a single station with the filter.
   * Before calling doCheck, the parameters have to be filled and the prepareCheck method called. 
   * @param unfilteredMeteoBuffer   [in] The buffer of unfiltered meteo data.
   * @param filteredMeteoBuffer   [in out] The buffer of filtered meteo data. Should already contain data
   *                              for the given timestamp. Filters should prefer using filteredMeteo
   *                              instead of unfilteredMeteo when possible, so that filters can be combined.
   * @param iFilteredElement   [in] The index in the filteredMeteoBuffer of the meteo data to filter,
   *                           thus specifying the date and time. 
   */
  virtual void doCheck(MeteoBuffer& unfilteredMeteoBuffer, MeteoBuffer& filteredMeteoBuffer, unsigned int iFilteredElement) = 0;

};

/** Type for functions that create concrete single-station filters (constructor-like). */
typedef FilterBase1Stn*(*Filter1StnCreator)(); 

/**
 * @class MeteoBufferIterator
 * @brief Iterator that provides meteo data in the past. It works over an unfiltered and a filtered meteo buffer.
 *        The iterator memorize the last returned data's time, so as to go in the past at each call of getPrevious. 
 * @author Florian Hof
 * @date   2009-03-16
 */
class MeteoBufferIterator {
  
 public: 
  
  /** 
   * Construct a new iterator that will provide meteo data in the past 
   * @param unfilteredMeteoBuffer   The unfiltered meteo buffer, as usually passed to the doCheck method. 
   * @param filteredMeteoBuffer   The filtered meteo buffer, as usually passed to the doCheck method. 
   * @param iFilteredElement   The index in the filteredMeteoBuffer of the current meteo data, 
   *                           thus specifying the date and time. 
   */
  MeteoBufferIterator(MeteoBuffer& unfilteredMeteoBuffer, MeteoBuffer& filteredMeteoBuffer, unsigned int iFilteredElement);
  
  /** Get the current (the last returned) meteo data. 
   * @throw NoAvailableDataException   When no data are available past the current date. 
   */
  MeteoData& getCurrent();
  
  /**
   * Get the previous meteo data, going backward at each call. 
   * Returns if possible data from the filtered buffer, otherwise from the unfiltered buffer. 
   * @throw NoAvailableDataException   When no data are available past the current date. 
   */
  MeteoData& getPrevious();
  
  /** 
   * Get the previous meteo data of the unfiltered buffer only, going backward at each call. 
   * @throw NoAvailableDataException   When no data are available past the current date. 
   */
  MeteoData& getPreviousUnfiltered();
  
 protected:
  
  MeteoBuffer& m_unfilteredMeteoBuffer;
  MeteoBuffer& m_filteredMeteoBuffer;
  unsigned int m_iUnfilteredElement;
  unsigned int m_iFilteredElement;
  
};

#endif
