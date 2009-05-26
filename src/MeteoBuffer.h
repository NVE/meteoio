#ifndef __METEOBUFFER_H__
#define __METEOBUFFER_H__

#include "Date.h"
#include "MeteoData.h"
#include "StationData.h"

#include <iostream>
#include <string>
#include <deque>

/**
 * @class MeteoBuffer
 * @brief A class to hold buffered MeteoData and StationData objects (corresponding)
 *        and allow efficient insert and delete operations. The underlying implementation is
 *        a ring buffer. The maximum size of the ring buffer is defined during construction.
 *
 * @author Thomas Egger
 * @date   2009-03-11
 */
class MeteoBuffer {

 public:

  /**
   * @brief A constructor that also sets the maximum size of the ring buffer
   * @param size unsigned int setting maximum size of internal buffer
   */  
  MeteoBuffer(const unsigned int& size); //fixed size ring buffer
  MeteoBuffer(const unsigned int& maxsize_in, const unsigned int& initsize_in); 
  ~MeteoBuffer();

  /**
   * @brief Delete all elements in ring buffer - effective size of buffer after clear() is 0
   *
   */  
  void clear();

  /**
   * @brief Delete a specific element with a specific index from the ring buffer
   * @param index unsigned int the index of the element to be deleted
   *
   */  
  void erase(const unsigned int& index);

  /**
   * @brief Return actual size of ring buffer
   * @return unsigned int representing actual size of the ring buffer
   *
   */  
  unsigned int size();

  unsigned int getMaxSize(void);


  //void get(const unsigned int& index, vector<MeteoData>* meteo_out, vector<StationData>* station_out);

  /**
   * @brief return a reference to buffered MeteoData object at specified index in the ring buffer
   * @param index unsigned int the index of the element to be referenced
   * @return MeteoData& reference
   *
   */  
  MeteoData& getMeteoData(const unsigned int& index);

  /**
   * @brief return a reference to buffered StationData object at specified index in the ring buffer
   * @param index unsigned int the index of the element to be referenced
   * @return StationData& reference
   *
   */  
  StationData& getStationData(const unsigned int& index);

  /**
   * @brief return a reference to buffered StationData object at specified index in the ring buffer
   * @param index unsigned int the index of the element before which the new element is to be inserted
   * @param meteo_in data to be copied into ring buffer
   * @param station_in data to be copied into ring buffer
   *
   */  
  void insert(const unsigned int& index, MeteoData& meteo_in, StationData& station_in);

  /**
   * @brief append data at the head of the ring buffer
   * @param meteo_in data to be copied into ring buffer
   * @param station_in data to be copied into ring buffer
   *
   */  
  void put(MeteoData& meteo_in, StationData& station_in); //copied into buffer, always at the head

  /**
   * @brief Find the index of the MeteoData object in the ring buffer, that matches the date_in parameter best
   * @param date_in Date of MeteoData object
   * @return unsigned int that represents the index of the element with the specified date, or the index closest to the Date
   *         if no element is found it returns MeteoBuffer::npos
   *
   */  
  unsigned int seek(const Date& date_in); //returns index of element, if element does not exist in 


  static const unsigned int npos = (unsigned int)-1;             ///<npos is the out-of-range value
  
 private:
  std::deque<MeteoData> meteobuffer;     ///<The ring buffer for MeteoData
  std::deque<StationData> stationbuffer; ///<The ring buffer for StationData
  unsigned int start, end, maxsize;
};

#endif
