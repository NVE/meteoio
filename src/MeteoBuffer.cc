#include "MeteoBuffer.h"

using namespace std;

MeteoBuffer::MeteoBuffer(const unsigned int& maxsize_in) : meteobuffer(0), stationbuffer(0){ //fixed size ring buffer
  maxsize = maxsize_in;
}

MeteoBuffer::MeteoBuffer(const unsigned int& maxsize_in, const unsigned int& initsize_in) : meteobuffer(0), stationbuffer(0){ //fixed size ring buffer
  maxsize = maxsize_in;

  if (initsize_in > maxsize)
    THROW IOException("In constructor of MeteoBuffer: cannot initilize more than maxsize_in objects", AT);

  MeteoData md;
  StationData sd;
  for (unsigned int jj=0; jj<initsize_in; jj++){
    put(md, sd);
  }
}

MeteoBuffer::~MeteoBuffer(){
  //Nothing here
}

void MeteoBuffer::clear(){
  meteobuffer.clear();
  stationbuffer.clear();
}

unsigned int MeteoBuffer::size(){
  return meteobuffer.size();
}

unsigned int MeteoBuffer::getMaxSize(){
  return maxsize;
}

MeteoData& MeteoBuffer::getMeteoData(const unsigned int& index){
  if (index >= meteobuffer.size())
    THROW IndexOutOfBoundsException("MeteoBuffer", AT);

  return meteobuffer[index];
}

StationData& MeteoBuffer::getStationData(const unsigned int& index){
  if (index >= stationbuffer.size())
    THROW IndexOutOfBoundsException("MeteoBuffer", AT);

  return stationbuffer[index];
}

void MeteoBuffer::insert(const unsigned int& index, MeteoData& meteo_in, StationData& station_in){
  if (index >= meteobuffer.size())
    THROW IndexOutOfBoundsException("Tried to insert into MeteoBuffer", AT);

  meteobuffer.insert(meteobuffer.begin()+index, meteo_in);
  stationbuffer.insert(stationbuffer.begin()+index, station_in);

  if (meteobuffer.size() > maxsize) {
    meteobuffer.pop_back();
    stationbuffer.pop_back();
  }
}

void MeteoBuffer::put(MeteoData& meteo_in, StationData& station_in){ //copied into buffer, always at the head
  meteobuffer.push_back(meteo_in);
  stationbuffer.push_back(station_in);

  if (meteobuffer.size() > maxsize) {
    meteobuffer.pop_front();
    stationbuffer.pop_front();
  }
}

void MeteoBuffer::erase(const unsigned int& index){
  if (index >= stationbuffer.size())
    THROW IndexOutOfBoundsException("MeteoBuffer", AT);
  
  meteobuffer.erase(meteobuffer.begin()+index, meteobuffer.end());
  stationbuffer.erase(stationbuffer.begin()+index, stationbuffer.end());
}

unsigned int MeteoBuffer::seek(const Date_IO& date_in){ //returns index of element, if element does not exist it returns closest index
  unsigned int ii = 1;

  if (size() <= 0) //no elements in buffer
    return MeteoBuffer::npos;

  //if we reach this point: at least one element in buffer
  if (meteobuffer[0].date > date_in)
    return MeteoBuffer::npos;

  if (meteobuffer[meteobuffer.size()-1].date < date_in) //last element is earlier, return npos
    return MeteoBuffer::npos;

  if (meteobuffer[0].date == date_in) //closest element
    return 0;

  //if we reach this point: the date is spanned by the buffer and there are at least two elements
  while ((ii < meteobuffer.size())){
    //cerr << "in search loop" << meteobuffer[ii].date.toString() << endl;
    if ((meteobuffer[ii].date >= date_in) && (meteobuffer[ii-1].date < date_in)) 
      return ii;
    
    ii++;
  }

  return MeteoBuffer::npos;
}

