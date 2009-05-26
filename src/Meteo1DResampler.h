#ifndef __METEO1DRESAMPLER_H__
#define __METEO1DRESAMPLER_H__

#include "MeteoBuffer.h"
#include "Date.h"
#include "libinterpol1D.h"
#include <string>

class Meteo1DResampler {

 public:
  Meteo1DResampler();
  void resample(const unsigned int& index_in, const Date& date_in, MeteoBuffer& mbuffer_out); 


 private:
  void seekIndices(MeteoBuffer& mbuffer, const std::string& parameter, unsigned int& leftindex, unsigned int& rightindex);
};

#endif
