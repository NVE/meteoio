#ifndef __METEO1DRESAMPLER_H__
#define __METEO1DRESAMPLER_H__

#include "MeteoBuffer.h"
#include "Date_IO.h"
#include "libinterpol1D.h"
#include <string>
#include "IOUtils.h"

class Meteo1DResampler {

	public:
		Meteo1DResampler();
		void resample(const unsigned int& index_in, const Date_IO& date_in, 
				    std::vector<MeteoData>& mbuffer_out, std::vector<StationData>& sbuffer_out); 


	private:
		void seekIndices(MeteoBuffer& mbuffer, const std::string& parameter, unsigned int& leftindex, unsigned int& rightindex);
};

#endif
