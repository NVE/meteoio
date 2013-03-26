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

#include <meteoio/DataGenerator.h>

using namespace std;

namespace mio {

/**
*
*/ //explain how to use the generators for the end user

DataGenerator::DataGenerator(const Config& i_cfg)
              : cfg(i_cfg), mapAlgorithms(), generators_defined(false)
{
	setAlgorithms();
}

DataGenerator::~DataGenerator()
{ //we have to deallocate the memory allocated by "new GeneratorAlgorithm()"
	std::map< size_t, std::vector<GeneratorAlgorithm*> >::iterator it;
	for(it=mapAlgorithms.begin(); it!=mapAlgorithms.end(); it++) {
		std::vector<GeneratorAlgorithm*> &vec = it->second;
		for(size_t ii=0; ii<vec.size(); ii++)
			delete vec[ii];
	}
}

DataGenerator& DataGenerator::operator=(const DataGenerator& source)
{
	if(this != &source) {
		mapAlgorithms = source.mapAlgorithms;
		generators_defined = source.generators_defined;
	}
	return *this;
}

/**
 * @brief generate data to fill missing data points.
 * This relies on data generators defined by the user for each meteo parameters.
 * This loops over the defined generators and stops as soon as all missing points
 * have been successfully replaced.
 * @param vecMeteo vector containing one point for each station
 */
void DataGenerator::fillMissing(METEO_SET& vecMeteo) const
{
	if(!generators_defined) return; //no generators defined by the end user

	for(size_t param=0; param < MeteoData::nrOfParameters; param++) { //loop over all possible meteo parameters
		const std::map< size_t, std::vector<GeneratorAlgorithm*> >::const_iterator it = mapAlgorithms.find(param);

		if (it != mapAlgorithms.end()){ //there is a data generator for this parameter
			const std::vector<GeneratorAlgorithm*> vecGenerators = it->second;

			for(size_t station=0; station<vecMeteo.size(); station++) { //process this parameter on all stations
				size_t jj=0;
				while (jj<vecGenerators.size() && vecGenerators[jj]->generate(param, vecMeteo[station]) != true) jj++;
			}
		}
	}
}

/**
 * @brief generate data to fill missing data points.
 * This relies on data generators defined by the user for each meteo parameters.
 * This loops over the defined generators and stops as soon as all missing points
 * have been successfully replaced.
 * @param vecVecMeteo vector containing a timeserie for each station
 */
void DataGenerator::fillMissing(std::vector<METEO_SET>& vecVecMeteo) const
{
	if(!generators_defined) return; //no generators defined by the end user

	for(size_t param=0; param < MeteoData::nrOfParameters; param++) { //loop over all possible meteo parameters
		const std::map< size_t, std::vector<GeneratorAlgorithm*> >::const_iterator it = mapAlgorithms.find(param);

		if (it != mapAlgorithms.end()){ //there is a data generator for this parameter
			const std::vector<GeneratorAlgorithm*> vecGenerators = it->second;

			for(size_t station=0; station<vecVecMeteo.size(); station++) { //process this parameter on all stations
				size_t jj=0;
				while (jj<vecGenerators.size() && vecGenerators[jj]->generate(param, vecVecMeteo[station]) != true) jj++;
			}
		}
	}
}

/** @brief build the generators for each meteo parameter
 * By reading the Config object build up a list of user configured algorithms
 * for each MeteoData::Parameters parameter (i.e. each member variable of MeteoData like ta, p, hnw, ...)
 * Concept of this constructor: loop over all MeteoData::Parameters and then look
 * for configuration of interpolation algorithms within the Config object.
 */
void DataGenerator::setAlgorithms()
{
	for (size_t ii=0; ii < MeteoData::nrOfParameters; ii++){ //loop over all MeteoData member variables
		std::vector<std::string> tmpAlgorithms;
		const std::string& parname = MeteoData::getParameterName(ii); //Current parameter name
		const size_t nrOfAlgorithms = getAlgorithmsForParameter(parname, tmpAlgorithms);

		std::vector<GeneratorAlgorithm*> vecGenerators(nrOfAlgorithms);
		for(size_t jj=0; jj<nrOfAlgorithms; jj++) {
			std::vector<std::string> vecArgs;
			getArgumentsForAlgorithm(parname, tmpAlgorithms[jj], vecArgs);
			vecGenerators[jj] = GeneratorAlgorithmFactory::getAlgorithm( tmpAlgorithms[jj], vecArgs);
		}

		if(nrOfAlgorithms>0) {
			mapAlgorithms[ii] = vecGenerators;
			generators_defined = true;
		}
	}
}

size_t DataGenerator::getAlgorithmsForParameter(const std::string& parname, std::vector<std::string>& vecAlgorithms)
{
	// This function retrieves the user defined generator algorithms for
	// parameter 'parname' by querying the Config object
	vecAlgorithms.clear();

	std::vector<std::string> vecKeys;
	cfg.findKeys(vecKeys, parname+"::generators", "Generators");

	if (vecKeys.size() > 1)
		throw IOException("Multiple definitions of " + parname + "::generators in config file", AT);;

	if (vecKeys.empty())
		return 0;


	cfg.getValue(vecKeys.at(0), "Generators", vecAlgorithms, IOUtils::nothrow);

	return vecAlgorithms.size();
}

size_t DataGenerator::getArgumentsForAlgorithm(const std::string& parname,
                                               const std::string& algorithm,
                                               std::vector<std::string>& vecArgs) const
{
	vecArgs.clear();
	cfg.getValue(parname+"::"+algorithm, "Generators", vecArgs, IOUtils::nothrow);

	return vecArgs.size();
}

std::ostream& operator<<(std::ostream &os, const DataGenerator &mi) {
	os << "<DataGenerator>\n";
	os << "Config& cfg = " << hex << &mi.cfg << dec << "\n";

	os << "User list of generators:\n";
	std::map< size_t, std::vector<GeneratorAlgorithm*> >::const_iterator iter = mi.mapAlgorithms.begin();
	for (; iter != mi.mapAlgorithms.end(); ++iter) {
		os << setw(10) << MeteoData::getParameterName(iter->first) << " :: ";
		for(size_t jj=0; jj<iter->second.size(); jj++) {
			os << iter->second[jj]->getAlgo() << " ";
		}
		os << "\n";
	}

	os << "</DataGenerator>\n";
	return os;
}


#ifdef _POPC_
#include "marshal_meteoio.h"
using namespace mio; //HACK for POPC
void DataGenerator::Serialize(POPBuffer &buf, bool pack)
{
	/*if (pack)
	{

	}
	else
	{

	}*/
}
#endif

} //namespace
