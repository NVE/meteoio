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

DataGenerator::DataGenerator(const Config& cfg)
              : mapAlgorithms(), generators_defined(false)
{
	setAlgorithms(cfg);
}

DataGenerator::~DataGenerator()
{ //we have to deallocate the memory allocated by "new GeneratorAlgorithm()"
	std::map< std::string, std::vector<GeneratorAlgorithm*> >::iterator it;
	for(it=mapAlgorithms.begin(); it!=mapAlgorithms.end(); ++it) {
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

	std::map< std::string, std::vector<GeneratorAlgorithm*> >::const_iterator it;
	for(it=mapAlgorithms.begin(); it!=mapAlgorithms.end(); ++it) {
		const std::vector<GeneratorAlgorithm*> vecGenerators = it->second;

		for(size_t station=0; station<vecMeteo.size(); station++) { //process this parameter on all stations
			const size_t param = vecMeteo[station].getParameterIndex(it->first);
			if(param==IOUtils::npos) continue;

			size_t jj=0;
			while (jj<vecGenerators.size() && vecGenerators[jj]->generate(param, vecMeteo[station]) != true) jj++;
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

	std::map< std::string, std::vector<GeneratorAlgorithm*> >::const_iterator it;
	for(it=mapAlgorithms.begin(); it!=mapAlgorithms.end(); ++it) {
		const std::vector<GeneratorAlgorithm*> vecGenerators = it->second;

		for(size_t station=0; station<vecVecMeteo.size(); station++) { //process this parameter on all stations
			const size_t param = vecVecMeteo[station][0].getParameterIndex(it->first);
			if(param==IOUtils::npos) continue;

			size_t jj=0;
			while (jj<vecGenerators.size() && vecGenerators[jj]->generate(param, vecVecMeteo[station]) != true) jj++;
		}
	}
}

/** @brief build the generators for each meteo parameter
 * By reading the Config object build up a list of user configured algorithms
 * for each MeteoData::Parameters parameter (i.e. each member variable of MeteoData like ta, p, hnw, ...)
 * Concept of this constructor: loop over all MeteoData::Parameters and then look
 * for configuration of interpolation algorithms within the Config object.
 * @param cfg configuration object to use for getting the algorithms configuration
 */
void DataGenerator::setAlgorithms(const Config& cfg)
{
	set<string> set_of_used_parameters;
	getParameters(cfg, set_of_used_parameters);

	set<string>::const_iterator it;
	for (it = set_of_used_parameters.begin(); it != set_of_used_parameters.end(); ++it) {
		std::vector<std::string> tmpAlgorithms;
		const std::string parname = *it;
		const size_t nrOfAlgorithms = getAlgorithmsForParameter(cfg, parname, tmpAlgorithms);

		std::vector<GeneratorAlgorithm*> vecGenerators(nrOfAlgorithms);
		for(size_t jj=0; jj<nrOfAlgorithms; jj++) {
			std::vector<std::string> vecArgs;
			getArgumentsForAlgorithm(cfg, parname, tmpAlgorithms[jj], vecArgs);
			vecGenerators[jj] = GeneratorAlgorithmFactory::getAlgorithm( tmpAlgorithms[jj], vecArgs);
		}

		if(nrOfAlgorithms>0) {
			mapAlgorithms[parname] = vecGenerators;
			generators_defined = true;
		}
	}
}

void DataGenerator::getParameters(const Config& cfg, std::set<std::string>& set_parameters)
{
	std::vector<std::string> vec_keys;
	cfg.findKeys(vec_keys, "::generators", "Generators", true); //search anywhere in key

	for (size_t ii=0; ii<vec_keys.size(); ii++) {
		const size_t found = vec_keys[ii].find_first_of(":");
		if (found != std::string::npos){
			const string tmp = vec_keys[ii].substr(0,found);
			set_parameters.insert( IOUtils::strToUpper(tmp) );
		}
	}
}

size_t DataGenerator::getAlgorithmsForParameter(const Config& cfg, const std::string& parname, std::vector<std::string>& vecAlgorithms)
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

size_t DataGenerator::getArgumentsForAlgorithm(const Config& cfg,
                                               const std::string& parname,
                                               const std::string& algorithm,
                                               std::vector<std::string>& vecArgs)
{
	vecArgs.clear();
	cfg.getValue(parname+"::"+algorithm, "Generators", vecArgs, IOUtils::nothrow);

	return vecArgs.size();
}

const std::string DataGenerator::toString() const {
	std::ostringstream os;
	os << "<DataGenerator>\n";
	os << "Generators defined: " << std::boolalpha << generators_defined << std::noboolalpha << "\n";
	os << "User list of generators:\n";
	std::map< std::string, std::vector<GeneratorAlgorithm*> >::const_iterator iter = mapAlgorithms.begin();
	for (; iter != mapAlgorithms.end(); ++iter) {
		os << setw(10) << iter->first << " :: ";
		for(size_t jj=0; jj<iter->second.size(); jj++) {
			os << iter->second[jj]->getAlgo() << " ";
		}
		os << "\n";
	}

	os << "</DataGenerator>\n";
	return os.str();
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
