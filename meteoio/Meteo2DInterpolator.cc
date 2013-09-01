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

#include <meteoio/Meteo2DInterpolator.h>

using namespace std;

namespace mio {

Meteo2DInterpolator::Meteo2DInterpolator(const Config& i_cfg, IOManager& i_iom)
                    : cfg(i_cfg), iomanager(&i_iom), algorithms_ready(false), mapAlgorithms()
{
	setAlgorithms();
}

Meteo2DInterpolator::Meteo2DInterpolator(const Config& i_cfg)
                    : cfg(i_cfg), iomanager(NULL), algorithms_ready(false), mapAlgorithms() {}

Meteo2DInterpolator::~Meteo2DInterpolator()
{
	std::map<std::string, std::vector<InterpolationAlgorithm*> >::iterator iter;
	for (iter = mapAlgorithms.begin(); iter != mapAlgorithms.end(); ++iter) {
		const vector<InterpolationAlgorithm*>& vecAlgs = iter->second;
		for(size_t ii=0; ii<vecAlgs.size(); ++ii)
			delete vecAlgs[ii];
	}
}

Meteo2DInterpolator& Meteo2DInterpolator::operator=(const Meteo2DInterpolator& source)
{
	//since this uses an IOManager on a given machine/node, since the pointers point to entry points
	//in the compiled code, they should remain valid and therefore can be copied
	if(this != &source) {
		//cfg: can not be copied
		iomanager = source.iomanager;
		mapAlgorithms = source.mapAlgorithms;
	}
	return *this;
}

void Meteo2DInterpolator::setIOManager(IOManager& i_iomanager) {
	iomanager = &i_iomanager;
}

/* By reading the Config object build up a list of user configured algorithms
* for each MeteoData::Parameters parameter (i.e. each member variable of MeteoData like ta, p, hnw, ...)
* Concept of this constructor: loop over all MeteoData::Parameters and then look
* for configuration of interpolation algorithms within the Config object.
*/
void Meteo2DInterpolator::setAlgorithms()
{
	set<string> set_of_used_parameters;
	get_parameters(cfg, set_of_used_parameters);

	set<string>::const_iterator it;
	for (it = set_of_used_parameters.begin(); it != set_of_used_parameters.end(); ++it) {
		const std::string parname = *it;
		std::vector<std::string> tmpAlgorithms;
		const size_t nrOfAlgorithms = getAlgorithmsForParameter(cfg, parname, tmpAlgorithms);

		std::vector<InterpolationAlgorithm*> vecAlgorithms(nrOfAlgorithms);
		for(size_t jj=0; jj<nrOfAlgorithms; jj++) {
			std::vector<std::string> vecArgs;
			getArgumentsForAlgorithm(parname, tmpAlgorithms[jj], vecArgs);
			vecAlgorithms[jj] = AlgorithmFactory::getAlgorithm( tmpAlgorithms[jj], *this, vecArgs, *iomanager);
		}

		if(nrOfAlgorithms>0) {
			mapAlgorithms[parname] = vecAlgorithms;
		}
	}
	algorithms_ready = true;
}

//get a list of all meteoparameters referenced in the Interpolations2D section
size_t Meteo2DInterpolator::get_parameters(const Config& cfg, std::set<std::string>& set_parameters)
{
	std::vector<std::string> vec_keys;
	cfg.findKeys(vec_keys, std::string(), "Interpolations2D");

	for (size_t ii=0; ii<vec_keys.size(); ii++) {
		const size_t found = vec_keys[ii].find_first_of(":");
		if (found != std::string::npos){
			const string tmp = vec_keys[ii].substr(0,found);
			set_parameters.insert( IOUtils::strToUpper(tmp) );
		}
	}

	return set_parameters.size();
}

void Meteo2DInterpolator::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                                      Grid2DObject& result)
{
	std::string InfoString;
	interpolate(date, dem, meteoparam, result, InfoString);
}

void Meteo2DInterpolator::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                                      Grid2DObject& result, std::string& InfoString)
{
	if(iomanager==NULL)
		throw IOException("No IOManager reference has been set!", AT);
	if(!algorithms_ready)
		setAlgorithms();

	//Show algorithms to be used for this parameter
	const string param_name = MeteoData::getParameterName(meteoparam);
	const map<string, vector<InterpolationAlgorithm*> >::iterator it = mapAlgorithms.find(param_name);
	if(it==mapAlgorithms.end()) {
		throw IOException("No interpolation algorithms configured for parameter "+param_name, AT);
	}

	//look for algorithm with the highest quality rating
	const vector<InterpolationAlgorithm*>& vecAlgs = it->second;
	double maxQualityRating = -1.;
	size_t bestalgorithm = 0;
	for (size_t ii=0; ii < vecAlgs.size(); ++ii){
		const double rating = vecAlgs[ii]->getQualityRating(date, meteoparam);
		if ((rating != 0.0) && (rating > maxQualityRating)) {
			//we use ">" so that in case of equality, the first choice will be kept
			bestalgorithm = ii;
			maxQualityRating = rating;
		}
	}

	//finally execute the algorithm with the best quality rating or throw an exception
	if(maxQualityRating<=0.0)
		throw IOException("No interpolation algorithm with quality rating >0 found for parameter "+param_name, AT);
	vecAlgs[bestalgorithm]->calculate(dem, result);
	InfoString = vecAlgs[bestalgorithm]->getInfo();

	//check that the output grid is using the same projection as the dem
	if(!result.llcorner.isSameProj(dem.llcorner)) {
		throw IOException("The output grid is not using the same geographic projection as the DEM", AT);
	}

	//Run soft min/max filter for RH, HNW and HS
	if (meteoparam == MeteoData::RH){
		Meteo2DInterpolator::checkMinMax(0.0, 1.0, result);
	} else if (meteoparam == MeteoData::HNW){
		Meteo2DInterpolator::checkMinMax(0.0, 10000.0, result);
	} else if (meteoparam == MeteoData::HS){
		Meteo2DInterpolator::checkMinMax(0.0, 10000.0, result);
	} else if (meteoparam == MeteoData::VW){
		Meteo2DInterpolator::checkMinMax(0.0, 10000.0, result);
	}
}

size_t Meteo2DInterpolator::getAlgorithmsForParameter(const Config& cfg, const std::string& parname, std::vector<std::string>& vecAlgorithms)
{
	// This function retrieves the user defined interpolation algorithms for
	// parameter 'parname' by querying the Config object
	vecAlgorithms.clear();
	std::vector<std::string> vecKeys;
	cfg.findKeys(vecKeys, parname+"::algorithms", "Interpolations2D");

	if (vecKeys.size() > 1)
		throw IOException("Multiple definitions of " + parname + "::algorithms in config file", AT);;

	if (vecKeys.empty())
		return 0;

	cfg.getValue(vecKeys[0], "Interpolations2D", vecAlgorithms, IOUtils::nothrow);
	return vecAlgorithms.size();
}

size_t Meteo2DInterpolator::getArgumentsForAlgorithm(const std::string& param,
                                                     const std::string& algorithm,
                                                     std::vector<std::string>& vecArgs) const
{
	vecArgs.clear();
	const string keyname = param +"::"+ algorithm;
	cfg.getValue(keyname, "Interpolations2D", vecArgs, IOUtils::nothrow);

	return vecArgs.size();
}

void Meteo2DInterpolator::checkMinMax(const double& minval, const double& maxval, Grid2DObject& gridobj)
{
	const unsigned int nxy = gridobj.getNx() * gridobj.getNy();

	for (unsigned int ii=0; ii<nxy; ii++){
		double& value = gridobj.grid2D(ii);
		if (value == IOUtils::nodata){
			continue;
		}
		if (value < minval) {
			value = minval;
		} else if (value > maxval) {
			value = maxval;
		}
	}
}

void Meteo2DInterpolator::check_projections(const DEMObject& dem, const std::vector<MeteoData>& vec_meteo)
{
	//check that the stations are using the same projection as the dem
	for (size_t ii=0; ii<vec_meteo.size(); ii++) {
		const StationData& meta = vec_meteo[ii].meta;
		if(!meta.position.isSameProj(dem.llcorner)) {
			std::stringstream os;
			std::string type, args;
			meta.position.getProj(type, args);
			os << "Station " << meta.stationID << " is using projection (" << type << " " << args << ") ";
			dem.llcorner.getProj(type, args);
			os << "while DEM is using projection ("<< type << " " << args << ") ";
			throw IOException(os.str(), AT);
		}
	}
}


const std::string Meteo2DInterpolator::toString() const {
	ostringstream os;
	os << "<Meteo2DInterpolator>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";
	os << "IOManager& iomanager = "  << hex << &iomanager << dec << "\n";

	os << "User list of algorithms:\n";
	std::map<std::string, std::vector<InterpolationAlgorithm*> >::const_iterator iter;
	for (iter = mapAlgorithms.begin(); iter != mapAlgorithms.end(); ++iter) {
		os << setw(10) << iter->first << " :: ";
		for(unsigned int jj=0; jj<iter->second.size(); jj++) {
			//os << iter->second[jj]-> << " ";
		}
		os << "\n";
	}

	os << "</Meteo2DInterpolator>\n";
	return os.str();
}


#ifdef _POPC_
#include "marshal_meteoio.h"
using namespace mio; //HACK for POPC
void Meteo2DInterpolator::Serialize(POPBuffer &buf, bool pack)
{
	/*if (pack)
	{
		Config cfg2(cfg);
		cfg2.Serialize(buf, true); //cfg will now be a full copy
		size_t nr=mapAlgorithms.size();
		buf.Pack(&nr, 1);
		for(std::map< std::string, std::vector<std::string> >::const_iterator it = mapAlgorithms.begin(); it != mapAlgorithms.end(); ++it) {
			buf.Pack(&(it->first), 1); //param name
			const size_t n = (it->second).size();
			buf.Pack(&n, 1);
			for(size_t ii=0; ii<n; ii++) buf.Pack(&(it->second[ii]), 1); //vector of algorithms' names
		}
	}
	else
	{
		cfg.Serialize(buf, false); //cfg will now be a full copy
		iomanager = new IOManager(cfg);

		size_t nr=0;
		std::string key;
		std::vector<std::string> value;
		buf.UnPack(&nr,1);
		mapAlgorithms.clear();
		for(size_t i=0; i<nr; i++) {
			buf.UnPack(&key,1);
			size_t vec_n;
			buf.UnPack(vec_n, 1);
			value.resize(vec_n);
			for(size_t ii=0; ii<n; ii++) buf.UnPack(&value[ii], 1); //vector of algorithms' names
			mapAlgorithms[key] = value;
		}

	}*/
}
#endif

} //namespace
