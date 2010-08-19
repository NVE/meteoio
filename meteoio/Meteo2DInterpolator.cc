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

Meteo2DInterpolator::Meteo2DInterpolator(const Config& _cfg, const DEMObject& _dem, 
								 const std::vector<MeteoData>& _vecMeteo, 
								 const std::vector<StationData>& _vecStation) 
	: cfg(_cfg), dem(_dem), vecMeteo(_vecMeteo), vecStation(_vecStation)
{
	/*
	 * By reading the Config object build up a list of user configured algorithms
	 * for each MeteoData::Parameters parameter (i.e. each member variable of MeteoData like ta, p, hnw, ...)
	 * Concept of this constructor: loop over all MeteoData::Parameters and then look
	 * for configuration of interpolation algorithms within the Config object.
	 */
	for (unsigned int ii=0; ii < MeteoData::nrOfParameters; ii++){ //loop over all MeteoData member variables
		std::vector<std::string> tmpAlgorithms;
		const std::string& parname = MeteoData::getParameterName(ii); //Current parameter name
		unsigned int nrOfAlgorithms = getAlgorithmsForParameter(parname, tmpAlgorithms);

		if (nrOfAlgorithms > 0)
			mapAlgorithms[parname] = tmpAlgorithms;
	}

	//check whether the size of the two vectors is equal
	if (vecMeteo.size() != vecStation.size())
		throw IOException("Size of vector<MeteoData> and vector<StationData> are not equal", AT);

	//check that the stations are using the same projection as the dem
	for (unsigned int i=0; i<(unsigned int)vecStation.size(); i++) {
		if(!vecStation[i].position.isSameProj(dem.llcorner)) {
			throw IOException("Some stations are not using the same geographic projection as the DEM", AT);
		}
	}
}

void Meteo2DInterpolator::interpolate(const MeteoData::Parameters& meteoparam, Grid2DObject& result) const
{
	//Show algorithms to be used for this parameter
	map<string, vector<string> >::const_iterator it = mapAlgorithms.find(MeteoData::getParameterName(meteoparam));
	
	if (it != mapAlgorithms.end()){
		double maxQualityRating = 0.0;
		auto_ptr<InterpolationAlgorithm> bestalgorithm(NULL);
		vector<string> vecArgs;

		for (unsigned int ii=0; ii < it->second.size(); ii++){
			const string& algoname = it->second.at(ii);
			getArgumentsForAlgorithm(meteoparam, algoname, vecArgs);
			
			//Get the configured algorithm
			auto_ptr<InterpolationAlgorithm> algorithm(AlgorithmFactory::getAlgorithm(algoname, *this, dem, 
			                                            vecMeteo, vecStation, vecArgs));
			//Get the quality rating and compare to previously computed quality ratings
			const double rating = algorithm->getQualityRating(meteoparam);
			if ((rating != 0.0) && (rating > maxQualityRating)){
				//we use ">" so that in case of equality, the first choice will be kept
				bestalgorithm = algorithm; //remember this algorithm: ownership belongs to bestalgorithm
				maxQualityRating = rating;
			}
		}

		//finally execute the algorithm with the best quality rating or throw an exception
		if (bestalgorithm.get() == NULL) {
			throw IOException("No interpolation algorithm with quality rating >0 found", AT);
		}
		bestalgorithm->calculate(meteoparam, result);
	} else {
		//Some default message, that interpolation for this parameter needs configuration
		throw IOException("You need to configure the interpolation algorithms for parameter " + 
					   MeteoData::getParameterName(meteoparam), AT);
		
	}

	//check that the output grid is using the same projection as the dem
	if(!result.llcorner.isSameProj(dem.llcorner)) {
		throw IOException("The output grid is not using the same geographic projection as the DEM", AT);
	}

	//Run soft min/max filter for RH and HNW
	if (meteoparam == MeteoData::RH){
		Meteo2DInterpolator::checkMinMax(0.0, 1.0, result);
	} else if (meteoparam == MeteoData::HNW){
		Meteo2DInterpolator::checkMinMax(0.0, 10000.0, result);
	}
}

unsigned int Meteo2DInterpolator::getAlgorithmsForParameter(const std::string& parname, std::vector<std::string>& vecAlgorithms)
{
	/* 
	 * This function retrieves the user defined interpolation algorithms for 
	 * parameter 'parname' by querying the Config object
	 */
	std::vector<std::string> vecKeys;

	vecAlgorithms.clear();
	cfg.findKeys(vecKeys, parname+"::algorithms", "Interpolations2D");

	if (vecKeys.size() > 1)
		throw IOException("Multiple definitions of " + parname + "::algorithms in config file", AT);;

	if (vecKeys.size() == 0)
		return 0;


	cfg.getValue(vecKeys.at(0), "Interpolations2D", vecAlgorithms, Config::nothrow);

	return vecAlgorithms.size();
}

unsigned int Meteo2DInterpolator::getArgumentsForAlgorithm(const MeteoData::Parameters& param, 
                                                           const std::string& algorithm,
                                                           std::vector<std::string>& vecArgs) const
{
	vecArgs.clear();
	const string keyname = MeteoData::getParameterName(param) +"::"+ algorithm;
	cfg.getValue(keyname, "Interpolations2D", vecArgs, Config::nothrow);

	return vecArgs.size();
}

void Meteo2DInterpolator::checkMinMax(const double& minval, const double& maxval, Grid2DObject& gridobj)
{
	unsigned int nx=0, ny=0;

	gridobj.grid2D.size(nx, ny);

	for (unsigned int jj=0; jj<ny; jj++){
		for (unsigned int ii=0; ii<nx; ii++){
			double& value = gridobj.grid2D(ii,jj);
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
}


#ifdef _POPC_
#include "marshal_meteoio.h"
void Meteo2DInterpolator::Serialize(POPBuffer &buf, bool pack)
{
//TODO: check this serialization!! It seems dubious that it would work at all...
	if (pack)
	{
		//buf.Pack(&cfg,1);
		/*buf.Pack(&dem,1);
		marshal_METEO_DATASET(buf, vecMeteo, 0, FLAG_MARSHAL, NULL);
		marshal_STATION_DATASET(buf, vecStation, 0, FLAG_MARSHAL, NULL);*/
		marshal_map_str_vecstr(buf, mapAlgorithms, 0, FLAG_MARSHAL, NULL);
	}
	else
	{
		//buf.UnPack(&cfg,1);
		/*buf.UnPack(&dem,1);
		marshal_METEO_DATASET(buf, vecMeteo, 0, !FLAG_MARSHAL, NULL);
		marshal_STATION_DATASET(buf, vecStation, 0, !FLAG_MARSHAL, NULL);*/
		marshal_map_str_vecstr(buf, mapAlgorithms, 0, !FLAG_MARSHAL, NULL);
	}
}
#endif

} //namespace
