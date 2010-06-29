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
#include <meteoio/Meteo1DInterpolator.h>

using namespace std;

namespace mio {

Meteo1DInterpolator::Meteo1DInterpolator(const ConfigReader& _cfg) : cfg(_cfg) {
	/*
	 * By reading the ConfigReader object build up a list of user configured resampling algorithm
	 * for each MeteoData::Parameters parameter (i.e. each member variable like ta, p, hnw, ...)
	 * Concept of this constructor: loop over all MeteoData::Parameters and then look
	 * for configuration of resampling algorithms within the ConfigReader object.
	 */
	for (unsigned int ii=0; ii<MeteoData::nrOfParameters; ii++){ //loop over all MeteoData member variables
		const std::string& parname = MeteoData::getParameterName(ii); //Current parameter name

		vector<string> vecResamplingArguments;
		string resamplingAlgorithm = getInterpolationForParameter(parname, vecResamplingArguments);

		tasklist.push_back(resamplingAlgorithm);
		taskargs.push_back(vecResamplingArguments);
	}

	/*//For debugging only: 	
	for (unsigned int jj=0; jj<tasklist.size(); jj++){
		cout << MeteoData::getParameterName(jj) << "::" << tasklist[jj] << endl;
		for (unsigned int ii=0; ii<taskargs[jj].size(); ii++){
			cout << "\tARGS: " << taskargs[jj][ii] << endl;
		}
	}
	*/
}

unsigned int Meteo1DInterpolator::resampleData(const Date& date, std::vector<MeteoData>& vecM, std::vector<StationData>& vecS)
{
	if (vecM.size() != vecS.size())
		throw IOException("Inconsistency between vecM and vecS detected", AT);

	//Find element in the vector, or insert it at the appropriate position
	unsigned int position = IOUtils::seek(date, vecM, false);

	if (position == IOUtils::npos){ //nothing found append new element at the left or right
		if (vecM.size() == 0){
			vecM.push_back(MeteoData(date));
			vecS.push_back(StationData());
			return 0; //nothing left to do
		}

		if (vecM.at(0).date > date){
			vecM.insert(vecM.begin(), MeteoData(date));
			vecS.insert(vecS.begin(), vecS.at(0)); //copy element
			position = 0;
		} else if (vecM.at(vecM.size()-1).date < date){
			vecM.push_back(MeteoData(date));
			vecS.push_back(vecS.at(vecS.size()-1)); //copy element
			position = vecM.size() - 1;
		}
	} else if ((position != IOUtils::npos) && (vecM[position].date != date)){//insert before position
		vecM.insert(vecM.begin()+position, MeteoData(date));
		vecS.insert(vecS.begin()+position, vecS[position]); //copy element
	}

	for (unsigned int ii=0; ii<tasklist.size(); ii++){ //For all meteo parameters
		//cout << "For parameter: " << MeteoData::getParameterName(ii) << ": " << tasklist[ii] << endl;

		if (tasklist[ii] != "no") //resampling can be disabled by stating e.g. TA::resample = no
			ResamplingAlgorithms::getAlgorithm(tasklist[ii])(position, MeteoData::Parameters(ii), 
                                                                taskargs[ii], vecM, vecS);
	}

	return position; //the position of the resampled MeteoData object within vecM
}

string Meteo1DInterpolator::getInterpolationForParameter(const std::string& parname, std::vector<std::string>& vecArguments)
{
	/*
	 * This function retrieves the resampling algorithm to be used for the 
	 * 1D interpolation of meteo parameters. It also extracts any possible 
	 * arguments for that specific algorithm.
	 */

	vecArguments.clear();
	cfg.getValue(parname+"::args", "Interpolations1D", vecArguments, ConfigReader::nothrow);

	std::string tmp = "";
	cfg.getValue(parname+"::resample", "Interpolations1D", tmp, ConfigReader::nothrow);

	if (tmp.length() > 0)
		return tmp;

	return "linear"; //the default resampling is linear
}

} //namespace
