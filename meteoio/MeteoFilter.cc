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
#include <meteoio/MeteoFilter.h>

using namespace std;

namespace mio {

MeteoFilter::MeteoFilter(const ConfigReader& _cfg) : cfg(_cfg) {
	/*
	 * By reading the ConfigReader object build up a list of user configured filters
	 * for each MeteoData::Parameters parameter (i.e. each member variable like ta, p, hnw, ...)
	 * Concept of this constructor: loop over all MeteoData::Parameters and then look
	 * for configuration of filters within the ConfigReader object.
	 */
	for (unsigned int ii=0; ii<MeteoData::nrOfParameters; ii++){ //loop over all MeteoData member variables
		std::vector<std::string> tmpFilters1;
		std::vector<std::string> tmpFilters2;
		std::vector< std::vector<std::string> > parArgs; //Arguments for each filter

		const std::string& parname = MeteoData::getParameterName(ii); //Current parameter name

		unsigned int nrOfFilters = getFiltersForParameter(parname, tmpFilters2);

		for (unsigned int ll=0; ll<nrOfFilters; ll++){
			//Get the arguments for the specific filter from the cfg object
			std::vector<std::string> filterArgs;
			std::stringstream tmp;
			tmp << parname << "::arg" << (ll+1);
			getArgumentsForFilter(tmp.str(), filterArgs); //Read arguments
			//cout << "ARGSEARCH: " << tmp.str() << "  found arguments: " << argnum << endl;
			
			tmpFilters1.push_back(tmpFilters2[ll]);
			parArgs.push_back(filterArgs);
		}

		//cout << "ParArgsSize: " << parArgs.size() << endl;
		
		tasklist.push_back(tmpFilters1);
		taskargs.push_back(parArgs);
	}

	/* //For debugging only:
	for (unsigned int jj=0; jj<tasklist.size(); jj++){
		cout << MeteoData::getParameterName(jj) << "::" << endl;
		for (unsigned int ii=0; ii<tasklist[jj].size(); ii++){
			cout << tasklist[jj][ii] << "  ARGS: " << taskargs[jj][ii].size() << endl;
		}
	}
	*/
}

void MeteoFilter::filterData(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
                             std::vector<MeteoData> vecWindowM, std::vector<StationData> vecWindowS, 
                             const bool& checkonly)
{
	//Loop through each meteo parameter, call the respective filter function	
	for (unsigned int ii=0; ii<tasklist.size(); ii++){ //For all meteo parameters
		//cout << "For parameter: " << MeteoData::getParameterName(ii) << endl;
		for (unsigned int jj=0; jj<tasklist[ii].size(); jj++){ //For eack activated filter
			//Call the appropriate filter function
			const bool& isCheckOnly = FilterAlgorithms::filterProperties(tasklist[ii][jj]).checkonly;
			if (checkonly && !isCheckOnly)
				continue;

			//cout << "\tExecuting: " << tasklist[ii][jj] << endl;
			FilterAlgorithms::filterProperties(tasklist[ii][jj]).filterfunc(vecM, vecS, 
															    taskargs.at(ii).at(jj), 
															    MeteoData::Parameters(ii),
															    vecWindowM, vecWindowS);
		}
	}


}

unsigned int MeteoFilter::getFiltersForParameter(const std::string& parname, std::vector<std::string>& vecFilters)
{
	/* 
	 * This function retrieves the filter sequence for parameter 'parname' 
	 * by querying the ConfigReader object
	 */
	std::vector<std::string> vecKeys;
	std::string tmp;
	cfg.findKeys(vecKeys, parname+"::filter", "Filters");

	for (unsigned int ii=0; ii<vecKeys.size(); ii++){
		cfg.getValue(vecKeys[ii], "Filters", tmp, ConfigReader::nothrow);
		vecFilters.push_back(tmp);
	}

	return vecFilters.size();
}

unsigned int MeteoFilter::getArgumentsForFilter(const std::string& keyname, std::vector<std::string>& vecArguments)
{
	/*
	 * Retrieve the values for a given 'keyname' and store them in a vector calles 'vecArguments'
	 */
	cfg.getValue(keyname, "Filters", vecArguments, ConfigReader::nothrow);
	return vecArguments.size();
}

} //namespace
