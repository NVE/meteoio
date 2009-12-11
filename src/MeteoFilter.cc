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
#include "MeteoFilter.h"

MeteoFilter::MeteoFilter(const ConfigReader& _cfg) : cfg(_cfg) {

	for (unsigned int ii=0; ii<MeteoData::nrOfParameters; ii++){
		std::vector<std::string> tmpFilters1;
		std::vector<std::string> tmpFilters2;
		std::vector< std::vector<std::string> > parArgs; //Arguments for each filter

		const std::string& parname = MeteoData::getParameterName(ii); //Current parameter name

		unsigned int nrOfFilters = getFiltersForParameter(parname, tmpFilters2);

		bool checkonly = false;
		for (unsigned int kk=0; kk<2; kk++){//Filtering occurs in two passes

			for (unsigned int ll=0; ll<nrOfFilters; ll++){
				//Get the arguments for the specific filter from the cfg object
				std::vector<std::string> filterArgs;
				std::stringstream tmp;
				tmp << parname << "::arg" << (ll+1);
				getArgumentsForFilter(tmp.str(), filterArgs); //Read arguments
				//cout << "ARGSEARCH: " << tmp.str() << "  found arguments: " << argnum << endl;

				if (tmpFilters2[ll] == "resample")
					throw InvalidArgumentException("Resampling not a valid filter", AT);
				
				if (checkonly){
					if (FilterAlgorithms::filterProperties(tmpFilters2[ll]).checkonly == true){
						tmpFilters1.push_back(tmpFilters2[ll]);
						parArgs.push_back(filterArgs);
					}
				} else {
					tmpFilters1.push_back(tmpFilters2[ll]);
					parArgs.push_back(filterArgs);
				}
			}

			//At the end of pass one of the filters, the resampling filter is added, unless disabled
			if (kk==0){
				std::string resamplingarg="";
				cfg.getValue(parname+"::resample", "Filters", resamplingarg, ConfigReader::nothrow);
				std::vector<std::string> filterArgs;
				getArgumentsForFilter(parname+"::resample", filterArgs); //Read arguments
				if (resamplingarg != "no"){
					tmpFilters1.push_back("resample");
					parArgs.push_back(filterArgs);
				
				}

				checkonly=true;
			}
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

bool MeteoFilter::filterData(const std::vector<MeteoData>& vecM, const std::vector<StationData>& vecS, 
					    const unsigned int& pos, const Date_IO& date,
					    MeteoData& md, StationData& sd)
{
	//No need to operate on the raw data, a copy of relevant data will be stored in these vectors:
	std::vector<MeteoData> vecFilteredM;   
	std::vector<StationData> vecFilteredS;

	if (vecM.at(pos).date == date){ //No resampling required, all other filters may apply
		vecFilteredM.push_back(vecM[pos]);
		vecFilteredS.push_back(vecS.at(pos));
	} else { 
		//resampling required
		vecFilteredM.push_back(vecM.at(pos-1));
		vecFilteredS.push_back(vecS.at(pos-1));
		vecFilteredM.push_back(vecM.at(pos));
		vecFilteredS.push_back(vecS.at(pos));
	}

	for (unsigned int ii=0; ii<tasklist.size(); ii++){ //For all meteo parameters
		//cout << "For parameter: " << MeteoData::getParameterName(ii) << endl;
		for (unsigned int jj=0; jj<tasklist[ii].size(); jj++){ //For eack activated filter
			//Call the appropriate filter function
			//cout << "\tExecuting: " << tasklist[ii][jj] << endl;
			if (!FilterAlgorithms::filterProperties(tasklist[ii][jj]).filterfunc(vecM, vecS, pos, date,
																    taskargs.at(ii).at(jj),
																    ii, vecFilteredM, vecFilteredS)){
				//break; //if one of the filters returns false, then stop filtering for this parameter
				;
			}
		}
	}

	if (vecFilteredM.size()==1){ //no resampling was required
		md = vecFilteredM[0];
		sd = vecFilteredS[0];		
	} else if (vecFilteredM.size()==3){ //resampling required and successfully executed
		md = vecFilteredM[1];
		sd = vecFilteredS[1];		
	} else { //Resampling was disabled, return a nodata data set
		md = MeteoData(date);
		sd = vecS.at(pos);
	}

	return true;
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

