#include "MeteoFilter.h"

using namespace std;

MeteoFilter::MeteoFilter(const ConfigReader& _cfg) : cfg(_cfg) {
	//Set up all the filters for each parameter
	//Init tasklist
	for (unsigned int ii=0; ii<MeteoData::nrOfParameters; ii++){
		vector<string> tmpFilters1;
		vector<string> tmpFilters2;
		vector< vector<string> > parArgs; //Arguments for each filter

		const string& parname = MeteoData::getParameterName(ii); //Current parameter name

		unsigned int nrOfFilters = getFiltersForParameter(parname, tmpFilters2);

		bool checkonly = false;
		for (unsigned int kk=0; kk<2; kk++){//Filtering occurs in two passes

			for (unsigned int ll=0; ll<nrOfFilters; ll++){
				//Get the arguments for the specific filter from the cfg object
				vector<string> filterArgs;
				stringstream tmp;
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
				string resamplingarg="";
				cfg.getValue(parname+"::resample", "Filters", resamplingarg, ConfigReader::nothrow);
				vector<string> filterArgs;
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

bool MeteoFilter::filterData(const vector<MeteoData>& vecM, const vector<StationData>& vecS, 
					    const unsigned int& pos, const Date_IO& date,
					    MeteoData& md, StationData& sd)
{
	//No need to operate on the raw data, a copy of relevant data will be stored in these vectors:
	vector<MeteoData> vecFilteredM;   
	vector<StationData> vecFilteredS;

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
		for (unsigned int jj=0; jj<tasklist[ii].size(); jj++){ //For eack activated filter
			//Call the appropriate filter function
			FilterAlgorithms::filterProperties(tasklist[ii][jj]).filterfunc(vecM,vecS,pos,date,taskargs.at(ii).at(jj),ii,vecFilteredM,vecFilteredS);
		}
	}

	if (vecFilteredM.size()==1){
		md = vecFilteredM[0];
		sd = vecFilteredS[0];		
	} else if (vecFilteredM.size()==3){
		md = vecFilteredM[1];
		sd = vecFilteredS[1];		
	} else {
		md = MeteoData(date);
		sd = vecS.at(pos);
	}

	return true;
}

unsigned int MeteoFilter::getFiltersForParameter(const string& parname, vector<string>& vecFilters)
{
	//get the vectors and associate them with the correct algorithms
	vector<string> vecKeys;
	string tmp;
	cfg.findKeys(vecKeys, parname+"::filter", "Filters");

	for (unsigned int ii=0; ii<vecKeys.size(); ii++){
		cfg.getValue(vecKeys[ii], "Filters", tmp, ConfigReader::nothrow);
		vecFilters.push_back(tmp);
	}

	return vecFilters.size();
}

unsigned int MeteoFilter::getArgumentsForFilter(const string& keyname, vector<string>& vecArguments)
{
	cfg.getValue(keyname, "Filters", vecArguments, ConfigReader::nothrow);
	return vecArguments.size();
}

