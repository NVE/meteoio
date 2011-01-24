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

Meteo1DInterpolator::Meteo1DInterpolator(const Config& _cfg) : cfg(_cfg) {
	/*
	 * By reading the Config object build up a list of user configured resampling algorithm
	 * for each MeteoData::Parameters parameter (i.e. each member variable like ta, p, hnw, ...)
	 * Concept of this constructor: loop over all MeteoData::Parameters and then look
	 * for configuration of resampling algorithms within the Config object.
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

void Meteo1DInterpolator::getWindowSize(ProcessingProperties& o_properties)
{
	o_properties.points_before = 1;
	o_properties.points_after  = 1;

	o_properties.time_before   = Date(0.8);
	o_properties.time_after    = Date(0.8);
}

unsigned int Meteo1DInterpolator::resampleData(const Date& date, std::vector<MeteoData>& vecM)
{
	if (vecM.size() == 0){ //Deal with case of the empty vector
		vecM.push_back(MeteoData(date));
		return 0; //nothing left to do
	}

	//Find element in the vector, or insert it at the appropriate position
	unsigned int position = IOUtils::seek(date, vecM, false);

	MeteoData tmpmd(vecM.at(0)); //create a clone of one of the elements
	tmpmd.reset(); //set all values to IOUtils::nodata
	tmpmd.setDate(date);
	tmpmd.setResampled(true);

	if (position == IOUtils::npos){ //nothing found append new element at the left or right
		if (vecM.at(0).date > date){
			vecM.insert(vecM.begin(), tmpmd);
			position = 0;
		} else if (vecM.at(vecM.size()-1).date < date){
			vecM.push_back(tmpmd);
			position = vecM.size() - 1;
		}
	} else if ((position != IOUtils::npos) && (vecM[position].date != date)){//insert before position
		vecM.insert(vecM.begin()+position, tmpmd);
	}

	unsigned int ii = 0;

	for (ii=0; ii<tasklist.size(); ii++){ //For all meteo parameters
		//cout << "For parameter: " << MeteoData::getParameterName(ii) << ": " << tasklist[ii] << endl;

		if (tasklist[ii] != "no") //resampling can be disabled by stating e.g. TA::resample = no
			ResamplingAlgorithms::getAlgorithm(tasklist[ii])(position, ii, taskargs[ii], vecM);
	}

	//There might be more parameters, interpolate them too
	const MeteoData& origmd = vecM.at(0); //this element must exist at this point
	for ( ; ii < origmd.getNrOfParameters(); ii++){
		string parametername = origmd.getNameForParameter(ii);

		//In order to parse the user config only once for this parameter we store
		//the algorithm and its arguments in a hash map calles extended_tasklist
		map<string, pair<string, vector<string> > >::const_iterator it = extended_tasklist.find(parametername);
		if (it == extended_tasklist.end()){
			vector<string> taskarg; //vector to be filled with taskarguments
			const string algo = getInterpolationForParameter(parametername, taskarg);			
			
			extended_tasklist[parametername] = pair<string, vector<string> > (algo, taskarg);
			it = extended_tasklist.find(parametername);
		}

		if (it->second.first != "no") //resampling can be disabled by stating e.g. TA::resample = no
			ResamplingAlgorithms::getAlgorithm(it->second.first)(position, ii, it->second.second, vecM);
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
	cfg.getValue(parname+"::args", "Interpolations1D", vecArguments, Config::nothrow);

	std::string tmp = "";
	cfg.getValue(parname+"::resample", "Interpolations1D", tmp, Config::nothrow);

	if (tmp.length() > 0)
		return tmp;

	return "linear"; //the default resampling is linear
}

std::ostream& operator<<(std::ostream& os, const Meteo1DInterpolator& Interpolator) {

	os << "<Meteo1DInterpolator>\n";
	os << "Config cfg = " << hex << &Interpolator.cfg << "\n";
	for (unsigned int jj=0; jj<Interpolator.tasklist.size(); jj++){
		os << setw(10) << MeteoData::getParameterName(jj) << "::" << Interpolator.tasklist[jj] << "\t";
		for (unsigned int ii=0; ii<Interpolator.taskargs[jj].size(); ii++){
			os << "ARGS: " << Interpolator.taskargs[jj][ii] << " ";
		}
		os << "\n";
	}
	os << "</Meteo1DInterpolator>\n";

	return os;
}

} //namespace
