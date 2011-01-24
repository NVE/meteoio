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
#include <meteoio/MeteoProcessor.h>

using namespace std;

namespace mio {

const unsigned int MeteoProcessor::window_half_size = 40; //org: 4

MeteoProcessor::MeteoProcessor(const Config& cfg) : mf(cfg), mi1d(cfg) 
{
	//Parse [Filters] section, create processing stack for each configured parameter
	set<string> set_of_used_parameters;
	//unsigned int nr_of_parameters = 
	get_parameters(cfg, set_of_used_parameters);

	for (set<string>::const_iterator it = set_of_used_parameters.begin(); it != set_of_used_parameters.end(); it++){
		//cout << "Creating stack for parameter: " << *it << endl;
		ProcessingStack* tmp = new ProcessingStack(cfg, *it);
		processing_stack[*it] = tmp;
	}
	//cout << "NrOfParameters: " << nr_of_parameters << endl;
}

MeteoProcessor::~MeteoProcessor()
{
	//clean up heap memory
	for (map<string, ProcessingStack*>::const_iterator it=processing_stack.begin(); it != processing_stack.end(); it++)
		delete it->second;
}

unsigned int MeteoProcessor::get_parameters(const Config& cfg, std::set<std::string>& set_parameters)
{
	std::vector<std::string> vec_keys;
	cfg.findKeys(vec_keys, "", "Filters");

	for (unsigned int ii=0; ii<vec_keys.size(); ii++){
		size_t found = vec_keys[ii].find_first_of(":");
		if (found != std::string::npos){
			string tmp = vec_keys[ii].substr(0,found);
			set_parameters.insert(tmp);
		}
	}

	return set_parameters.size();
}

void MeteoProcessor::getWindowSize(ProcessingProperties& o_properties)
{
	ProcessingProperties tmp;

	for (map<string, ProcessingStack*>::const_iterator it=processing_stack.begin(); it != processing_stack.end(); it++){
		(*(it->second)).getWindowSize(tmp);

		compareProperties(tmp, o_properties);
	}

	//Also take the Meteo1DInterpolator into account:
	mi1d.getWindowSize(tmp);
	compareProperties(tmp, o_properties);
}

void MeteoProcessor::compareProperties(const ProcessingProperties& newprop, ProcessingProperties& current)
{
	current.points_before = MAX(current.points_before, newprop.points_before);
	current.points_after = MAX(current.points_after, newprop.points_after);
	
	if (newprop.time_before > current.time_before)
		current.time_before = newprop.time_before;
	
	if (newprop.time_after > current.time_after)
		current.time_after = newprop.time_after;
}

void MeteoProcessor::process(const std::vector< std::vector<MeteoData> >& ivec, 
                             std::vector< std::vector<MeteoData> >& ovec, const bool& second_pass)
{
	//call the different processing stacks
	std::vector< std::vector<MeteoData> > vec_tmp = ivec;
	//cout << "Calling processing stacks for all parameters" << endl;
	for (map<string, ProcessingStack*>::const_iterator it=processing_stack.begin(); it != processing_stack.end(); it++){
		//cout << "Calling processing stack for parameter " << it->first << endl;
		(*(it->second)).process(vec_tmp, ovec, second_pass);
		vec_tmp = ovec;
	}

	if (processing_stack.size() == 0)
		ovec = ivec;
}

unsigned int MeteoProcessor::resample(const Date& date, std::vector<MeteoData>& ivec)
{
	return mi1d.resampleData(date, ivec);
}

void MeteoProcessor::processData(const Date& date, const std::vector<MeteoData>& vecM, MeteoData& md)
{
	unsigned int currentpos = IOUtils::seek(date, vecM, false); 
	unsigned int startindex = IOUtils::npos, endindex = IOUtils::npos;
	
	//No need to operate on the raw data, a copy of relevant data will be stored in these vectors:
	std::vector<MeteoData> vecWindowM;

	/*
	 * Cut out a window of data, on which the filtering and the resampling will occur
	 */
	bool windowexists = false;
	for (int ii=(int)(currentpos-window_half_size-1); ii<=(int)(currentpos+window_half_size); ii++){
		if ((ii>=0) && (ii<(int)vecM.size())){
			
			if (!windowexists){
				windowexists = true;
				startindex = ii;
			}
			endindex = ii;

			vecWindowM.push_back(vecM.at(ii));
			//cout << "Added " << vecM[ii].date.toString(Date::ISO) << endl;
		}
	}

	mf.filterData(vecM, vecWindowM, false); //first pass

	unsigned int position = mi1d.resampleData(date, vecWindowM); //resampling

	mf.filterData(vecM, vecWindowM, true); //checkonly, second filter pass

	md = vecWindowM[position];
}

std::ostream& operator<<(std::ostream& os, const MeteoProcessor& data)
{
	os << "<MeteoProcessor>\n";
	os << data.mf;
	os << data.mi1d;
	os << "window_half_size = " << data.window_half_size << "\n";
	os << "</MeteoProcessor>\n";
	return os;
}

} //namespace
