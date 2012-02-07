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
#include <meteoio/meteofilters/ProcessingStack.h>

using namespace std;

namespace mio {

ProcessingStack::ProcessingStack(const Config& cfg, const std::string& parname) : param_name(parname)
{
	vector<string> vec_filters;
	const size_t nr_of_filters = getFiltersForParameter(cfg, param_name, vec_filters);
	for (size_t ii=0; ii<nr_of_filters; ii++){
		//create a processing block for each filter
		string block_name = vec_filters[ii];
		IOUtils::toUpper(block_name);
		std::vector<std::string> vec_args;
		std::stringstream tmp;
		tmp << param_name << "::arg" << (ii+1);

		getArgumentsForFilter(cfg, tmp.str(), vec_args); //Read arguments
		filter_stack.push_back(BlockFactory::getBlock(block_name, vec_args));
	}
}

ProcessingStack::~ProcessingStack() {
	//this is suboptimal, shared_ptr<> would be the preference
	//it's unfortunately a part of boost only
	for (size_t ii=0; ii<filter_stack.size(); ii++)
		delete filter_stack[ii];
}

void ProcessingStack::getWindowSize(ProcessingProperties& o_properties) {
	o_properties.points_before = 0;
	o_properties.points_after = 0;
	o_properties.time_after = Duration(0.0, 0.);
	o_properties.time_before = Duration(0.0, 0.);

	for (size_t jj=0; jj<filter_stack.size(); jj++){
		const ProcessingProperties& properties = (*filter_stack[jj]).getProperties();

		o_properties.points_before = MAX(o_properties.points_before, properties.points_before);
		o_properties.points_after = MAX(o_properties.points_after, properties.points_after);

		if (properties.time_before > o_properties.time_before)
			o_properties.time_before = properties.time_before;

		if (properties.time_after > o_properties.time_after)
			o_properties.time_after = properties.time_after;
	}
}

size_t ProcessingStack::getFiltersForParameter(const Config& cfg, const std::string& parname, std::vector<std::string>& vecFilters)
{
	/*
	 * This function retrieves the filter sequence for parameter 'parname'
	 * by querying the Config object
	 */
	std::vector<std::string> vecKeys;
	std::string tmp;
	cfg.findKeys(vecKeys, parname+"::filter", "Filters");

	for (size_t ii=0; ii<vecKeys.size(); ii++){
		cfg.getValue(vecKeys[ii], "Filters", tmp, Config::nothrow);
		vecFilters.push_back(tmp);
	}

	return vecFilters.size();
}

size_t ProcessingStack::getArgumentsForFilter(const Config& cfg, const std::string& keyname,
                                                    std::vector<std::string>& vecArguments)
{
	// Retrieve the values for a given 'keyname' and store them in a vector calles 'vecArguments'
	cfg.getValue(keyname, "Filters", vecArguments, Config::nothrow);
	return vecArguments.size();
}

void ProcessingStack::process(const std::vector< std::vector<MeteoData> >& ivec,
                              std::vector< std::vector<MeteoData> >& ovec, const bool& second_pass)
{
	ovec.clear();
	ovec.insert(ovec.begin(), ivec.size(), vector<MeteoData>());

	for (size_t ii=0; ii<ivec.size(); ii++){ //for every station
		if (ivec[ii].size() > 0){
			//pick one element and check whether the param_name parameter exists
			const size_t index = ivec[ii][0].getParameterIndex(param_name);

			if (index != IOUtils::npos){
				std::vector<MeteoData> tmp = ivec[ii];

				//Now call the filters in a row
				bool appliedFilter = false;
				for (size_t jj=0; jj<filter_stack.size(); jj++){
					//cout << param_name << ": processing filter " << (*filter_stack[jj]).getName() << endl;
					if (second_pass){
						if ((*filter_stack[jj]).getProperties().stage==ProcessingProperties::first
						    || (*filter_stack[jj]).getProperties().stage==ProcessingProperties::none)
							continue;
					}
					if (!second_pass){
						if ((*filter_stack[jj]).getProperties().stage==ProcessingProperties::second
						    || (*filter_stack[jj]).getProperties().stage==ProcessingProperties::none)
							continue;
					}
					appliedFilter = true;

					(*filter_stack[jj]).process(index, tmp, ovec[ii]);

					if (tmp.size() == ovec[ii].size()){
						if ((jj+1) != filter_stack.size()){//after the last filter not necessary
							for (size_t jj=0; jj<ovec[ii].size(); jj++){
								tmp[jj](index) = ovec[ii][jj](index);
							}
						}
					} else {
						tmp = ovec[ii];
					}
				}

				if (!appliedFilter) //if not a single filter was applied
					ovec[ii] = ivec[ii]; //just copy input to output
			} else {
				ovec[ii] = ivec[ii]; //just copy input to output
			}
		}
	}
}

std::ostream& operator<<(std::ostream& os, const ProcessingStack& data) {
	//os << "<ProcessingStack>";
	os << setw(10) << data.param_name << "::";

	for(size_t ii=0; ii<data.filter_stack.size(); ii++) {
		os << setw(10) << *(data.filter_stack[ii]);
	}

	//os << "</ProcessingStack>";
	os << "\n";
	return os;
}

} //end namespace
