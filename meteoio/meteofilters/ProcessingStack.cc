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

ProcessingStack::ProcessingStack(const Config& cfg, const std::string& parname) : filter_stack(), param_name(parname)
{
	vector<string> vecFilters;
	cfg.getValues(parname+"::filter", "Filters", vecFilters);
	const size_t nr_of_filters = vecFilters.size();
	for (size_t ii=0; ii<nr_of_filters; ii++){
		//create a processing block for each filter
		const string block_name = IOUtils::strToUpper( vecFilters[ii] );
		std::vector<std::string> vec_args;
		std::ostringstream tmp;
		tmp << param_name << "::arg" << (ii+1);

		getArgumentsForFilter(cfg, tmp.str(), vec_args); //Read arguments
		filter_stack.push_back(BlockFactory::getBlock(block_name, vec_args));
	}
}

ProcessingStack::~ProcessingStack()
{
	//this is suboptimal, shared_ptr<> would be the preference
	//it's unfortunately a part of boost only
	for (size_t ii=0; ii<filter_stack.size(); ii++)
		delete filter_stack[ii];
}

void ProcessingStack::getWindowSize(ProcessingProperties& o_properties)
{
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

size_t ProcessingStack::getArgumentsForFilter(const Config& cfg, const std::string& keyname,
                                                    std::vector<std::string>& vecArguments)
{
	// Retrieve the values for a given 'keyname' and store them in a vector calles 'vecArguments'
	cfg.getValue(keyname, "Filters", vecArguments, IOUtils::nothrow);
	return vecArguments.size();
}

void ProcessingStack::process(const std::vector< std::vector<MeteoData> >& ivec,
                              std::vector< std::vector<MeteoData> >& ovec, const bool& second_pass)
{
	ovec.resize( ivec.size() );

	for (size_t ii=0; ii<ivec.size(); ii++){ //for every station
		if (!ivec[ii].empty()){
			//pick one element and check whether the param_name parameter exists
			const size_t param = ivec[ii].front().getParameterIndex(param_name);
			if (param != IOUtils::npos){
				std::vector<MeteoData> tmp = ivec[ii];

				//Now call the filters in a row
				bool appliedFilter = false;
				for (size_t jj=0; jj<filter_stack.size(); jj++){
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

					(*filter_stack[jj]).process(param, tmp, ovec[ii]);

					if (tmp.size() == ovec[ii].size()){
						if ((jj+1) != filter_stack.size()){//after the last filter not necessary
							for (size_t kk=0; kk<ovec[ii].size(); kk++){
								tmp[kk](param) = ovec[ii][kk](param);
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

const std::string ProcessingStack::toString() const
{
	std::stringstream os;
	//os << "<ProcessingStack>";
	os << setw(10) << param_name << "::";

	for(size_t ii=0; ii<filter_stack.size(); ii++) {
		os << setw(10) << (*filter_stack[ii]).toString();
	}

	//os << "</ProcessingStack>";
	os << "\n";
	return os.str();
}

} //end namespace
