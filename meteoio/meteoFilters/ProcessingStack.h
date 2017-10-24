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
#ifndef PROCESSINGSTACK_H
#define PROCESSINGSTACK_H

#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <meteoio/Config.h>
#include <memory>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  ProcessingStack
 * @brief This builds and runs through a filter stack for filtering a given parameter.
 * @author Thomas Egger
 * @date   2011-01-11
 */
class ProcessingStack {
	public:
		/**
		 * @brief Constructor parses cfg and builds up a filter stack for param_name
		 */
		ProcessingStack(const Config& cfg, const std::string& param_name);
		virtual ~ProcessingStack() {for (size_t ii=0; ii<filter_stack.size(); ii++) delete filter_stack[ii];}

		void process(const std::vector< std::vector<MeteoData> >& ivec,
		             std::vector< std::vector<MeteoData> >& ovec, const bool& second_pass=false);

		void getWindowSize(ProcessingProperties& o_properties) const;

		const std::string toString() const;
		
	private:
		virtual bool filterStation(std::vector<MeteoData> ivec, std::vector< std::vector<MeteoData> >& ovec, const bool& second_pass, const size_t& stat_idx);
		
	protected:
		virtual void addFilter(const std::string& block_name, const std::vector< std::pair<std::string, std::string> >& vecArgs, const Config& cfg);
		
		std::vector<ProcessingBlock*> filter_stack; //for now: strictly linear chain of processing blocks
		const std::string param_name;
};

/**
 * @class  TimeProcStack
 * @brief Since the time filters are quite specific to TIME (and need to be applied before), they have their own
 * ProcessingStack.
 */
class TimeProcStack : public ProcessingStack {
	public:
		TimeProcStack(const Config& cfg) : ProcessingStack(cfg, "TIME") {}
		virtual ~TimeProcStack() {for (size_t ii=0; ii<filter_stack.size(); ii++) delete filter_stack[ii];}

	private:
		virtual void addFilter(const std::string& block_name, const std::vector< std::pair<std::string, std::string> >& vecArgs, const Config& cfg);
		virtual bool filterStation(std::vector<MeteoData> ivec, std::vector< std::vector<MeteoData> >& ovec, const bool& second_pass, const size_t& stat_idx);
};

} //end namespace

#endif
