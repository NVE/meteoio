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
#ifndef __PROCESSINGBLOCK_H__
#define __PROCESSINGBLOCK_H__

#include <meteoio/MeteoData.h>
#include <vector>
#include <string>
#include <set>

namespace mio {

class ProcessingProperties {
	public:
		ProcessingProperties() : time_before(0., 0.), time_after(0., 0.),
		                         points_before(0), points_after(0),
		                         for_second_pass(false) {}

		friend std::ostream& operator<<(std::ostream& os, const ProcessingProperties& data);

		Duration time_before;
		Duration time_after;

		unsigned int points_before;
		unsigned int points_after;

		bool for_second_pass;
};

/**
 * @class  ProcessingBlock
 * @brief  An abstract class 
 * @author Thomas Egger
 * @date   2011-01-02
 */
class ProcessingBlock {
	public:
		virtual ~ProcessingBlock();
		
		virtual void process(const unsigned int& index, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec) = 0;

		std::string getName() const;
		const ProcessingProperties& getProperties() const;
		friend std::ostream& operator<<(std::ostream& os, const ProcessingBlock& data);

	protected:
		ProcessingBlock(const std::string& name); ///< protected constructor only to be called by children

		ProcessingProperties properties;
		std::string block_name;
};

class BlockFactory {
	public:
		static ProcessingBlock* getBlock(const std::string& blockname, const std::vector<std::string>& vec_args);

		static std::set<std::string> availableBlocks; ///<all blocks that exist
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static availableBlocks
};

} //end namespace

#endif
