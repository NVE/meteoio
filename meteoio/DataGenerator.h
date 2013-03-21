/***********************************************************************************/
/*  Copyright 2013 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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

#ifndef __DATAGENERATOR_H__
#define __DATAGENERATOR_H__

#include <meteoio/Config.h>
#include <meteoio/MeteoData.h>
#include <meteoio/GeneratorAlgorithms.h>

#include <vector>
#include <map>

namespace mio {

/**
*
*/ //explain how to write a new generator algorithm here

/**
 * @class DataGenerator
 * @brief A class to generate meteo data from user-selected models or parametrizations.
 * This class sits in between the actual implementation of the various methods and the IOManager in
 * order to offer some high level interface. It basically reads the arguments and creates the objects for
 * the various data generators in its constructor and loop through the parameters and stations when called to fill the data.
 *
 * @ingroup meteolaws
 * @author Mathias Bavay
 * @date   2013-03-20
 */

#ifdef _POPC_
#include <paroc_base.h>
class DataGenerator : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class DataGenerator {
#endif
 	public:
		DataGenerator(const Config& i_cfg);
		DataGenerator(const DataGenerator& c) : cfg(c.cfg), mapAlgorithms(c.mapAlgorithms), generators_defined(c.generators_defined) {};

		void fillMissing(METEO_SET& vecMeteo) const;
		void fillMissing(std::vector<METEO_SET>& vecVecMeteo) const;

		DataGenerator& operator=(const DataGenerator& source);

		friend std::ostream& operator<<(std::ostream& os, const DataGenerator& mi);

	private:
		void setAlgorithms();
		size_t getAlgorithmsForParameter(const std::string& parname, std::vector<std::string>& vecAlgorithms);
		size_t getArgumentsForAlgorithm(const std::string& parname,
		                                const std::string& algorithm,
		                                std::vector<std::string>& vecArgs) const;

		const Config& cfg; ///< Reference to Config object, initialized during construction

		std::map< size_t, std::vector<GeneratorAlgorithm*> > mapAlgorithms; //per parameter interpolation algorithms
		bool generators_defined; //if true, there are some generators to run. if false, nothing to do
};

} //end namespace

#endif
