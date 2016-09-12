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
#include <meteoio/dataGenerators/GeneratorAlgorithms.h>

#include <meteoio/dataGenerators/AllSkyLWGenerator.h>
#include <meteoio/dataGenerators/AllSkySWGenerator.h>
#include <meteoio/dataGenerators/ClearSkyLWGenerator.h>
#include <meteoio/dataGenerators/ConstGenerator.h>
#include <meteoio/dataGenerators/ESOLIPGenerator.h>
#include <meteoio/dataGenerators/IswrAlbedoGenerator.h>
#include <meteoio/dataGenerators/PPHASEGenerator.h>
#include <meteoio/dataGenerators/RelHumGenerator.h>
#include <meteoio/dataGenerators/SinGenerator.h>
#include <meteoio/dataGenerators/StdPressGenerator.h>
#include <meteoio/dataGenerators/TauCLDGenerator.h>
#include <meteoio/dataGenerators/TsGenerator.h>

using namespace std;

namespace mio {

const double GeneratorAlgorithm::soil_albedo = .23; //grass
const double GeneratorAlgorithm::snow_albedo = .85; //snow
const double GeneratorAlgorithm::snow_thresh = .1; //if snow height greater than this threshold -> snow albedo

GeneratorAlgorithm* GeneratorAlgorithmFactory::getAlgorithm(const std::string& i_algoname, const std::vector<std::string>& vecArgs)
{
	std::string algoname(i_algoname);
	IOUtils::toUpper(algoname);

	if (algoname == "CST"){
		return new ConstGenerator(vecArgs, i_algoname);
	} else if (algoname == "SIN"){
		return new SinGenerator(vecArgs, i_algoname);
	} else if (algoname == "STD_PRESS"){
		return new StandardPressureGenerator(vecArgs, i_algoname);
	} else if (algoname == "RELHUM"){
		return new RhGenerator(vecArgs, i_algoname);
	} else if (algoname == "TAU_CLD"){
		return new TauCLDGenerator(vecArgs, i_algoname);
	} else if (algoname == "TS_OLWR"){
		return new TsGenerator(vecArgs, i_algoname);
	} else if (algoname == "ISWR_ALBEDO"){
		return new IswrAlbedoGenerator(vecArgs, i_algoname);
	} else if (algoname == "CLEARSKY_LW"){
		return new ClearSkyLWGenerator(vecArgs, i_algoname);
	} else if (algoname == "ALLSKY_LW"){
		return new AllSkyLWGenerator(vecArgs, i_algoname);
	} else if (algoname == "ALLSKY_SW"){
		return new AllSkySWGenerator(vecArgs, i_algoname);
	} else if (algoname == "ESOLIP"){
		return new ESOLIPGenerator(vecArgs, i_algoname);
	} else if (algoname == "PPHASE"){
		return new PPhaseGenerator(vecArgs, i_algoname);
	} else {
		throw IOException("The generator algorithm '"+algoname+"' is not implemented" , AT);
	}
}

std::string GeneratorAlgorithm::getAlgo() const {
	return algo;
}

void GeneratorAlgorithm::parse_args(const std::vector<std::string>& vecArgs)
{
	if (!vecArgs.empty()) { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

} //namespace

