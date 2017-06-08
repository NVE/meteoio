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
#include <meteoio/dataGenerators/ClearSkySWGenerator.h>
#include <meteoio/dataGenerators/ConstGenerator.h>
#include <meteoio/dataGenerators/ESOLIPGenerator.h>
#include <meteoio/dataGenerators/IswrAlbedoGenerator.h>
#include <meteoio/dataGenerators/PPHASEGenerator.h>
#include <meteoio/dataGenerators/PrecUnsplit.h>
#include <meteoio/dataGenerators/RelHumGenerator.h>
#include <meteoio/dataGenerators/SinGenerator.h>
#include <meteoio/dataGenerators/StdPressGenerator.h>
#include <meteoio/dataGenerators/TauCLDGenerator.h>
#include <meteoio/dataGenerators/TsGenerator.h>

namespace mio {

/**
 * @page generators Data generators
 * New data can be generated based on some parametrizations at two very different stages:
 *    + in raw data editing, when calling a data creator;
 *    + when the requested data could not be provided as last resort as data generator.
 *
 * In the first case, the goal is to create new parameters fully based on parametrizations. In such a case, the "generator" is called
 * a "creator" and behaves the same way as a generator, except that it creates an additional parameter. It is declared as
 * {new_parameter}::%create = {data generators} in the [Input] section (see \ref data_creation "data creation" in the
 * \ref data_manipulations "Raw data editing" section).
 *
 * The second case takes place once the data has been read, filtered and resampled, if some data points are still missing.
 * These are either a few isolated periods (a sensor was not functioning) that are too large for performing
 * a statistical temporal interpolation or that a meteorological parameter was not even measured. In such a case,
 * we generate data, generally relying on some parametrization using other meteorological parameters. In a few
 * cases, even fully arbitrary data might be helpful (replacing missing value by a given constant so a model can
 * run over the data gap).
 *
 * @note it is generally not advised to use data generators in combination with spatial interpolations as this would
 * potentially mix measured and generated values in the resulting grid. It is therefore advised to turn the data generators
 * off and let the spatial interpolations algorithms adjust to the amount of measured data.
 * @note it is also possible to make a copy of a given parameter under a different name. This is explained in section
 * \ref data_manipulations "Raw data editing".
 *
 * @section generators_section Data generators section
 * The data generators are defined per meteorological parameter. They are applied to all stations
 * (if using multiple meteorological stations). If multiple dat generators are specified for each parameter,
 * they would be used in the order of declaration, meaning that only the data points that could not be
 * generated by the first generator would be tentatively generated by the second generator, etc. Please also keep
 * in mind that at this stage, all data <b>must be in SI</b> units!
 * @code
 * [Generators]
 * RH::generators = CST
 * RH::Cst        = .7
 *
 * P::generators  = STD_PRESS
 *
 * ILWR::generators = AllSky_LW ClearSky_LW
 *
 * TAU_CLD::create  = CST
 * TA_CLD::Cst      = 0.5
 * @endcode
 *
 * @section generators_keywords Available generators
 * The keywords defining the algorithms are the following:
 * - STD_PRESS: standard atmospheric pressure as a function of the elevation of each station (see StandardPressureGenerator)
 * - RELHUM: relative humidity from other humidity measurements (see RhGenerator)
 * - TS_OLWR: surface temperature from Outgoing Long Wave Radiation (see TsGenerator)
 * - ISWR_ALBEDO: ISWR from RSWR or RSWR from ISWR with either a snow or a soil albedo, depending on HS (see IswrAlbedoGenerator)
 * - CST: constant value as provided in argument (see ConstGenerator)
 * - SIN: sinusoidal variation (see SinGenerator)
 * - CLEARSKY_LW: use a clear sky model to generate ILWR from TA, RH (see ClearSkyLWGenerator)
 * - ALLSKY_LW: use an all sky model to generate ILWR from TA, RH and cloudiness (see AllSkyLWGenerator)
 * - CLEARSKY_SW: use a clear sky model to generate ISWR from TA, RH (see ClearSkySWGenerator)
 * - ALLSKY_SW: generate the incoming short wave radiation from the potential radiation, corrected for cloudiness if possible (see AllSkySWGenerator)
 * - TAU_CLD: generate the atmospheric transmissivity based on cloud cover fraction (see TauCLDGenerator)
 * - ESOLIP: generate precipitation from snow height changes (see ESOLIPGenerator)
 * - PPHASE: generate precipitation phase with a user-selected method (see PPhaseGenerator)
 * - PrecUnsplit: generate the precipitation amount and/or phase from split precipitation (see PrecUnsplit)
 *
 * @section generators_biblio Bibliography
 * The data generators have been inspired by the following papers:
 * - Brutsaert -- <i>"On a Derivable Formula for Long-Wave Radiation From Clear Skies"</i>, Journal of Water Resources
 * Research, <b>11</b>, No. 5, October 1975, pp 742-744.
 * - Prata -- <i>"A new long-wave formula for estimating downward clear-sky radiation at the surface"</i>, Q. J. R. Meteorolo. Soc., <b>122</b>, 1996, pp 1127-1151.
 * - Dilley and O'Brien -- <i>"Estimating downward clear sky long-wave irradiance at the surface from screen temperature and precipitable water"</i>, Q. J. R. Meteorolo. Soc., Vol. 124, 1998, doi:10.1002/qj.49712454903
 * - Clark & Allen -- <i>"The estimation of atmospheric radiation for clear and cloudy skies"</i>, Proceedings of the second national passive solar conference, <b>2</b>, 1978, p 676.
 * - Tang et al. -- <i>"Estimates of clear night sky emissivity in the Negev Highlands, Israel"</i>, Energy Conversion and Management, <b>45.11</b>, 2004, pp 1831-1843.
 * - Idso -- <i>"A set of equations for full spectrum and 8 to 14 um and 10.5 to 12.5 um thermal radiation from cloudless skies"</i>, Water Resources Research, <b>17</b>, 1981, pp 295-304.
 * - Kasten and Czeplak -- <i>"Solar and terrestrial radiation dependent on the amount and type of cloud"</i>, 1980, Solar energy, 24.2, pp 177-189
 * - Omstedt -- <i>"A coupled one-dimensional sea ice-ocean model applied to a semi-enclosed basin"</i>, Tellus, <b>42 A</b>, 568-582, 1990, DOI:10.1034/j.1600-0870.1990.t01-3-00007.
 * - Konzelmann et al. -- <i>"Parameterization of global and longwave incoming radiation for the Greenland Ice Sheet."</i> Global and Planetary change <b>9.1</b> (1994): 143-164.
 * - Crawford and Duchon -- <i>"An Improved Parametrization for Estimating Effective Atmospheric Emissivity for Use in Calculating Daytime
 * Downwelling Longwave Radiation"</i>, Journal of Applied Meteorology, <b>38</b>, 1999, pp 474-480
 * - Unsworth and Monteith -- <i>"Long-wave radiation at the ground"</i>, Q. J. R. Meteorolo. Soc., Vol. 101, 1975, pp 13-24
 * - Meeus -- <i>"Astronomical Algorithms"</i>, second edition, 1998, Willmann-Bell, Inc., Richmond, VA, USA
 * - Mair et al. -- <i>" ESOLIP–estimate of solid and liquid precipitation at sub-daily time resolution by combining snow height
 * and rain gauge measurements"</i>, Hydrology and Earth System Sciences Discussions, <b>10(7)</b>, 8683-8714, 2013.
 *
 *
 * @author Mathias Bavay
 * @date   2013-03-20
 */

const double GeneratorAlgorithm::soil_albedo = .23; //grass
const double GeneratorAlgorithm::snow_albedo = .85; //snow
const double GeneratorAlgorithm::snow_thresh = .1; //if snow height greater than this threshold -> snow albedo

GeneratorAlgorithm* GeneratorAlgorithmFactory::getAlgorithm(const Config& /*cfg*/, const std::string& i_algoname, const std::vector< std::pair<std::string, std::string> >& vecArgs)
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
	} else if (algoname == "CLEARSKY_SW"){
		return new ClearSkySWGenerator(vecArgs, i_algoname);
	} else if (algoname == "ESOLIP"){
		return new ESOLIPGenerator(vecArgs, i_algoname);
	} else if (algoname == "PPHASE"){
		return new PPhaseGenerator(vecArgs, i_algoname);
	} else if (algoname == "PRECUNSPLIT"){
		return new PrecUnsplit(vecArgs, i_algoname);
	} else {
		throw IOException("The generator algorithm '"+algoname+"' is not implemented" , AT);
	}
}

void GeneratorAlgorithm::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	if (!vecArgs.empty()) { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

} //namespace

