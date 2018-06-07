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

#include <meteoio/dataGenerators/RelHumGenerator.h>
#include <meteoio/meteoLaws/Atmosphere.h>

namespace mio {

void HumidityGenerator::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	const std::string where( "generators::"+algo );

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="TYPE") {
			const std::string user_type( IOUtils::strToUpper(vecArgs[ii].second) );

			if (user_type=="AH") type = GEN_AH;
			else if (user_type=="RH") type = GEN_RH;
			else if (user_type=="TD") type = GEN_TD;
			else if (user_type=="QI") type = GEN_QI;
			else
				throw InvalidArgumentException("Unknown humidity parameter \""+user_type+"\" supplied for "+where, AT);
		}
	}
}

bool HumidityGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value != IOUtils::nodata) return true;

	if (type==GEN_RH) return generateRH(value, md);
	else throw IOException("Not implemented yet", AT);
}

bool HumidityGenerator::create(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if (vecMeteo.empty()) return true;

	bool all_filled = true;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		if (!generate(param, vecMeteo[ii])) all_filled=false;
	}

	return all_filled;
}

bool HumidityGenerator::generateRH(double& value, MeteoData& md)
{//HACK handle AH
	const double TA = md(MeteoData::TA);
	if (TA==IOUtils::nodata) return false;//nothing else we can do here

	//first chance to compute RH
	if (md.param_exists("TD")) {
		const double TD = md("TD");
		if (TD!=IOUtils::nodata) {
			value = Atmosphere::DewPointtoRh(TD, TA, false);
			return true;
		}
	}

	//second chance to try to compute RH
	if (md.param_exists("QI")) {
		const double QI = md("QI");
		const double altitude = md.meta.position.getAltitude();
		if (QI!=IOUtils::nodata && altitude!=IOUtils::nodata) {
			value = Atmosphere::specToRelHumidity(altitude, TA, QI);
			return true;
		}
	}

	return false;
}

} //namespace
