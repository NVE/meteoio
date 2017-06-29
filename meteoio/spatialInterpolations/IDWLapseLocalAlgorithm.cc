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

#include <meteoio/spatialInterpolations/IDWLapseLocalAlgorithm.h>
#include <meteoio/meteoStats/libinterpol2D.h>

namespace mio {

LocalIDWLapseAlgorithm::LocalIDWLapseAlgorithm(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_param, TimeSeriesManager& i_tsm)
                      : InterpolationAlgorithm(vecArgs, i_algo, i_param, i_tsm), scale(1e3), alpha(1.), nrOfNeighbors(0)
{
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="NEIGHBORS") {
			parseArg(vecArgs[ii], nrOfNeighbors);
		} else if (vecArgs[ii].first=="SCALE") {
			parseArg(vecArgs[ii], scale);
		} else if (vecArgs[ii].first=="ALPHA") {
			parseArg(vecArgs[ii], alpha);
		}
	}

	if (nrOfNeighbors==0) throw InvalidArgumentException("Please provide the number of nearest neighbors to use for the "+algo+" algorithm", AT);
}

double LocalIDWLapseAlgorithm::getQualityRating(const Date& i_date)
{
	date = i_date;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (nrOfMeasurments == 0)
		return 0.0;

	return 0.7;
}

void LocalIDWLapseAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	if (nrOfMeasurments == 0)
		throw IOException("Interpolation FAILED for parameter " + param, AT);

	Interpol2D::LocalLapseIDW(vecData, vecMeta, dem, nrOfNeighbors, grid, scale, alpha);
	info << "using nearest " << nrOfNeighbors << " neighbors";
}

} //namespace
