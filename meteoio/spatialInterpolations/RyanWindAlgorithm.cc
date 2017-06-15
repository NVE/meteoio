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

#include <meteoio/spatialInterpolations/RyanWindAlgorithm.h>
#include <meteoio/meteoStats/libinterpol2D.h>

namespace mio {

RyanAlgorithm::RyanAlgorithm(Meteo2DInterpolator& i_mi, const std::vector< std::pair<std::string, std::string> >& vecArgs,
                                const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
                                : InterpolationAlgorithm(i_mi, vecArgs, i_algo, i_tsmanager, i_gridsmanager), vecDataVW(), vecDataDW(), scale(1e3), alpha(1.), inputIsAllZeroes(false)
{
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="SCALE") {
			parseArg(vecArgs[ii], scale);
		} else if (vecArgs[ii].first=="ALPHA") {
			parseArg(vecArgs[ii], alpha);
		}
	}
}

double RyanAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	vecMeta.clear();
	vecDataVW.clear(); vecDataDW.clear();

	tsmanager.getMeteoData(date, vecMeteo);
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		if ((vecMeteo[ii](MeteoData::VW) != IOUtils::nodata) && (vecMeteo[ii](MeteoData::DW) != IOUtils::nodata)){
			vecDataVW.push_back(vecMeteo[ii](MeteoData::VW));
			vecDataDW.push_back(vecMeteo[ii](MeteoData::DW));
			vecMeta.push_back(vecMeteo[ii].meta);
		}
	}

	nrOfMeasurments = vecMeta.size();

	if (nrOfMeasurments==0) return 0.0;

	if (Interpol2D::allZeroes(vecDataVW)) {
		inputIsAllZeroes = true;
		return 0.9;
	}

	if (nrOfMeasurments<2) return 0.6;

	return 0.9;
}

void RyanAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");

	//if all data points are zero, simply fill the grid with zeroes
	if (inputIsAllZeroes) {
		Interpol2D::constant(0., dem, grid);
		return;
	}

	if (param==MeteoData::VW) {
		Grid2DObject DW;
		simpleWindInterpolate(dem, vecDataVW, vecDataDW, grid, DW, scale, alpha);
		Interpol2D::RyanWind(dem, grid, DW);
	}
	if (param==MeteoData::DW) {
		Grid2DObject VW;
		simpleWindInterpolate(dem, vecDataVW, vecDataDW, VW, grid, scale, alpha);
		Interpol2D::RyanWind(dem, VW, grid);
	}
}

} //namespace
