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

#include <meteoio/spatialInterpolations/AvgLapseAlgorithm.h>
#include <meteoio/meteoStats/libinterpol1D.h>
#include <meteoio/meteoStats/libinterpol2D.h>

namespace mio {

AvgLapseRateAlgorithm::AvgLapseRateAlgorithm(Meteo2DInterpolator& i_mi,
                                         const std::vector< std::pair<std::string, std::string> >& vecArgs,
                                         const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager, const std::string& i_param)
                                         : InterpolationAlgorithm(i_mi, vecArgs, i_algo, i_tsmanager, i_gridsmanager, i_param)
{
	setTrendParams(vecArgs);
}

double AvgLapseRateAlgorithm::getQualityRating(const Date& i_date)
{
	date = i_date;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (nrOfMeasurments == 0) {
		return 0.0;
	} else if (nrOfMeasurments == 1) {
		if (user_lapse!=IOUtils::nodata)
			return 0.9; //the laspe rate is provided
		else
			return 0.0; //no lapse rate is provided and it can not be computed
	} else if (nrOfMeasurments == 2) {
		return 0.7; // as good as IDW_LAPSE
	} else if (nrOfMeasurments>2) {
		return 0.2;
	}

	return 0.2;
}

void AvgLapseRateAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	const std::vector<double> vecAltitudes( getStationAltitudes(vecMeta) );
	if (vecAltitudes.empty())
		throw IOException("Not enough data for spatially interpolating parameter " + param, AT);

	const Fit1D trend( getTrend(vecAltitudes, vecData) );
	info << trend.getInfo();
	detrend(trend, vecAltitudes, vecData);
	Interpol2D::constant(Interpol1D::arithmeticMean(vecData), dem, grid);
	retrend(dem, trend, grid);
}

} //namespace
