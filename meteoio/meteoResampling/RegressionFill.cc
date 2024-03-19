// SPDX-License-Identifier: LGPL-3.0-or-later
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

#include <meteoio/meteoResampling/RegressionFill.h>

#include <sstream>

namespace mio {

RegressionFill::RegressionFill(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs)
             : ResamplingAlgorithms(i_algoname, i_parname, dflt_window_size, vecArgs)
{
	/*implement here the arguments parsing
	const std::string where( "Interpolations1D::"+i_parname+"::"+i_algoname );
	if (!vecArgs.empty()) //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments for \""+where+"\"", AT);*/
}

std::string RegressionFill::toString() const
{
	//this should help when debugging, so output relevant parameters for your algorithm
	std::ostringstream ss;
	ss << std::right << std::setw(10) << parname << "::"  << std::left << std::setw(15) << algo << "[ ]";
	return ss.str();
}

double linear(double julian_date, const std::vector<double>& coefficients)
{
	return coefficients[0] + coefficients[1] * julian_date;
}

void RegressionFill::resample(const std::string& /*stationHash*/, const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
                            const std::vector<MeteoData>& vecM, MeteoData& md) {
								throw IOException("The Regression Fill needs additional stations to work properly", AT);
							}

void RegressionFill::resample(const std::string& /*stationHash*/, const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
                            const std::vector<MeteoData>& vecM, MeteoData& md, const std::vector<METEO_SET>& additional_stations)
{
	if (additional_stations.empty()) throw IOException("The Regression Fill needs additional stations to work properly. Make sure to use the correct Station IDs", AT);
	if (index >= vecM.size()) throw IOException("The index of the element to be resampled is out of bounds", AT);

	if (position == ResamplingAlgorithms::exact_match) {
		const double value = vecM[index](paramindex);
		if (value != IOUtils::nodata) {
			md(paramindex) = value;
			return;
		}
	}

	if (regression_coefficients.find(index) != regression_coefficients.end()) {
		md(paramindex) = linear(md.date.getJulian(true), regression_coefficients[index]); // TODO: do i need to convert to gmt?
		return;
	}

	// get the regression data before the missing value
	std::vector<double> x, y;

	return;
}

} //namespace
