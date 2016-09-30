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

#include <meteoio/dataGenerators/ClearSkySWGenerator.h>
#include <meteoio/meteoLaws/Atmosphere.h>

namespace mio {

void ClearSkySWGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	//Get the optional arguments for the algorithm: constant value to use
	if (!vecArgs.empty()) { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

bool ClearSkySWGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value == IOUtils::nodata) {
		const double ISWR=md(MeteoData::ISWR), RSWR=md(MeteoData::RSWR), HS=md(MeteoData::HS);
		double TA=md(MeteoData::TA), RH=md(MeteoData::RH);

		const double lat = md.meta.position.getLat();
		const double lon = md.meta.position.getLon();
		const double alt = md.meta.position.getAltitude();
		if (lat==IOUtils::nodata || lon==IOUtils::nodata || alt==IOUtils::nodata) return false;

		double albedo = .5;
		if (RSWR==IOUtils::nodata || ISWR==IOUtils::nodata) {
			if (HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
				albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;
		} else if (ISWR>0. && RSWR>0.) { //this could happen if the user calls this generator for a copy parameter, etc
			albedo = std::max(0.01, std::min(0.99, RSWR / ISWR));
		}

		if (TA==IOUtils::nodata || RH==IOUtils::nodata) {
			//set TA & RH so the reduced precipitable water will get an average value
			TA=274.98;
			RH=0.666;
		}

		sun.setLatLon(lat, lon, alt);
		sun.setDate(md.date.getJulian(true), 0.);

		const double P=md(MeteoData::P);
		if (P==IOUtils::nodata)
			sun.calculateRadiation(TA, RH, albedo);
		else
			sun.calculateRadiation(TA, RH, P, albedo);

		double toa, direct, diffuse;
		sun.getHorizontalRadiation(toa, direct, diffuse);
		if (param!=MeteoData::RSWR)
			value = (direct+diffuse); //ISWR
		else
			value = (direct+diffuse)*albedo; //RSWR
	}

	return true; //all missing values could be filled
}

bool ClearSkySWGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if (vecMeteo.empty()) return true;

	const double lat = vecMeteo.front().meta.position.getLat();
	const double lon = vecMeteo.front().meta.position.getLon();
	const double alt = vecMeteo.front().meta.position.getAltitude();
	if (lat==IOUtils::nodata || lon==IOUtils::nodata || alt==IOUtils::nodata) return false;
	sun.setLatLon(lat, lon, alt);

	bool all_filled = true;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		double &value = vecMeteo[ii](param);
		if (value == IOUtils::nodata) {
			const double ISWR=vecMeteo[ii](MeteoData::ISWR), RSWR=vecMeteo[ii](MeteoData::RSWR), HS=vecMeteo[ii](MeteoData::HS);
			double TA=vecMeteo[ii](MeteoData::TA), RH=vecMeteo[ii](MeteoData::RH);

			double albedo = .5;
			if (RSWR==IOUtils::nodata || ISWR==IOUtils::nodata) {
				if (HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
					albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;
			} else if (ISWR>0. && RSWR>0.) { //this could happen if the user calls this generator for a copy parameter, etc
				albedo = RSWR / ISWR;
				if (albedo>=1.) albedo=0.99;
				if (albedo<=0.) albedo=0.01;
			}

			if (TA==IOUtils::nodata || RH==IOUtils::nodata) {
				//set TA & RH so the reduced precipitable water will get an average value
				TA=274.98;
				RH=0.666;
			}

			sun.setDate(vecMeteo[ii].date.getJulian(true), 0.);

			const double P=vecMeteo[ii](MeteoData::P);
			if (P==IOUtils::nodata)
				sun.calculateRadiation(TA, RH, albedo);
			else
				sun.calculateRadiation(TA, RH, P, albedo);

			double toa, direct, diffuse;
			sun.getHorizontalRadiation(toa, direct, diffuse);
			if (param!=MeteoData::RSWR)
				value = (direct+diffuse); //ISWR
			else
				value = (direct+diffuse)*albedo; //RSWR
		}
	}

	return all_filled;
}

} //namespace
