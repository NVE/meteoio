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

#include <meteoio/dataGenerators/TauCLDGenerator.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <algorithm>

namespace mio {

TauCLDGenerator::TauCLDGenerator(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& i_algo, const std::string& i_section, const double& TZ, const bool& parse_args)
                              : GeneratorAlgorithm(vecArgs, i_algo, i_section, TZ), last_cloudiness(), cloudiness_model(KASTEN), use_rswr(false), use_rad_threshold(false)
{
	if (!parse_args) return;
	
	const std::string where( section+"::"+algo );
	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="TYPE") {
			const std::string user_algo( IOUtils::strToUpper(vecArgs[ii].second) );
			
			if (user_algo=="LHOMME") cloudiness_model = CLF_LHOMME;
			else if (user_algo=="KASTEN") cloudiness_model = KASTEN;
			else if (user_algo=="CRAWFORD") cloudiness_model = CLF_CRAWFORD;
			else
				throw InvalidArgumentException("Unknown parametrization \""+user_algo+"\" supplied for "+where, AT);
		}
		if (vecArgs[ii].first=="USE_RSWR") {
			IOUtils::parseArg(vecArgs[ii], where, use_rswr);
		}
		if (vecArgs[ii].first=="USE_RAD_THRESHOLD") {
			IOUtils::parseArg(vecArgs[ii], where, use_rad_threshold);
		}
	}
}

/**
 * @brief Compute the clearness index from an atmospheric cloudiness value
 * @details
 * This is a convenience method that helps process the same way various types of inputs: if a cloudiness
 * is provided (which is quite rare), it can be converted to a clearness index (ie the ratio of the incoming 
 * short wave radiation over the ground potential radiation, projected on the horizontal) and then processed 
 * the same way as more traditional measurements (ie only ISWR provided) where it will be re-converted
 * to a cloudiness (thus falling abck to the same cloudiness as originally provided).

 * @param[in] cloudiness cloudiness (between 0 and 1)
 * @return clearness index (between 0 and 1)
 */
double TauCLDGenerator::getClearness(const double& cloudiness) const
{
	if (cloudiness_model==CLF_LHOMME) {
		return Atmosphere::Lhomme_cloudiness( cloudiness/8. );
	} else if (cloudiness_model==KASTEN) {
		return Atmosphere::Kasten_cloudiness( cloudiness/8. );
	} else if (cloudiness_model==CLF_CRAWFORD) {
		return Atmosphere::Lhomme_cloudiness( cloudiness/8. );
	} else
		return IOUtils::nodata; //this should never happen
}

/**
 * @brief Compute the atmospheric cloudiness from the available measurements
 * @details
 * The clearness index (ie the ratio of the incoming short wave radiation over the ground potential radiation, projected on the horizontal) 
 * is computed and used to evaluate the cloudiness, based on the chosen parametrization.
 * @param[in] clf_model cloudiness parametrization
 * @param[in] md MeteoData
 * @param[in] i_use_rswr if set to true, in case of no iswr measurements, a ground albedo is assumed and used to compute iswr. Based on HS, this albedo can either
 * be a soil ro a snow albedo
 * @param[in] i_use_rad_threshold use a radiation threshold to force resampling of cloudiness over prediods of low ISWR?
 * @param sun For better efficiency, the SunObject for this location (so it can be cached)
 * @param[out] is_night set to TRUE if it is night time
 * @return cloudiness (between 0 and 1)
 */
double TauCLDGenerator::getCloudiness(const MeteoData& md, SunObject& sun, bool &is_night) const
{
	//we know that TA and RH are available, otherwise we would not get called
	const double TA=md(MeteoData::TA), RH=md(MeteoData::RH), HS=md(MeteoData::HS), RSWR=md(MeteoData::RSWR);
	double ISWR=md(MeteoData::ISWR);

	//at sunrise or sunset, we might get very wrong results -> return nodata in order to use interpolation instead
	//obviously, when it is really night neither can we compute anything here...
	is_night = (sun.position.getSolarElevation() <= 5.);
	if (is_night) return IOUtils::nodata;
	
	double albedo = .5;
	if (RSWR!=IOUtils::nodata && ISWR!=IOUtils::nodata && RSWR>Atmosphere::day_iswr_thresh && ISWR>Atmosphere::day_iswr_thresh) {
		albedo = std::min( 0.99, std::max(0.01, RSWR / ISWR) );
	} else { //so some measurements are missing
		if (HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
			albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;

		if (ISWR==IOUtils::nodata) { //ISWR is missing, trying to compute it
			if (!use_rswr) return IOUtils::nodata;
			if (RSWR!=IOUtils::nodata && HS!=IOUtils::nodata)
				ISWR = RSWR / albedo;
			else
				return IOUtils::nodata; //no way to get ISWR, aborting
		}
	}

	sun.calculateRadiation(TA, RH, albedo);
	double toa, direct, diffuse;
	sun.getHorizontalRadiation(toa, direct, diffuse);
	const double iswr_clear_sky = direct+diffuse;
	
	if (use_rad_threshold) {
		if (ISWR<Atmosphere::day_iswr_thresh || iswr_clear_sky<Atmosphere::day_iswr_thresh || iswr_clear_sky<ISWR) {
			is_night = true;
			return IOUtils::nodata;
		}
	}

	const double clearness_index = std::min(ISWR/iswr_clear_sky, 1.);
	if (cloudiness_model==CLF_LHOMME) {
		const double clf = Atmosphere::Lhomme_cloudiness( clearness_index );
		if (clf<0. || clf>1.) return IOUtils::nodata;
		return clf;
	} else if (cloudiness_model==KASTEN) {
		const double clf = Atmosphere::Kasten_cloudiness( clearness_index );
		if (clf<0. || clf>1.) return IOUtils::nodata;
		return clf;
	} else if (cloudiness_model==CLF_CRAWFORD) {
		const double clf = Atmosphere::Lhomme_cloudiness( clearness_index );
		if (clf<0. || clf>1.) return IOUtils::nodata;
		return clf;
	} else
		return IOUtils::nodata; //this should never happen
}

bool TauCLDGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value == IOUtils::nodata) {
		double cld = (md.param_exists("CLD"))? md("CLD") : IOUtils::nodata;
		if (cld!=IOUtils::nodata) {
			if (cld==9) cld=8.; //Synop sky obstructed from view -> fully cloudy
			if (cld>8. || cld<0.) throw InvalidArgumentException("Cloud cover CLD should be between 0 and 8!", AT);
			value = getClearness( cld );
			return true;
		}

		const double TA=md(MeteoData::TA), RH=md(MeteoData::RH);
		if (TA==IOUtils::nodata || RH==IOUtils::nodata) return false;

		const std::string station_hash = md.meta.stationID + ":" + md.meta.stationName;
		const double julian_gmt = md.date.getJulian(true);
		bool cloudiness_from_cache = false;

		const double lat = md.meta.position.getLat();
		const double lon = md.meta.position.getLon();
		const double alt = md.meta.position.getAltitude();
		if (lat==IOUtils::nodata || lon==IOUtils::nodata || alt==IOUtils::nodata) return false;
		SunObject sun;
		sun.setLatLon(lat, lon, alt);
		sun.setDate(julian_gmt, 0.);

		bool is_night;
		double cloudiness = TauCLDGenerator::getCloudiness(md, sun, is_night);
		if (cloudiness==IOUtils::nodata && !is_night) return false;

		if (is_night) { //interpolate the cloudiness over the night
			const std::map< std::string, std::pair<double, double> >::const_iterator it = last_cloudiness.find(station_hash);
			if (it==last_cloudiness.end()) return false;

			cloudiness_from_cache = true;
			const double last_cloudiness_julian = it->second.first;
			const double last_cloudiness_value = it->second.second;
			if ((julian_gmt - last_cloudiness_julian) < 1.) cloudiness = last_cloudiness_value;
			else return false;
		}

		//save the last valid cloudiness
		if (!cloudiness_from_cache)
			last_cloudiness[station_hash] = std::pair<double,double>( julian_gmt, cloudiness );

		value = 1.- cloudiness;
	}

	return true; //all missing values could be filled
}

bool TauCLDGenerator::create(const size_t& param, const size_t& ii_min, const size_t& ii_max, std::vector<MeteoData>& vecMeteo)
{
	if (vecMeteo.empty()) return true;

	bool status = true;
	for (size_t ii=ii_min; ii<ii_max; ii++) {
		if (!generate(param, vecMeteo[ii]))
			status = false;
	}

	return status;
}

} //namespace
