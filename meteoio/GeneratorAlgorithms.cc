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
#include <meteoio/GeneratorAlgorithms.h>
#include <meteoio/MathOptim.h>
#include <meteoio/meteolaws/Atmosphere.h>

using namespace std;

namespace mio {

GeneratorAlgorithm* GeneratorAlgorithmFactory::getAlgorithm(const std::string& i_algoname, const std::vector<std::string>& vecArgs)
{
	std::string algoname(i_algoname);
	IOUtils::toUpper(algoname);

	if (algoname == "CST"){
		return new ConstGenerator(vecArgs, i_algoname);
	} else if (algoname == "STD_PRESS"){
		return new StandardPressureGenerator(vecArgs, i_algoname);
	} else if (algoname == "UNSWORTH"){
		return new UnsworthGenerator(vecArgs, i_algoname);
	} else if (algoname == "POT_RADIATION"){
		return new PotRadGenerator(vecArgs, i_algoname);
	} else {
		throw IOException("The generator algorithm '"+algoname+"' is not implemented" , AT);
	}
}

std::string GeneratorAlgorithm::getAlgo() const {
	return algo;
}

////////////////////////////////////////////////////////////////////////

void ConstGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	//Get the optional arguments for the algorithm: constant value to use
	if(vecArgs.size()==1) {
		IOUtils::convertString(constant, vecArgs[0]);
	} else { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

bool ConstGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if(value == IOUtils::nodata)
		value = constant;

	return true; //all missing values could be filled
}

bool ConstGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		generate(param, vecMeteo[ii]);
	}

	return true; //all missing values could be filled
}

void StandardPressureGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	//Get the optional arguments for the algorithm: constant value to use
	if(!vecArgs.empty()) { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

bool StandardPressureGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if(value == IOUtils::nodata) {
		const double altitude = md.meta.position.getAltitude();
		if(altitude==IOUtils::nodata) return false;
		value = Atmosphere::stdAirPressure(altitude);
	}

	return true; //all missing values could be filled
}

bool StandardPressureGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	const double altitude = vecMeteo.front().meta.position.getAltitude(); //if the stations move, this has to be in the loop
	if(altitude==IOUtils::nodata) return false;

	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		double &value = vecMeteo[ii](param);
		if(value == IOUtils::nodata)
			value = Atmosphere::stdAirPressure(altitude);
	}

	return true; //all missing values could be filled
}

const double UnsworthGenerator::soil_albedo = .23; //grass
const double UnsworthGenerator::snow_albedo = .56; //white surface
const double UnsworthGenerator::snow_thresh = .1; //if snow height greater than this threshold -> snow albedo

void UnsworthGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	//Get the optional arguments for the algorithm: constant value to use
	if(!vecArgs.empty()) { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

bool UnsworthGenerator::generate(const size_t& param, MeteoData& md)
{
	const double lat = md.meta.position.getLat();
	const double lon = md.meta.position.getLon();
	const double alt = md.meta.position.getAltitude();

	double &value = md(param);
	if(value==IOUtils::nodata) {
		const double TA=md(MeteoData::TA), RH=md(MeteoData::RH), HS=md(MeteoData::HS), RSWR=md(MeteoData::RSWR);
		double ISWR=md(MeteoData::ISWR);
		if(TA==IOUtils::nodata || RH==IOUtils::nodata) return false;

		if(ISWR==IOUtils::nodata && (RSWR!=IOUtils::nodata && HS!=IOUtils::nodata)) {
			ISWR = (HS>=snow_thresh)? RSWR*snow_albedo : RSWR*soil_albedo;
		}

		const double julian = md.date.getJulian(true);
		const double ilwr_dilley = Atmosphere::Dilley_ilwr(RH, TA);
		const double ilwr_no_iswr = ((julian - last_cloudiness_julian) < 1.)? ilwr_dilley*last_cloudiness_ratio : ilwr_dilley;

		if(ISWR==IOUtils::nodata || ISWR<5.) {
			value = ilwr_no_iswr;
		} else {
			const double ilwr_uns = Atmosphere::Unsworth_ilwr(lat, lon, alt, julian, 0., RH, TA, ISWR);
			if(ilwr_uns==IOUtils::nodata || ilwr_uns<=0.) {
				value = ilwr_no_iswr;
				return true;
			}
			last_cloudiness_ratio = ilwr_uns / ilwr_dilley;
			last_cloudiness_julian = julian;
			value = ilwr_uns;
		}
	}

	return true; //all missing values could be filled
}

bool UnsworthGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	const double lat = vecMeteo.front().meta.position.getLat();
	const double lon = vecMeteo.front().meta.position.getLon();
	const double alt = vecMeteo.front().meta.position.getAltitude();

	bool all_filled = true;
	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		double &value = vecMeteo[ii](param);
		if(value==IOUtils::nodata) {
			const double TA=vecMeteo[ii](MeteoData::TA), RH=vecMeteo[ii](MeteoData::RH), HS=vecMeteo[ii](MeteoData::HS), RSWR=vecMeteo[ii](MeteoData::RSWR);
			double ISWR=vecMeteo[ii](MeteoData::ISWR);
			if(TA==IOUtils::nodata || RH==IOUtils::nodata) {
				all_filled = false;
				continue;
			}

			if(ISWR==IOUtils::nodata && (RSWR!=IOUtils::nodata && HS!=IOUtils::nodata)) {
				ISWR = (HS>=snow_thresh)? RSWR*snow_albedo : RSWR*soil_albedo;
			}

			const double julian = vecMeteo[ii].date.getJulian(true);
			const double ilwr_dilley = Atmosphere::Dilley_ilwr(RH, TA);
			const double ilwr_no_iswr = ((julian - last_cloudiness_julian) < 1.)? ilwr_dilley*last_cloudiness_ratio : ilwr_dilley;

			if(ISWR==IOUtils::nodata || ISWR<5.) {
				value = ilwr_no_iswr;
			} else {
				const double ilwr_uns = Atmosphere::Unsworth_ilwr(lat, lon, alt, julian, 0., RH, TA, ISWR);
				if(ilwr_uns==IOUtils::nodata || ilwr_uns<=0.) {
					value = ilwr_no_iswr;
					continue;
				}
				last_cloudiness_ratio = ilwr_uns / ilwr_dilley;
				last_cloudiness_julian = julian;
				value = ilwr_uns;
			}
		}
	}

	return all_filled;
}

const double PotRadGenerator::soil_albedo = .23; //grass
const double PotRadGenerator::snow_albedo = .56; //white surface
const double PotRadGenerator::snow_thresh = .1; //if snow height greater than this threshold -> snow albedo
void PotRadGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	//Get the optional arguments for the algorithm: constant value to use
	if(!vecArgs.empty()) { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

bool PotRadGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if(value == IOUtils::nodata) {
		const double TA=md(MeteoData::TA), RH=md(MeteoData::RH);
		if(TA==IOUtils::nodata || RH==IOUtils::nodata) return false;

		const double lat = md.meta.position.getLat();
		const double lon = md.meta.position.getLon();
		const double alt = md.meta.position.getAltitude();
		if(lat==IOUtils::nodata || lon==IOUtils::nodata || alt==IOUtils::nodata) return false;

		sun.setLatLon(lat, lon, alt);
		sun.setDate(md.date.getJulian(true), 0.);

		double albedo = .5;
		const double HS=md(MeteoData::HS);
		if(HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
			albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;

		double solarIndex = 1.;
		const double ILWR=md(MeteoData::ILWR);
		if(ILWR!=IOUtils::nodata)
			solarIndex = getSolarIndex(TA, RH, ILWR);

		const double P=md(MeteoData::P);
		if(P==IOUtils::nodata)
			sun.calculateRadiation(TA, RH, albedo);
		else
			sun.calculateRadiation(TA, RH, P, albedo);

		double toa, direct, diffuse;
		sun.getHorizontalRadiation(toa, direct, diffuse);
		if(param!=MeteoData::RSWR)
			value = (direct+diffuse)*solarIndex; //ISWR
		else
			value = (direct+diffuse)*albedo*solarIndex; //RSWR
	}

	return true; //all missing values could be filled
}

bool PotRadGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	const double lat = vecMeteo.front().meta.position.getLat();
	const double lon = vecMeteo.front().meta.position.getLon();
	const double alt = vecMeteo.front().meta.position.getAltitude();
	if(lat==IOUtils::nodata || lon==IOUtils::nodata || alt==IOUtils::nodata) return false;
	sun.setLatLon(lat, lon, alt);

	bool all_filled = true;
	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		double &value = vecMeteo[ii](param);
		if(value == IOUtils::nodata) {
			const double TA=vecMeteo[ii](MeteoData::TA), RH=vecMeteo[ii](MeteoData::RH);
			if(TA==IOUtils::nodata || RH==IOUtils::nodata) {
				all_filled = false;
				continue;
			}

			sun.setDate(vecMeteo[ii].date.getJulian(true), 0.);

			double albedo = .5;
			const double HS=vecMeteo[ii](MeteoData::HS);
			if(HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
				albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;

			double solarIndex = 1.;
			const double ILWR=vecMeteo[ii](MeteoData::ILWR);
			if(ILWR!=IOUtils::nodata)
				solarIndex = getSolarIndex(TA, RH, ILWR);

			const double P=vecMeteo[ii](MeteoData::P);
			if(P==IOUtils::nodata)
				sun.calculateRadiation(TA, RH, albedo);
			else
				sun.calculateRadiation(TA, RH, P, albedo);

			double toa, direct, diffuse;
			sun.getHorizontalRadiation(toa, direct, diffuse);
			if(param!=MeteoData::RSWR)
				value = (direct+diffuse)*solarIndex; //ISWR
			else
				value = (direct+diffuse)*albedo*solarIndex; //RSWR
		}
	}

	return all_filled;
}

double PotRadGenerator::getSolarIndex(const double& ta, const double& rh, const double& ilwr)
{// this is based on Kartsen cloudiness, Dilley clear sky emissivity and Unsworth ILWR
//this means that this solar index is the ratio of iswr for clear sky on a horizontal
//surface and the measured iswr
	const double epsilon_clear = Atmosphere::Dilley_emissivity(rh, ta);
	const double ilwr_clear = Atmosphere::blkBody_Radiation(1., ta);

	double cloudiness = (ilwr/ilwr_clear - epsilon_clear) / (.84 * (1.-epsilon_clear));
	if(cloudiness>1.) cloudiness=1.;
	if(cloudiness<0.) cloudiness=0.;

	const double b1 = 0.75, b2 = 3.4;
	const double karsten_Si = 1. - (b1 * pow(cloudiness, b2));
	return karsten_Si;
}

} //namespace

