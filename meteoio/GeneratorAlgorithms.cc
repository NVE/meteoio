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
#include <meteoio/meteolaws/Meteoconst.h>
#include <meteoio/meteofilters/ProcHNWDistribute.h> //for the precipitation distribution

using namespace std;

namespace mio {

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
	} else if (algoname == "CLEARSKY"){
		return new ClearSkyGenerator(vecArgs, i_algoname);
	} else if (algoname == "ALLSKY"){
		return new AllSkyGenerator(vecArgs, i_algoname);
	} else if (algoname == "UNSWORTH"){
		return new UnsworthGenerator(vecArgs, i_algoname);
	} else if (algoname == "POT_RADIATION"){
		return new PotRadGenerator(vecArgs, i_algoname);
	} else if (algoname == "HS_SWE"){
		return new HSSweGenerator(vecArgs, i_algoname);
	} else {
		throw IOException("The generator algorithm '"+algoname+"' is not implemented" , AT);
	}
}

std::string GeneratorAlgorithm::getAlgo() const {
	return algo;
}

void GeneratorAlgorithm::parse_args(const std::vector<std::string>& vecArgs)
{
	if(!vecArgs.empty()) { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
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

void SinGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	//Get the optional arguments for the algorithm: constant value to use
	if(vecArgs.size()==4) {
		const string type_str=IOUtils::strToUpper(vecArgs[0]);
		if( type_str=="YEARLY" ) type='y';
		else if( type_str=="DAILY" ) type='d';
		else
			throw InvalidArgumentException("Invalid period \""+type_str+"\" specified for the "+algo+" generator", AT);

		double min, max;
		IOUtils::convertString(min, vecArgs[1]);
		IOUtils::convertString(max, vecArgs[2]);
		amplitude = 0.5*(max-min); //the user provides min, max
		offset = min+amplitude;
		IOUtils::convertString(phase, vecArgs[3]);
	} else { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

bool SinGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if(value == IOUtils::nodata) {
		double t; //also, the minimum must occur at 0 if phase=0
		if(type=='y') {
			t = (static_cast<double>(md.date.getJulianDayNumber()) - phase*365.25) / 366.25 - .25;
		} else if(type=='d') {
			const double julian = md.date.getJulian();
			t = (julian - Optim::intPart(julian) - phase) + .25; //watch out: julian day starts at noon!
		} else {
			std::ostringstream ss;
			ss << "Invalid period \"" << type << "\" specified for the " << algo << " generator";
			throw InvalidArgumentException(ss.str(), AT);
		}

		const double w = 2.*Cst::PI;
		value = amplitude * sin(w*t) + offset;
	}

	return true; //all missing values could be filled
}

bool SinGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		generate(param, vecMeteo[ii]);
	}

	return true; //all missing values could be filled
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


void ClearSkyGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	//Get the optional arguments for the algorithm: constant value to use
	if(vecArgs.size()==1) {
		const std::string user_algo = IOUtils::strToUpper(vecArgs[0]);

		if (user_algo=="BRUTSAERT") model = BRUTSAERT;
		else if (user_algo=="DILLEY") model = DILLEY;
		else if (user_algo=="PRATA") model = PRATA;
		else
			throw InvalidArgumentException("Unknown parametrization \""+user_algo+"\" supplied for the "+algo+" generator", AT);
	} else { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

bool ClearSkyGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value==IOUtils::nodata) {
		const double TA=md(MeteoData::TA), RH=md(MeteoData::RH);
		if (TA==IOUtils::nodata || RH==IOUtils::nodata) return false;

		if (model==BRUTSAERT)
			value = Atmosphere::Brutsaert_ilwr(RH, TA);
		else if (model==DILLEY)
			value = Atmosphere::Dilley_ilwr(RH, TA);
		else if (model==PRATA)
			value = Atmosphere::Prata_ilwr(RH, TA);
	}

	return true; //all missing values could be filled
}

bool ClearSkyGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	bool all_filled = true;
	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		const bool status = generate(param, vecMeteo[ii]);
		if(status==false) all_filled=false;
	}

	return all_filled;
}


void AllSkyGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	//Get the optional arguments for the algorithm: constant value to use
	if(vecArgs.size()==1) {
		const std::string user_algo = IOUtils::strToUpper(vecArgs[0]);

		if (user_algo=="OMSTEDT") model = OMSTEDT;
		else if (user_algo=="KONZELMANN") model = KONZELMANN;
		else if (user_algo=="UNSWORTH") model = UNSWORTH;
		else if (user_algo=="CRAWFORD") model = CRAWFORD;
		else
			throw InvalidArgumentException("Unknown parametrization \""+user_algo+"\" supplied for the "+algo+" generator", AT);
	} else { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

bool AllSkyGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value==IOUtils::nodata) {
		const double TA=md(MeteoData::TA), RH=md(MeteoData::RH), cloudiness=0.5; //HACK: read cloudiness from meteoData!
		if (TA==IOUtils::nodata || RH==IOUtils::nodata || cloudiness==IOUtils::nodata) return false;

		if (model==OMSTEDT)
			value = Atmosphere::Omstedt_ilwr(RH, TA, cloudiness);
		else if (model==KONZELMANN)
			value = Atmosphere::Konzelmann_ilwr(RH, TA, cloudiness);
		else if (model==UNSWORTH)
			value = Atmosphere::Unsworth_ilwr(RH, TA, IOUtils::nodata, IOUtils::nodata, cloudiness);
		else if (model==CRAWFORD) {
			int year, month, day;
			md.date.getDate(year, month, day);
			value = Atmosphere::Crawford_ilwr(RH, TA, IOUtils::nodata, IOUtils::nodata, static_cast<unsigned char>(month), cloudiness);
		}
	}

	return true; //all missing values could be filled
}

bool AllSkyGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	bool all_filled = true;
	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		const bool status = generate(param, vecMeteo[ii]);
		if(status==false) all_filled=false;
	}

	return all_filled;
}


const double UnsworthGenerator::soil_albedo = .23; //grass
const double UnsworthGenerator::snow_albedo = .85; //snow
const double UnsworthGenerator::snow_thresh = .1; //if snow height greater than this threshold -> snow albedo

bool UnsworthGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if(value==IOUtils::nodata) {
		const double lat = md.meta.position.getLat();
		const double lon = md.meta.position.getLon();
		const double alt = md.meta.position.getAltitude();

		const double TA=md(MeteoData::TA), RH=md(MeteoData::RH), HS=md(MeteoData::HS), RSWR=md(MeteoData::RSWR);
		double ISWR=md(MeteoData::ISWR);
		if(TA==IOUtils::nodata || RH==IOUtils::nodata) return false;

		double albedo = .5;
		if(RSWR==IOUtils::nodata || ISWR==IOUtils::nodata || RSWR<=0 || ISWR<=0) {
			if(HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
				albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;

			if(ISWR==IOUtils::nodata && (RSWR!=IOUtils::nodata && HS!=IOUtils::nodata)) {
				ISWR = RSWR / albedo;
			}
		} else {
			albedo = RSWR / ISWR;
			if(albedo>=1.) albedo=0.99;
			if(albedo<=0.) albedo=0.01;
		}

		const double julian = md.date.getJulian(true);
		const double ilwr_dilley = Atmosphere::Dilley_ilwr(RH, TA);
		const double ilwr_no_iswr = ((julian - last_cloudiness_julian) < 1.)? ilwr_dilley*last_cloudiness_ratio : IOUtils::nodata;

		if(ISWR==IOUtils::nodata || ISWR<5.) {
			value = ilwr_no_iswr;
			if(value==IOUtils::nodata) return false;
		} else {
			sun.setLatLon(lat, lon, alt);
			sun.setDate(julian, 0.);

			sun.calculateRadiation(TA, RH, albedo);
			double toa, direct, diffuse;
			sun.getHorizontalRadiation(toa, direct, diffuse);
			const double ilwr_uns = Atmosphere::Unsworth_ilwr(RH, TA, ISWR, direct+diffuse);

			if(ilwr_uns==IOUtils::nodata || ilwr_uns<=0.) {
				value = ilwr_no_iswr;
				if(value==IOUtils::nodata) return false;
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
	sun.setLatLon(lat, lon, alt);

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

			double albedo = .5;
			if(RSWR==IOUtils::nodata || ISWR==IOUtils::nodata || RSWR<=0 || ISWR<=0) {
				if(HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
					albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;

				if(ISWR==IOUtils::nodata && (RSWR!=IOUtils::nodata && HS!=IOUtils::nodata)) {
					ISWR = RSWR / albedo;
				}
			} else {
				albedo = RSWR / ISWR;
				if(albedo>=1.) albedo=0.99;
				if(albedo<=0.) albedo=0.01;
			}

			const double julian = vecMeteo[ii].date.getJulian(true);
			const double ilwr_dilley = Atmosphere::Dilley_ilwr(RH, TA);
			const double ilwr_no_iswr = ((julian - last_cloudiness_julian) < 1.)? ilwr_dilley*last_cloudiness_ratio : IOUtils::nodata;

			if(ISWR==IOUtils::nodata || ISWR<5.) {
				value = ilwr_no_iswr;
				if(value==IOUtils::nodata) all_filled=false;
			} else {
				sun.setDate(julian, 0.);
				sun.calculateRadiation(TA, RH, albedo);
				double toa, direct, diffuse;
				sun.getHorizontalRadiation(toa, direct, diffuse);
				const double ilwr_uns = Atmosphere::Unsworth_ilwr(RH, TA, ISWR, direct+diffuse);

				if(ilwr_uns==IOUtils::nodata || ilwr_uns<=0.) {
					value = ilwr_no_iswr;
					if(value==IOUtils::nodata) all_filled=false;
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
const double PotRadGenerator::snow_albedo = .85; //snow
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
		const double ISWR=md(MeteoData::ISWR), RSWR=md(MeteoData::RSWR), HS=md(MeteoData::HS);
		double TA=md(MeteoData::TA), RH=md(MeteoData::RH), ILWR=md(MeteoData::ILWR);

		const double lat = md.meta.position.getLat();
		const double lon = md.meta.position.getLon();
		const double alt = md.meta.position.getAltitude();
		if(lat==IOUtils::nodata || lon==IOUtils::nodata || alt==IOUtils::nodata) return false;

		double albedo = .5;
		if(RSWR==IOUtils::nodata || ISWR==IOUtils::nodata) {
			if(HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
				albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;
		} else if(ISWR>0. && RSWR>0.) { //this could happen if the user calls this generator for a copy parameter, etc
			albedo = RSWR / ISWR;
			if(albedo>=1.) albedo=0.99;
			if(albedo<=0.) albedo=0.01;
		}

		if(TA==IOUtils::nodata || RH==IOUtils::nodata) {
			//set TA & RH so the reduced precipitable water will get an average value
			TA=274.98;
			RH=0.666;
			ILWR=IOUtils::nodata; //skip solarIndex correction
		}

		sun.setLatLon(lat, lon, alt);
		sun.setDate(md.date.getJulian(true), 0.);
		const double solarIndex = (ILWR!=IOUtils::nodata)? getSolarIndex(TA, RH, ILWR) : 1.;

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
			const double ISWR=vecMeteo[ii](MeteoData::ISWR), RSWR=vecMeteo[ii](MeteoData::RSWR), HS=vecMeteo[ii](MeteoData::HS);
			double TA=vecMeteo[ii](MeteoData::TA), RH=vecMeteo[ii](MeteoData::RH), ILWR=vecMeteo[ii](MeteoData::ILWR);

			double albedo = .5;
			if(RSWR==IOUtils::nodata || ISWR==IOUtils::nodata) {
				if(HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
					albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;
			} else if(ISWR>0. && RSWR>0.) { //this could happen if the user calls this generator for a copy parameter, etc
				albedo = RSWR / ISWR;
				if(albedo>=1.) albedo=0.99;
				if(albedo<=0.) albedo=0.01;
			}

			if(TA==IOUtils::nodata || RH==IOUtils::nodata) {
				//set TA & RH so the reduced precipitable water will get an average value
				TA=274.98;
				RH=0.666;
				ILWR=IOUtils::nodata; //skip solarIndex correction
			}

			sun.setDate(vecMeteo[ii].date.getJulian(true), 0.);
			const double solarIndex = (ILWR!=IOUtils::nodata)? getSolarIndex(TA, RH, ILWR) : 1.;

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


const bool HSSweGenerator::soft = true;
void HSSweGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	if(vecArgs.size()>0) { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

bool HSSweGenerator::generate(const size_t& /*param*/, MeteoData& /*md*/)
{//HACK: modify prototype so we can get the full vector + the index of the replacement
	return false; //all missing values could be filled
}

//when we can not guarantee HNW=0, we leave it at nodata. Therefore, it is highly recommended to
//run through a Cst=0 data generator afterward
bool HSSweGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(param!=MeteoData::HNW)
		throw InvalidArgumentException("Trying to use "+algo+" generator on " + MeteoData::getParameterName(param) + " but it can only be applied to HNW!!", AT);

	if(vecMeteo.empty()) return true;

	//Find first point that is not IOUtils::nodata
	size_t last_good = IOUtils::npos;
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		if (vecMeteo[ii](MeteoData::HS) != IOUtils::nodata){
			last_good = ii;
			break;
		}
	}

	if (last_good == IOUtils::npos) //can not find a good point to start
		return false;

	bool all_filled = (last_good>0)? false : true;
	for(size_t ii=last_good+1; ii<vecMeteo.size(); ii++) {
		const double HS_curr = vecMeteo[ii](MeteoData::HS);
		if(HS_curr==IOUtils::nodata) continue;

		const size_t start_idx = last_good+1;
		const double HS_prev = vecMeteo[last_good](MeteoData::HS);
		const double HS_delta = HS_curr - HS_prev;

		if(HS_delta>0.) {
			const double rho = newSnowDensity(vecMeteo[ii]);
			const double precip = HS_delta * rho; //in kg/m2 or mm
			ProcHNWDistribute::SmartDistributeHNW(precip, start_idx, ii, param, vecMeteo);
		} else {
			all_filled = false;
		}

		last_good=ii;
	}

	return all_filled;
}

double HSSweGenerator::newSnowDensity(const MeteoData& md) const
{ //C. Zwart, "Significance of new-snow properties for snowcover development",
//master's thesis, 2007, Institute for Marine and Atmospheric Research, University of Utrecht, 78 pp.
	const double vw = max(2., md(MeteoData::VW));
	const double rh = md(MeteoData::RH);
	const double ta = md(MeteoData::TA) - Cst::t_water_triple_pt;
	const double beta01=3.28, beta1=0.03, beta02=-0.36, beta2=-0.75, beta3=0.3;

	double arg = beta01 + beta1*ta + beta2*asin(sqrt(rh)) + beta3*log10(vw);
	if(ta>=-14.)
		arg += beta02; // += beta2*ta;

	return min( max(30., pow(10., arg)), 250. ); //limit the density to the [30, 250] kg/m3 range
}

} //namespace

