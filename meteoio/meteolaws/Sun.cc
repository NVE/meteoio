/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <cmath>
#include <string>
#include <iostream> //for "fixed"
#include <iomanip> //for "setprecision"
#include <meteoio/meteolaws/Sun.h>
#include <meteoio/meteolaws/Atmosphere.h>
#include <meteoio/meteolaws/Meteoconst.h>

namespace mio {

const double SunObject::elevation_dftlThreshold = 5.; //in degrees

SunObject::SunObject(SunObject::position_algo /*alg*/) {
	julian_gmt = IOUtils::nodata;
	latitude = longitude = altitude = IOUtils::nodata;
	beam_toa = beam_direct = beam_diffuse = IOUtils::nodata;
	elevation_threshold = elevation_dftlThreshold;
}

SunObject::SunObject(const double& _latitude, const double& _longitude, const double& _altitude) {
	julian_gmt = IOUtils::nodata;
	latitude = _latitude;
	longitude = _longitude;
	altitude = _altitude;
	beam_toa = beam_direct = beam_diffuse = IOUtils::nodata;
	elevation_threshold = elevation_dftlThreshold;
}

SunObject::SunObject(const double& _latitude, const double& _longitude, const double& _altitude, const double& _julian, const double& TZ) {
	julian_gmt = _julian - TZ/24.;
	latitude = _latitude;
	longitude = _longitude;
	altitude = _altitude;
	beam_toa = beam_direct = beam_diffuse = IOUtils::nodata;
	elevation_threshold = elevation_dftlThreshold;
	update();
}

void SunObject::setDate(const double& _julian, const double& TZ) {
	position.reset();
	julian_gmt = _julian - TZ/24.;
	beam_toa = beam_direct = beam_diffuse = IOUtils::nodata;
	if(latitude!=IOUtils::nodata && longitude!=IOUtils::nodata && altitude!=IOUtils::nodata) {
		update();
	}
}

void SunObject::setLatLon(const double& _latitude, const double& _longitude, const double& _altitude) {
	position.reset();
	latitude = _latitude;
	longitude = _longitude;
	altitude = _altitude;
	beam_toa = beam_direct = beam_diffuse = IOUtils::nodata;
	if(julian_gmt!=IOUtils::nodata) {
		update();
	}
}

void SunObject::setElevationThresh(const double& _elevation_threshold) {
	elevation_threshold = _elevation_threshold;
	beam_toa = beam_direct = beam_diffuse = IOUtils::nodata;
}

double SunObject::getJulian(const double& TZ) {
	return (julian_gmt+TZ/24.);
}

void SunObject::calculateRadiation(const double& ta, const double& rh, const double& pressure, const double& mean_albedo) {
//ta in K, rh in [0,1], pressure in Pa, altitude in m
	double azimuth, elevation, eccentricity;
	position.getHorizontalCoordinates(azimuth, elevation, eccentricity);
	getBeamPotential(elevation, eccentricity, ta, rh, pressure, mean_albedo, beam_toa, beam_direct, beam_diffuse);
}

void SunObject::calculateRadiation(const double& ta, const double& rh, const double& mean_albedo) {
//ta in K, rh in [0,1], pressure in Pa, altitude in m
	double azimuth, elevation, eccentricity;
	position.getHorizontalCoordinates(azimuth, elevation, eccentricity);

	const double p = Atmosphere::stdAirPressure(altitude);
	getBeamPotential(elevation, eccentricity, ta, rh, p, mean_albedo, beam_toa, beam_direct, beam_diffuse);
}

void SunObject::update() {
	position.setAll(latitude, longitude, julian_gmt);
}

void SunObject::getBeamPotential(const double& sun_elevation, const double& Eccentricity_corr,
                                 const double& ta, const double& rh, const double& pressure, const double& mean_albedo,
                                 double& R_toa, double& R_direct, double& R_diffuse)
{
	const double olt = 0.32;      // ozone layer thickness (cm) U.S.standard = 0.34 cm
	const double w0 = 0.9;        // fraction of energy scattered to total attenuation by aerosols (Bird and Hulstrom(1981))
	const double fc = 0.84;       // fraction of forward scattering to total scattering (Bird and Hulstrom(1981))
	const double alpha = 1.3;     // wavelength exponent (Iqbal(1983) p.118)
	const double beta = 0.03;
	const double to_rad = M_PI/180.;

	if(sun_elevation<0.) { //the Sun is below the horizon, our formulas don't apply
		R_toa = 0.;
		R_direct = 0.;
		R_diffuse = 0.;
	} else {
		const double zenith = 90. - sun_elevation; //this is the TRUE zenith because the elevation is the TRUE elevation
		const double cos_zenith = cos(zenith*to_rad); //this uses true zenith angle

		// relative optical air mass Young (1994), see http://en.wikipedia.org/wiki/Airmass
		//const double mr = 1. / (cos_zenith + 0.50572 * pow( 96.07995-zenith , -1.6364 )); //pbl: this should use apparent zenith angle, and we only get true zenith angle here...
		// relative optical air mass, Young, A. T. 1994. Air mass and refraction. Applied Optics. 33:1108–1110.
		const double mr = ( 1.002432*cos_zenith*cos_zenith + 0.148386*cos_zenith + 0.0096467) /
		                  ( cos_zenith*cos_zenith*cos_zenith + 0.149864*cos_zenith*cos_zenith
		                  + 0.0102963*cos_zenith +0.000303978);

		// actual air mass: because mr is applicable for standard pressure
		// it is modified for other pressures (in Iqbal (1983), p.100)
		// pressure in Pa
		const double ma = mr * (pressure / 101325.);

		// the equations for all the transmittances of the individual atmospheric constituents
		// are from Bird and Hulstrom (1980, 1981) and can be found summarized in Iqbal (1983)
		// on the quoted pages

		// broadband transmittance by Rayleigh scattering (Iqbal (1983), p.189)
		const double taur = exp( -0.0903 * pow( ma,0.84 ) * (1 + ma - pow( ma,1.01 )) );

		// broadband transmittance by ozone (Iqbal (1983), p.189)
		const double u3 = olt * mr; // ozone relative optical path length
		const double tauoz = 1. - (0.1611 * u3 * pow( 1. + 139.48 * u3,-0.3035 ) -
		                     0.002715 * u3 * pow( 1. + 0.044  * u3 + 0.0003 * u3 * u3,-1 ));

		// broadband transmittance by uniformly mixed gases (Iqbal (1983), p.189)
		const double taug = exp( -0.0127 * pow( ma,0.26 ) );

		// saturation vapor pressure in Pa
		const double e_stern = Atmosphere::waterSaturationPressure(ta);

		// Leckner (1978); pressure and temperature correction not necessary since it is
		// included in its numerical constant (in Iqbal (1983), p.94)
		const double precw = 0.493 * rh * e_stern / ta;

		// pressure corrected relative optical path length of precipitable water (Iqbal (1983), p.176)
		const double u1 = precw * mr;

		// broadband transmittance by water vapor (in Iqbal (1983), p.189)
		const double tauw = 1 - 2.4959 * u1 * (1. / (pow( 1.0 + 79.034 * u1,0.6828 ) + 6.385 * u1));

		// broadband total transmittance by aerosols (in Iqbal (1983), pp.189-190)
		// using Angstroem's turbidity formula Angstroem (1929, 1930) for the aerosol thickness
		// in Iqbal (1983), pp.117-119
		// aerosol optical depth at wavelengths 0.38 and 0.5 micrometer
		const double ka1 = beta * pow( 0.38,-alpha );
		const double ka2 = beta * pow( 0.5 ,-alpha );

		// broadband aerosol optical depth:
		const double ka  = 0.2758 * ka1 + 0.35 * ka2;

		// total aerosol transmittance function for the two wavelengths 0.38 and 0.5 micrometer:
		const double taua = exp( -pow( ka,0.873 ) * (1. + ka - pow( ka,0.7088 )) * pow( ma,0.9108 ) );

		// broadband transmittance by aerosols due to absorption only (Iqbal (1983) p. 190)
		const double tauaa = 1. - (1. - w0) * (1. - ma + pow( ma,1.06 )) * (1. - taua);

		// broadband transmittance function due to aerosols scattering only
		// Iqbal (1983) p. 146 (Bird and Hulstrom (1981))
		const double tauas = taua / tauaa;

		// cloudless sky albedo Bird and Hulstrom (1980, 1981) (in Iqbal (1983) p. 190)
		// alb_sky = alb_rayleighscattering + alb_aerosolscattering
		const double alb_sky = 0.0685 + (1. - fc) * (1. - tauas);

		// direct normal solar irradiance in range 0.3 to 3.0 micrometer (Iqbal (1983) ,p.189)
		// 0.9751 is for the wavelength range ??
		// Bintanja (1996) (see Corripio (2002)) introduced a correction beta_z for increased
		// transmittance with altitude that is linear up to 3000 m and than fairly constant up to 5000 - 6000 m
		double beta_z;
		if( altitude < 3000. ) {
			beta_z = 2.2 * 1.e-5 * altitude;
		} else {
			beta_z = 2.2 * 1.e-5 * 3000.;
		}

		//Now calculating the radiation
		//Top of atmosphere radiation (it will always be positive, because we check for sun elevation before)
		R_toa = Cst::solcon * (1.+Eccentricity_corr) * cos_zenith;

		R_direct = 0.9751*( taur * tauoz * taug * tauw * taua + beta_z )
		           * Cst::solcon * (1.+Eccentricity_corr) * cos_zenith;

		// Diffuse radiation from the sky
		// Rayleigh-scattered diffuse radiation after the first pass through atmosphere (Iqbal (1983), p.190)
		R_diffuse = 0.79 * R_toa * tauoz * taug * tauw * tauaa
		            * 0.5 * (1. - taur ) / (1. - ma + pow( ma,1.02 ));

		// aerosol scattered diffuse radiation after the first pass through atmosphere (Iqbal (1983), p.190)
		R_diffuse += 0.79 * R_toa * tauoz * taug * tauw * tauaa
		             * fc  * (1. - tauas) / (1. - ma + pow( ma,1.02 ));

		// multiple reflected diffuse radiation between surface and sky (Iqbal (1983), p.154)
		R_diffuse += (R_diffuse + R_direct * cos_zenith) * mean_albedo * alb_sky /
		             (1. - mean_albedo * alb_sky);

		if( sun_elevation < elevation_threshold ) {
			//if the Sun is too low on the horizon, we put all the radiation as diffuse
			//the splitting calculation that might take place later on will reflect this
			R_diffuse += R_direct;
			R_direct = 0.;
		}
	}
}

void SunObject::getBeamRadiation(double& R_toa, double& R_direct, double& R_diffuse) const
{
	R_toa = beam_toa;
	R_direct = beam_direct;
	R_diffuse = beam_diffuse;
}

//diffuse remains beam_diffuse
void SunObject::getHorizontalRadiation(double& R_toa, double& R_direct, double& R_diffuse) const
{
	R_toa = position.getRadiationOnHorizontal(beam_toa);
	R_direct = position.getRadiationOnHorizontal(beam_direct);
	R_diffuse = beam_diffuse;
}

void SunObject::getSlopeRadiation(const double& slope_azi, const double& slope_elev, double& R_toa, double& R_direct, double& R_diffuse) const
{
	R_toa = position.getRadiationOnSlope(slope_azi, slope_elev, beam_toa);
	R_direct = position.getRadiationOnSlope(slope_azi, slope_elev, beam_direct);
	R_diffuse = beam_diffuse;
}

/**
 * @brief Evaluate the splitting coefficient between direct and diffuse components of the
 * incoming short wave radiation. Splitting is based on "clearness of the sky", i.e. the ratio of
 * measured incoming global radiation to top of the atmosphere radiation toa_h.g
 * Erbs et al., Iqbal p.269e
 * @param iswr_modeled modelled radiation, it should be Top Of Atmosphere Radiation (W/m²)
 * @param iswr_measured measured Incoming Short Wave Radiation on the ground (W/m²)
 * @return splitting coefficient (between 0 and 1, 1 being 100% diffuse radiation)
 */
double SunObject::getSplitting(const double& iswr_modeled, const double& iswr_measured) const
{
	const double to_rad = M_PI/180.;
	double splitting_coef;
	double azimuth, elevation;
	position.getHorizontalCoordinates(azimuth, elevation);

	if( elevation < elevation_threshold ) {
		//when the Sun is low above the horizon, Mt is getting abnormaly too large pretending 
		// this is a clear sky day when almost all the radiation should be diffuse
		// no matter how the sky is
		splitting_coef = 1.0;
	} else {
		// clear sky index (ratio global measured to top of atmosphere radiation)
		const double Mt = iswr_measured / iswr_modeled; // should be <=1.2
		const double clear_sky = 0.147;
		
		// diffuse fraction: hourly ratio of diffuse to global radiation incident on a horizontal surface
		// splitting according to a combination of Reindl et al.(1990)'s models (Mt-model and Mt&Psolar->elev-model):
		if( Mt >= 0.78 ) { // Mt in [0.78;1] -> clear day
			 splitting_coef = clear_sky;
		} else {
			if( Mt <= 0.3 ) { // Mt in [0;0.3] -> overcast
				splitting_coef = 1.02 - 0.248*Mt;
				if(splitting_coef>1.) splitting_coef=1.;
			} else {           // Mt in ]0.3;0.78[ -> cloudy
				splitting_coef = 1.4 - 1.749*Mt + 0.177*sin(elevation*to_rad);
				if(splitting_coef>0.97) splitting_coef = 0.97;
				if(splitting_coef<clear_sky) splitting_coef = clear_sky;
			}
		}
	}

	return (splitting_coef); // should be <=1.1; diff/toa should be <=0.8
}

std::ostream& operator<<(std::ostream &os, const SunObject& data)
{
	os << "<SunObject>\n";
	os << data.position;
	os << std::fixed << std::setprecision(4);
	os << "Julian (gmt)\t" << data.julian_gmt << "\n";
	os << "Lat/Long/Alt\t" << std::setw(7) << data.latitude << "° " << std::setw(7) << data.longitude << "° " << std::setprecision(0) << std::setw(4) << data.altitude << "\n";
	os << "Elev. thresh.\t" << std::setprecision(1) << data.elevation_threshold << "°\n";

	const unsigned int colw=10;
	os << std::setw(colw) << "" << std::fixed << std::setw(colw) << std::setprecision(1) << "toa";
	os << std::fixed << std::setw(colw) << std::setprecision(1) << "direct";
	os << std::fixed << std::setw(colw) << std::setprecision(1) << "diffuse";
	os << std::fixed << std::setw(colw) << std::setprecision(1) << "sum\n";

	os << std::setw(colw) << "Beam" << std::fixed << std::setw(colw) << std::setprecision(1) << data.beam_toa;
	os << std::fixed << std::setw(colw) << std::setprecision(1) << data.beam_direct;
	os << std::fixed << std::setw(colw) << std::setprecision(1) << data.beam_diffuse;
	os << std::fixed << std::setw(colw) << std::setprecision(1) << data.beam_direct+data.beam_diffuse << "\n";

	double R_toa, R_direct, R_diffuse;
	data.getHorizontalRadiation(R_toa, R_direct, R_diffuse);
	os << std::setw(colw) << "Horizontal" << std::fixed << std::setw(colw) << std::setprecision(1) << R_toa;
	os << std::fixed << std::setw(colw) << std::setprecision(1) << R_direct;
	os << std::fixed << std::setw(colw) << std::setprecision(1) << R_diffuse;
	os << std::fixed << std::setw(colw) << std::setprecision(1) << R_direct+R_diffuse << "\n";

	os << "</SunObject>\n";
	return os;
}

} //end namespace
