/***********************************************************************************/
/*  Copyright 2010 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef INTERPOLATIONALGORITHMS_H
#define INTERPOLATIONALGORITHMS_H

#include <meteoio/dataClasses/DEMObject.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/meteoStats/libinterpol1D.h>
#include <meteoio/meteoStats/libinterpol2D.h>
#include <meteoio/meteoStats/libfit1D.h>
#include <meteoio/TimeSeriesManager.h>
#include <meteoio/GridsManager.h>
#include <meteoio/Meteo2DInterpolator.h>

#include <vector>
#include <string>
#include <set>

namespace mio {

class Meteo2DInterpolator; // forward declaration, cyclic header include

/**
 * @page interpol2d Spatial interpolations
 * Using the vectors of MeteoData and StationData as filled by the IOInterface::readMeteoData call
 * as well as a grid of elevations (DEM, stored as a DEMObject), it is possible to get spatially
 * interpolated parameters.
 *
 * First, an interpolation method has to be selected for each variable which needs interpolation. Then the class computes
 * the interpolation for each 2D grid point, combining the inputs provided by the available data sources.
 * Any parameter of MeteoData can be interpolated, using the names given by \ref meteoparam. One has to keep
 * in mind that the interpolations are time-independent: each interpolation is done at a given time step and no
 * memory of (eventual) previous time steps is kept. This means that all parameters and variables that are
 * automatically calculated get recalculated anew for each time step.
 *
 * @section interpol2D_section Spatial interpolations section
 * Practically, the user
 * has to specify in his configuration file (typically io.ini), for each parameter to be interpolated, which
 * spatial interpolations algorithms should be considered, in the [Interpolations2D] section. This is provided as a space separated list of keywords
 * (one per interpolation algorithm). Please notice that some algorithms may require extra arguments.
 * Then, each algorithm will be evaluated (through the use of its rating method) and receive a grade (that might
 * depend on the number of available data, the quality of the data, etc). The algorithm that receives the higher
 * score within the user list, will be used for interpolating the selected variable at the given timestep. This means that at another
 * timestep, the same parameter might get interpolated by a different algorithm.
 * An example of such section is given below:
 * @code
 * [Interpolations2D]
 * TA::algorithms = IDW_LAPSE AVG_LAPSE
 * TA::avg_lapse = -0.008
 *
 * RH::algorithms = RH IDW_LAPSE AVG_LAPSE AVG
 *
 * PSUM::algorithms = PSUM_SNOW IDW_LAPSE AVG_LAPSE AVG CST
 * PSUM::psum_snow = avg_lapse
 * PSUM::avg_lapse = 0.0005 frac
 * PSUM::cst        = 0
 *
 * VW::algorithms = IDW_LAPSE AVG_LAPSE
 *
 * P::algorithms = STD_PRESS
 * @endcode
 *
 * @section interpol2D_keywords Available algorithms
 * The keywords defining the algorithms are the following:
 * - NONE: returns a nodata filled grid (see NoneAlgorithm)
 * - STD_PRESS: standard atmospheric pressure as a function of the elevation of each cell (see StandardPressureAlgorithm)
 * - CST: constant value in each cell (see ConstAlgorithm)
 * - AVG: average of the measurements in each cell (see AvgAlgorithm)
 * - AVG_LAPSE: constant value reprojected to the elevation of the cell (see AvgLapseRateAlgorithm)
 * - IDW: Inverse Distance Weighting averaging (see IDWAlgorithm)
 * - IDW_LAPSE: Inverse Distance Weighting averaging with reprojection to the elevation of the cell (see IDWLapseAlgorithm)
 * - LIDW_LAPSE: IDW_LAPSE restricted to a local scale (n neighbor stations, see LocalIDWLapseAlgorithm)
 * - RH: the dew point temperatures are interpolated using IDW_LAPSE, then reconverted locally to relative humidity (see RHAlgorithm)
 * - ILWR: the incoming long wave radiation is converted to emissivity and then interpolated (see ILWRAlgorithm)
 * - SWRAD: The atmospheric attenuation and splitting coefficients are evaluated and used to compute the short wave radiation with topographic shading (see SWRadInterpolation)
 * - LISTON_WIND: the wind field (VW and DW) is interpolated using IDW_LAPSE and then altered depending on the local curvature and slope (taken from the DEM, see ListonWindAlgorithm)
 * - RYAN: the wind direction is interpolated using IDW and then altered depending on the local slope (see RyanAlgorithm)
 * - WINSTRAL: the solid precipitation is redistributed by wind according to (Winstral, 2002) (see WinstralAlgorithm)
 * - PSUM_SNOW: precipitation interpolation according to (Magnusson, 2011) (see SnowPSUMInterpolation)
 * - PPHASE: precipitation phase parametrization performed at each cell (see PPHASEInterpolation)
 * - ODKRIG: ordinary kriging (see OrdinaryKrigingAlgorithm)
 * - ODKRIG_LAPSE: ordinary kriging with lapse rate (see LapseOrdinaryKrigingAlgorithm)
 * - USER: user provided grids to be read from disk (if available, see USERInterpolation)
 * - ALS_SCALING: scaling from Airborn Laser Scan data (see ALS_Interpolation)
 *
 * @section interpol2D_trends Altitudinal trends
 * Several algorithms use elevation trends, all of them relying on the same principles: the lapse rates are recomputed at each time steps
 * (see section \ref interpol2D_lapse), all stations' data are detrended with this lapse rate, the residuals are spatially interpolated
 * with the algorithm as configured by the user and finally, the values at each cell are retrended (ie the lapse rates are re-applied
 * using the cell's elevation).
 *
 * @subsection interpol2D_lapse Lapse rates
 * The altitudinal trends are currently modelled as a linear relation. The slope of this linear relation can
 * sometimes be provided by the end user (through his io.ini configuration file), otherwise it is computed from the data.
 * In order to bring slightly more robustness, if the correlation between the input data and the computed linear regression
 * is not good enought (below 0.7, as defined in Interpol2D::LinRegression), the same regression will get re-calculated
 * with one point less (cycling throught all the points). The best result (ie: highest correlation coefficient) will be
 * kept. If the final correlation coefficient is less than 0.7, a warning is displayed.
 *
 * @section interpol2D_dev_use Developer usage
 * From the developer's point of view, all that has to be done is instantiate an IOManager object and call its
 * IOManager::interpolate method.
 * @code
 * 	Config cfg("io.ini");
 * 	IOManager io(cfg);
 *
 * 	//reading the dem (necessary for several spatial interpolations algoritms)
 * 	DEMObject dem;
 * 	io.readDEM(dem);
 *
 * 	//performing spatial interpolations
 * 	Grid2DObject param;
 * 	io.interpolate(date, dem, MeteoData::TA, param);
 *
 * @endcode
 *
 * @section interpol2D_biblio Bibliography
 * The interpolation algorithms have been inspired by the following papers:
 * - <i>"A Meteorological Distribution System for High-Resolution Terrestrial Modeling (MicroMet)"</i>, Liston and Elder, Journal of Hydrometeorology <b>7</b> (2006), 217-234.
 * - <i>"Simulating wind ﬁelds and snow redistribution using terrain-based parameters to model snow accumulation and melt over a semi-arid mountain catchment"</i>, Adam Winstral and Danny Marks, Hydrological Processes <b>16</b> (2002), 3585– 3603. DOI: 10.1002/hyp.1238
 * - <i>"Quantitative evaluation of different hydrological modelling approaches in a partly glacierized Swiss watershed"</i>, Jan Magnusson, Daniel Farinotti, Tobias Jonas and Mathias Bavay, Hydrological Processes, 2010, under review.
 * - <i>"Modelling runoff from highly glacierized alpine catchments in a changing climate"</i>, Matthias Huss, Daniel Farinotti, Andreas Bauder and Martin Funk, Hydrological Processes, <b>22</b>, 3888-3902, 2008.
 * - <i>"Geostatistics for Natural Resources Evaluation"</i>, Pierre Goovaerts, Oxford University Press, Applied Geostatistics Series, 1997, 483 p., ISBN 0-19-511538-4
 * - <i>"Statistics for spatial data"</i>, Noel A. C. Cressie, John Wiley & Sons, revised edition, 1993, 900 p.
 *
 * @author Mathias Bavay
 * @date   2010-04-12
 */

/**
 * @class InterpolationAlgorithm
 * @brief A class to perform 2D spatial interpolations. For more, see \ref interpol2d
 *
 * @ingroup stats
 * @author Thomas Egger
 * @date   2010-04-01
*/
class InterpolationAlgorithm {

	public:
		InterpolationAlgorithm(Meteo2DInterpolator& i_mi,
		                       const std::vector<std::string>& i_vecArgs,
		                       const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager) :
		                      algo(i_algo), mi(i_mi), tsmanager(i_tsmanager), gridsmanager(i_gridsmanager), date(0., 0), vecArgs(i_vecArgs), vecMeteo(), vecData(),
		                      vecMeta(), info(), param(MeteoData::firstparam), nrOfMeasurments(0) {}
		virtual ~InterpolationAlgorithm() {}
		//if anything is not ok (wrong parameter for this algo, insufficient data, etc) -> return zero
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param) = 0;
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid) = 0;
		std::string getInfo() const;
		const std::string algo;

 	protected:
		size_t getData(const Date& i_date, const MeteoData::Parameters& i_param, std::vector<double>& o_vecData);
		size_t getData(const Date& i_date, const MeteoData::Parameters& i_param,
		               std::vector<double>& o_vecData, std::vector<StationData>& o_vecMeta);
		static size_t getStationAltitudes(const std::vector<StationData>& i_vecMeta, std::vector<double>& o_vecData);
		void getTrend(const std::vector<double>& vecAltitudes, const std::vector<double>& vecDat, Fit1D &trend) const;
		static void detrend(const Fit1D& trend, const std::vector<double>& vecAltitudes, std::vector<double> &vecDat, const double& min_alt=-1e4, const double& max_alt=1e4);
		static void retrend(const DEMObject& dem, const Fit1D& trend, Grid2DObject &grid, const double& min_alt=-1e4, const double& max_alt=1e4);
		void simpleWindInterpolate(const DEMObject& dem, const std::vector<double>& vecDataVW, const std::vector<double>& vecDataDW, Grid2DObject &VW, Grid2DObject &DW);

		Meteo2DInterpolator& mi;
		TimeSeriesManager& tsmanager;
		GridsManager& gridsmanager;
		Date date;
		const std::vector<std::string> vecArgs; //we must keep our own copy, it is different for each algorithm!

		std::vector<MeteoData> vecMeteo;
		std::vector<double> vecData; ///<store the measurement for the given parameter
		std::vector<StationData> vecMeta; ///<store the station data for the given parameter
		std::ostringstream info; ///<to store some extra information about the interplation process
		MeteoData::Parameters param; ///<the parameter that we will interpolate
		size_t nrOfMeasurments; ///<the available number of measurements
};

class AlgorithmFactory {
	public:
		static InterpolationAlgorithm* getAlgorithm(const std::string& i_algoname,
		                                            Meteo2DInterpolator& i_mi,
		                                            const std::vector<std::string>& i_vecArgs, TimeSeriesManager& tsm, GridsManager& gdm);
};

/**
 * @class NoneAlgorithm
 * @brief Returns a nodata filled grid
 * This allows to tolerate missing data, which can be usefull if an alternate strategy could
 * later be used to generate the data (ie. a parametrization). This algorithm will only run
 * after all others failed.
 */
class NoneAlgorithm : public InterpolationAlgorithm {
	public:
		NoneAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
};

/**
 * @class StandardPressureAlgorithm
 * @brief Standard atmospheric pressure interpolation algorithm.
 * This first fills the grid with the standard atmosphere's pressure, depending on the local elevation. Then, depending on the available data:
 *     - if there are no measured atmospheric pressure, nothing else happens;
 *     - if one station has measured local atmospheric pressure, its offset to the standard atmospheric pressure is computed and applied to
 *       the computed grid;
 *     - if multiple stations have measured local atmospheric pressure:
 *                - default: the average offset will be applied to the computed grid;
 *                - USE_RESIDUALS option: the residuals are computed at each station, spatially distributed (with IDW) and applied to the computed grid;
 * 
 * @code
 * P::algorithms = STD_PRESS
 * P::Std_Press = USE_RESIDUALS
 * @endcode
 */
class StandardPressureAlgorithm : public InterpolationAlgorithm {
	public:
		StandardPressureAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), use_residuals(false) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		bool use_residuals; ///< should we compute residuals ate each station and distribute them spatially?
};

/**
 * @class ConstAlgorithm
 * @brief Constant filling interpolation algorithm.
 * Fill the grid with a user provided constant.
 * @code
 * PSUM::algorithms = CST
 * PSUM::cst        = 0.
 * @endcode
 */
class ConstAlgorithm : public InterpolationAlgorithm {
	public:
		ConstAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), user_cst(0.) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		double user_cst;
};

/**
 * @class AvgAlgorithm
 * @brief Average filling interpolation algorithm.
 * Fill the grid with the average of the inputs for this parameter.
 * @code
 * PSUM::algorithms = AVG
 * @endcode
 */
class AvgAlgorithm : public InterpolationAlgorithm {
	public:
		AvgAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
};

/**
 * @class AvgLapseRateAlgorithm
 * @brief Average filling with elevation lapse rate interpolation algorithm.
 * The grid is filled with the average of the detrended measured values and then re-trended. Or to put it
 * differently, the following operations are performed: detrending - averaging - re-trending.
 * The lapse rate is either calculated from the data
 * (if no extra argument is provided), or given by the user-provided the optional argument <i>"avg_lapse"</i>.
 * If followed by <i>"soft"</i>, then an attempt to calculate the lapse rate from the data is made, any only if
 * unsuccessful, then user provided lapse rate is used as a fallback. If the optional user given lapse rate is
 * followed by <i>"frac"</i>, then the lapse rate is understood as a fractional lapse rate, that is a relative change
 * of the value as a function of the elevation (for example, +0.05% per meters given as 0.0005). In this case, no attempt to calculate
 * the fractional lapse from the data is made.
 * @code
 * PSUM::algorithms = AVG_LAPSE
 * PSUM::avg_lapse   = soft 0.05 frac
 * @endcode
 */
class AvgLapseRateAlgorithm : public InterpolationAlgorithm {
	public:
		AvgLapseRateAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
};

/**
 * @class IDWAlgorithm
 * @brief Inverse Distance Weighting interpolation algorithm.
 * Each cell receives the weighted average of the whole data set with weights being 1/r²
 * (r being the distance of the current cell to the contributing station) and renormalized
 * (so that the sum of the weights is equal to 1.0).
 */
class IDWAlgorithm : public InterpolationAlgorithm {
	public:
		IDWAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
};

/**
 * @class IDWLapseAlgorithm
 * @brief Inverse Distance Weighting interpolation algorithm with elevation detrending/reprojection.
 * The input data is detrended and the residuals are spatially interpolated using an Inverse Distance
 * Weighting interpolation algorithm (see IDWAlgorithm). Then, each value is reprojected to the real
 * elevation of the relative cell (re-trending). The lapse rate is either calculated from the data
 * (if no extra argument is provided), or given by the user-provided the optional argument <i>"idw_lapse"</i>.
 * If followed by <i>"soft"</i>, then an attempt to calculate the lapse rate from the data is made, any only if
 * unsuccessful or too bad (r^2<0.6), then the user provided lapse rate is used as a fallback.
 * If the optional user given lapse rate is
 * followed by <i>"frac"</i>, then the lapse rate is understood as a fractional lapse rate, that is a relative change
 * of the value as a function of the elevation (for example, +0.05% per meters given as 0.0005). In this case, no attempt to calculate
 * the fractional lapse from the data is made.
 */
class IDWLapseAlgorithm : public InterpolationAlgorithm {
	public:
		IDWLapseAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager);
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		bool lapse_rate_provided; ///< when giving a lapse rate, the requirements on the number of stations are relaxed
};


/**
 * @class LocalIDWLapseAlgorithm
 * @brief Inverse Distance Weighting interpolation algorithm with elevation detrending/reprojection.
 * The closest n stations (n being given as an extra argument of <i>"lidw_lapse"</i>) to each pixel are
 * used to compute the local lapse rate, allowing to project the contributions of these n stations to the
 * local pixel with an inverse distance weight. Beware, this method sometimes produces very sharp transitions
 * as it spatially moves from one station's area of influence to another one!
 */
class LocalIDWLapseAlgorithm : public InterpolationAlgorithm {
	public:
		LocalIDWLapseAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager);
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		size_t nrOfNeighbors;
};

/**
 * @class RHAlgorithm
 * @brief Relative humidity interpolation algorithm.
 * This is an implementation of the method described in (Liston & Elder, 2006): for each input point, the dew
 * point temperature is calculated. Then, the dew point temperatures are spatially interpolated using IDWLapseAlgorithm.
 * Finally, each local dew point temperature is converted back to a local relative humidity.
 *
 * As a side effect, the user must have defined algorithms to be used for air temperature (since this is needed for dew
 * point to RH conversion)
 */
class RHAlgorithm : public InterpolationAlgorithm {
	public:
		RHAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), vecDataTA(), vecDataRH() {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		std::vector<double> vecDataTA, vecDataRH; ///<vectors of extracted TA and RH
};

/**
 * @class ILWRAlgorithm
 * @brief Incoming Long Wave Radiation interpolation algorithm.
 * Each ILWR is converted to an emissivity (using the local air temperature), interpolated using AVG_LAPSE or IDW_LAPSE with
 * a fixed lapse rate and reconverted to ILWR.
 *
 * As a side effect, the user must have defined algorithms to be used for air temperature (since this is needed for
 * emissivity to ILWR conversion)
 */
class ILWRAlgorithm : public InterpolationAlgorithm {
	public:
		ILWRAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), vecDataEA() {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		std::vector<double> vecDataEA; ///<vectors of extracted emissivities
};

/**
 * @class ListonWindAlgorithm
 * @brief Curvature/slope influenced wind interpolation algorithm.
 * This is an implementation of the method described in G. E. Liston and K. Elder,
 * <i>"A meteorological distribution system for high-resolution terrestrial modeling (MicroMet)"</i>, Journal of Hydrometeorology, <b>7.2</b>, 2006.
 * The wind speed and direction are spatially interpolated using IDWLapseAlgorithm. Then, the wind speed and
 * direction fields are altered by wind weighting factors and wind diverting factors (respectively) calculated
 * from the local curvature and slope (as taken from the DEM, see DEMObject). The wind diverting factor is
 * actually the same as in RyanAlgorithm.
 */
class ListonWindAlgorithm : public InterpolationAlgorithm {
	public:
		ListonWindAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), vecDataVW(), vecDataDW(), inputIsAllZeroes(false) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		std::vector<double> vecDataVW, vecDataDW; ///<vectors of extracted VW and DW
		bool inputIsAllZeroes;
};

/**
 * @class RyanAlgorithm
 * @brief DEM-based wind direction interpolation algorithm.
 * This is an implementation of the method described in Ryan,
 * <i>"a mathematical model for diagnosis and prediction of surface winds in mountainous terrain"</i>,
 * 1977, journal of applied meteorology, <b>16</b>, 6.
 * The DEM is used to compute wind drection changes that are used to alter the wind direction fields.
 * @code
 * DW::algorithms    = RYAN
 * @endcode
 */
class RyanAlgorithm : public InterpolationAlgorithm {
	public:
		RyanAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), vecDataVW(), vecDataDW(), inputIsAllZeroes(false) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		std::vector<double> vecDataVW, vecDataDW; ///<vectors of extracted VW and DW
		bool inputIsAllZeroes;
};

/**
 * @class WinstralAlgorithm
 * @brief DEM-based wind-exposure interpolation algorithm.
 * This is an implementation of the method described in Winstral, Elder, & Davis,
 * <i>"Spatial snow modeling of wind-redistributed snow using terrain-based parameters"</i>, 2002,
 * Journal of Hydrometeorology, <b>3(5)</b>, 524-538.
 * The DEM is used to compute wind exposure factors that are used to alter the precipitation fields.
 * It is usually a good idea to provide a DEM that also contain the accumulated snow height in order
 * to get a progressive softening of the terrain features.
 *
 * This method must therefore first use another algorithm to generate an initial precipitation field,
 * and then modifies this field accordingly. By default, this base method is "idw_lapse" and switches to
 * "avg" if only one station can provide the precipitation at a given time step.
 *
 * Then it requires a synoptic wind direction that can be provided by different means:
 *  - without any extra argument, the stations are located in the DEM and their wind shading (or exposure)
 * is computed. If at least one station is found that is not sheltered from the wind (in every direction), it
 * provides the synoptic wind (in case of multiple stations, the vector average is used). Please note that
 * the stations that are not included in the DEM are considered to be sheltered. If no such station
 * is found, the vector average of all the available stations is used.
 *  - by providing a fixed synoptic wind bearing that is used for all time steps
 *  - by providing the station_id of the station to get the wind direction from. In this case, the base algorithm
 * for generating the initial wind field must be specified in the first position.
 *
 * @remarks Only cells with an air temperature below freezing participate in the redistribution
 * @code
 * PSUM::algorithms    = WINSTRAL
 * PSUM::winstral = idw_lapse 180
 * @endcode
 */
class WinstralAlgorithm : public InterpolationAlgorithm {
	public:
		WinstralAlgorithm(Meteo2DInterpolator& i_mi,
		                  const std::vector<std::string>& i_vecArgs,
		                  const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager);
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		void initGrid(const DEMObject& dem, Grid2DObject& grid);
		static bool windIsAvailable(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station);
		static bool isExposed(const DEMObject& dem, Coords location);
		static double getSynopticBearing(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station);
		static double getSynopticBearing(const std::vector<MeteoData>& vecMeteo);
		static double getSynopticBearing(const DEMObject& dem, const std::vector<MeteoData>& vecMeteo);

		std::string base_algo, ref_station;
		double user_synoptic_bearing;
		bool inputIsAllZeroes;
		static const double dmax;
};

class WinstralListonAlgorithm : public InterpolationAlgorithm {
	public:
		WinstralListonAlgorithm(Meteo2DInterpolator& i_mi,
		                  const std::vector<std::string>& i_vecArgs,
		                  const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager);
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		void initGrid(const DEMObject& dem, Grid2DObject& grid);
		static bool windIsAvailable(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station);
		static void getSynopticWind(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station, double& VW, double& DW);

		std::string base_algo, ref_station;
		bool inputIsAllZeroes;
		static const double dmax;
};

/**
 * @class USERInterpolation
 * @brief Reads user provided gridded data on the disk.
 * The grids are all in the GRID2DPATH directory given in the [Input] section or in one of 
 * its sub-directories that is given as the algorithm's argument (optional). By default, the file extension is assumed to 
 * be ".asc" but it is possible to provide as second argument another file extension (then it is mandatory to 
 * also provide a sub-directory argument in first position).
 * The files must be named according to the following schema: <b>{numeric date with second resolution}_{capitalized meteo parameter}.{ext}</b>, for example 20081201150000_TA.asc
 * The meteo parameters can be found in \ref meteoparam "MeteoData". Example of use:
 * @code
 * TA::algorithms = USER	# read grids from GRID2DPATH using the GRID2D plugin
 * 
 * VW::algorithms = USER	# read grids from GRID2DPATH/wind
 * VW::user       = wind
 * 
 * HNW::algorithms = USER	# read grids from GRID2DPATH/precip with the ".dat" extension
 * HNW::user       = precip .dat
 * @endcode
 *
 * If no grid exists for a given timestamp and parameter, the algorithm returns a zero rating so any other interpolation algorithm can pickup 
 * and provide a fallback. Therefore, it is not necessary to provide grids for all time steps but one can focuss on only the relevant and interesting
 * time steps.
 */
class USERInterpolation : public InterpolationAlgorithm {
	public:
		USERInterpolation(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), filename(), grid2d_path() {nrOfMeasurments=0;}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		std::string getGridFileName() const;
		std::string filename, grid2d_path;
};

/**
 * @class ALS_Interpolation
 * @brief Scale and distribute the precipitation according to Airborn Laser Scans (ALS) grids.
 * This needs two arguments: first the base method to fill the grid (for example, idw_lapse) and
 * then the name of the file (in GRID2DPATH) containing the gridded ALS data (relying on the GRID2D plugin).
 * If there are some time steps when only one station provides the necessary parameter, the base method will
 * automatically switch to "AVG". A third (optional) argument can be provided that is the air temperature
 * threshold (in K) below which such redistribution occurs (so liquid precipitation is not redistributed). 
 * 
 * @code
 * PSUM::algorithms = ALS_SCALING
 * PSUM::als_scaling = idw_lapse als_20150213.asc
 * @endcode
 */
class ALS_Interpolation : public InterpolationAlgorithm {
	public:
		ALS_Interpolation(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager);
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		void initGrid(const DEMObject& dem, Grid2DObject& grid);
		Grid2DObject ALS_scan;
		std::string filename, grid2d_path, base_algo, base_algo_user;
		double ta_thresh, als_mean; ///< the air temperature must be below a given threshold for the scaling to be applied
		bool inputIsAllZeroes;
};

/**
 * @class PPHASEInterpolation
 * @brief Precipitation phase splitting generation
 * This does not interpolate any measured precipitation phase but generates it for each point based on parametrizations, similarly to the PPHASE generator
 * (see PPhaseGenerator).
 *
 * The methods that are offered are currently the following:
 * - THRESH: a provided fixed air temperature threshold splits precipitation as either fully solid or fully liquid
 * - RANGE: two air temperature thresholds provide the lower and upper range for fully solid / fully liquid precipitation.
 *                 Within the provided range, a linear transition is assumed.
 * @code
 * PSUM::algorithms = PPHASE
 * PSUM::pphase = THRESH 274.35
 * @endcode
 */
class PPHASEInterpolation : public InterpolationAlgorithm {
	public:
		PPHASEInterpolation(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
  			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager),
  			model(THRESH), fixed_thresh(IOUtils::nodata), range_start(IOUtils::nodata), range_norm(IOUtils::nodata) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		typedef enum PARAMETRIZATION {
				THRESH,
				RANGE
			} parametrization;
		parametrization model;
		double fixed_thresh, range_start, range_norm;
};

/**
 * @class SnowPSUMInterpolation
 * @brief Precipitation distribution according to the local slope and curvature.
 * The precipitation distribution is initialized using a specified algorithm (IDW_LAPSE by default, see IDWLapseAlgorithm).
 * An optional parameter can be given to specify which algorithm has to be used for initializing the grid.
 * Please do not forget to provide the arguments of the chosen algorithm itself if necessary!
 *
 * After this initialization, the pixels whose air temperatures are below or at freezing are modified according
 * to the method described in <i>"Quantitative evaluation of different hydrological modelling approaches
 * in a partly glacierized Swiss watershed"</i>, Magnusson et Al., Hydrological Processes, <b>25</b>, 2071-2084, 2011 and
 * <i>"Modelling runoff from highly glacierized alpine catchments in a changing climate"</i>, Huss et All., Hydrological Processes, <b>22</b>, 3888-3902, 2008.
 *
 * An example using this algorithm, initializing the grid with a constant lapse rate fill using +0.05% precipitation increase per meter of elevation, is given below:
 * @code
 * PSUM::algorithms = PSUM_SNOW
 * PSUM::psum_snow = avg_lapse
 * PSUM::avg_lapse = 0.0005 frac
 * @endcode
 *
 * @author Florian Kobierska, Jan Magnusson and Mathias Bavay
 */
class SnowPSUMInterpolation : public InterpolationAlgorithm {
	public:
		SnowPSUMInterpolation(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
  			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
};

/**
 * @class SWRadInterpolation
 * @brief %Solar radiation interpolation with optional terrain shading. 
 * The splitting coefficients and an atmospheric losses factors are computed at each station that provides ISWR and spatially interpolated
 * with an Inverse Distance Weighting scheme. Then the potential radiation is computed at each pixel and scaled appropriately with the
 * atmospheric loss factor for this pixel. When applying topographic shading (default), the local splitting coefficient is used. The global, horizontal 
 * short wave radiation is then returned. To turn off the topographic shading, provide the "no_shading" argument.
 * 
 * @code
 * ISWR::algorithms = SWRad
 * ISWR::SWRad = no_shading
 * @endcode
 *
 * @note For this method to work, you also need to define spatial interpolations algorithms for TA, RH and P (a basic STD_PRESS algorithm
 * is usually enough)
 * @note This algorithm is quite time consuming (specially the topographic shading) and therefore not appropriate for very large domains.
 */
class SWRadInterpolation : public InterpolationAlgorithm {
	public:
		SWRadInterpolation(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
  			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), Sun(), vecIdx(), shading(true) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		SunObject Sun;
		std::vector<size_t> vecIdx;
		bool shading; ///<sould we also compute the shading?
		static const double soil_albedo, snow_albedo, snow_thresh;
};

/**
 * @class OrdinaryKrigingAlgorithm
 * @brief Ordinary kriging.
 * This implements ordinary krigging (see https://secure.wikimedia.org/wikipedia/en/wiki/Kriging)
 * with user-selectable variogram model (see https://secure.wikimedia.org/wikipedia/en/wiki/Variogram).
 * More details about the specific computation steps of kriging are provided in Interpol2D::ODKriging.
 *
 * The variogram is currently computed with the current data (as 1/2*(X1-X2)^2), which makes it quite
 * uninteresting... The next improvement will consist in calculating the covariances (used to build the
 * variogram) from time series (thus reflecting the time-correlation between stations).
 *
 * Please note that the variogram and krigging coefficients are re-computed fresh for each new grid (or time step).
 * The available variogram models are found in Fit1D::regression and given as optional arguments
 * (by default, LINVARIO is used). Several models can be given, the first that can fit the data will be used
 * for the current timestep:
 * @code
 * TA::algorithms = ODKRIG
 * TA::odkrig = SPHERICVARIO linvario
 * @endcode
 *
 * @author Mathias Bavay
 */
class OrdinaryKrigingAlgorithm : public InterpolationAlgorithm {
	public:
		OrdinaryKrigingAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), variogram() {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	protected:
		size_t getTimeSeries(const bool& detrend_data, std::vector< std::vector<double> > &vecVecData) const;
		void getDataForEmpiricalVariogram(std::vector<double> &distData, std::vector<double> &variData) const;
		void getDataForVariogram(std::vector<double> &distData, std::vector<double> &variData, const bool& detrend_data=false) const;
		bool computeVariogram(const bool& detrend_data=false);
		Fit1D variogram;
};


/**
 * @class LapseOrdinaryKrigingAlgorithm
 * @brief Ordinary kriging with detrending.
 * This is very similar to OrdinaryKrigingAlgorithm but performs detrending on the data.
 * @code
 * TA::algorithms = ODKRIG_LAPSE
 * TA::odkrig_lapse = SPHERICVARIO
 * @endcode
 *
 * @author Mathias Bavay
 */
class LapseOrdinaryKrigingAlgorithm : public OrdinaryKrigingAlgorithm {
	public:
		LapseOrdinaryKrigingAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: OrdinaryKrigingAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager) {}
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
};

} //end namespace mio

#endif
