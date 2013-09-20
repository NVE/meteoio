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
#ifndef __INTERPOLATIONALGORITHMS_H__
#define __INTERPOLATIONALGORITHMS_H__

#include <meteoio/DEMObject.h>
#include <meteoio/MeteoData.h>
#include <meteoio/meteostats/libinterpol1D.h>
#include <meteoio/meteostats/libinterpol2D.h>
#include <meteoio/meteostats/libfit1D.h>

#include <vector>
#include <string>
#include <set>

namespace mio {

class IOManager;
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
 * score within the user list, will be used for interpolating the selected variable. An example of such section is given below:
 * @code
 * [Interpolations2D]
 * TA::algorithms = IDW_LAPSE CST_LAPSE
 * TA::cst_lapse = -0.008
 *
 * RH::algorithms = RH IDW_LAPSE CST_LAPSE CST
 *
 * HNW::algorithms = HNW_SNOW IDW_LAPSE CST_LAPSE CST
 * HNW::hnw_snow = cst_lapse
 * HNW::cst_lapse = 0.0005 frac
 *
 * VW::algorithms = IDW_LAPSE CST_LAPSE
 *
 * P::algorithms = STD_PRESS
 * @endcode
 *
 * @section interpol2D_keywords Available algorithms
 * The keywords defining the algorithms are the following:
 * - STD_PRESS: standard atmospheric pressure as a function of the elevation of each cell (see StandardPressureAlgorithm)
 * - CST: constant value in each cell (see ConstAlgorithm)
 * - CST_LAPSE: constant value reprojected to the elevation of the cell (see ConstLapseRateAlgorithm)
 * - IDW: Inverse Distance Weighting averaging (see IDWAlgorithm)
 * - IDW_LAPSE: Inverse Distance Weighting averaging with reprojection to the elevation of the cell (see IDWLapseAlgorithm)
 * - RH: the dew point temperatures are interpolated using IDW_LAPSE, then reconverted locally to relative humidity (see RHAlgorithm)
 * - ILWR: the incoming long wave radiation is converted to emissivity and then interpolated (see ILWRAlgorithm)
 * - WIND_CURV: the wind field (VW and DW) is interpolated using IDW_LAPSE and then altered depending on the local curvature and slope (taken from the DEM, see SimpleWindInterpolationAlgorithm)
 * - HNW_SNOW: precipitation interpolation according to (Magnusson, 2011) (see SnowHNWInterpolation)
 * - ODKRIG: ordinary kriging (see OrdinaryKrigingAlgorithm)
 * - ODKRIG_LAPSE: ordinary kriging with lapse rate (see LapseOrdinaryKrigingAlgorithm)
 * - USER: user provided grids to be read from disk (if available, see USERInterpolation)
 * <!-- - LIDW_LAPSE: IDW_LAPSE restricted to a local scale (n neighbor stations, see LocalIDWLapseAlgorithm) -->
 *
 * @section interpol2D_lapse Lapse rates
 * Several algorithms use elevation trends, currently modelled as a linear relation. The slope of this linear relation can
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
 *	//performing spatial interpolations
 * 	Grid2DObject param;
 *	io.interpolate(date, dem, MeteoData::TA, param);
 *
 * @endcode
 *
 * @section interpol2D_biblio Bibliography
 * The interpolation algorithms have been inspired by the following papers:
 * - <i>"A Meteorological Distribution System for High-Resolution Terrestrial Modeling (MicroMet)"</i>, Liston and Elder, Journal of Hydrometeorology <b>7</b> (2006), 217-234.
 * - <i>"Simulating wind ﬁelds and snow redistribution using terrain-based parameters to model snow accumulation and melt over a semi-arid mountain catchment"</i>, Adam Winstral and Danny Marks, Hydrological Processes <b>16</b> (2002), 3585– 3603. DOI: 10.1002/hyp.1238 [NOT YET IMPLEMENTED]
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
		                       const std::string& i_algo, IOManager& iom) :
		                      algo(i_algo), mi(i_mi), date(0.), vecArgs(i_vecArgs), vecMeteo(), vecData(),
		                      vecMeta(), info(), param(MeteoData::firstparam), nrOfMeasurments(0), iomanager(iom) {};
		virtual ~InterpolationAlgorithm() {};
		//if anything is not ok (wrong parameter for this algo, insufficient data, etc) -> return zero
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param) = 0;
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid) = 0;
		std::string getInfo() const;
		const std::string algo;

 	protected:
		size_t getData(const Date& i_date, const MeteoData::Parameters& i_param, std::vector<double>& o_vecData);
		size_t getData(const Date& i_date, const MeteoData::Parameters& i_param,
		               std::vector<double>& o_vecData, std::vector<StationData>& o_vecMeta);
		size_t getStationAltitudes(const std::vector<StationData>& i_vecMeta, std::vector<double>& o_vecData) const;
		void getTrend(const std::vector<double>& vecAltitudes, const std::vector<double>& vecDat, Fit1D &trend) const;
		void detrend(const Fit1D& trend, const std::vector<double>& vecAltitudes, std::vector<double> &vecDat) const;
		void retrend(const DEMObject& dem, const Fit1D& trend, Grid2DObject &grid) const;

		Meteo2DInterpolator& mi;
		Date date;
		const std::vector<std::string> vecArgs; //we must keep our own copy, it is different for each algorithm!

		std::vector<MeteoData> vecMeteo;
		std::vector<double> vecData; ///<store the measurement for the given parameter
		std::vector<StationData> vecMeta; ///<store the station data for the given parameter
		std::ostringstream info; ///<to store some extra information about the interplation process
		MeteoData::Parameters param; ///<the parameter that we will interpolate
		size_t nrOfMeasurments; ///<the available number of measurements
		IOManager& iomanager;
};

class AlgorithmFactory {
	public:
		static InterpolationAlgorithm* getAlgorithm(const std::string& i_algoname,
		                                            Meteo2DInterpolator& i_mi,
		                                            const std::vector<std::string>& i_vecArgs, IOManager& iom);
};

/**
 * @class ConstAlgorithm
 * @brief Constant filling interpolation algorithm.
 * Fill the grid with the average of the inputs for this parameter.
 */
class ConstAlgorithm : public InterpolationAlgorithm {
	public:
		ConstAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, IOManager& iom)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, iom) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
};

/**
 * @class StandardPressureAlgorithm
 * @brief Standard atmospheric pressure interpolation algorithm.
 * Fill the grid with the standard atmosphere's pressure, depending on the local elevation.
 */
class StandardPressureAlgorithm : public InterpolationAlgorithm {
	public:
		StandardPressureAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, IOManager& iom)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, iom) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
};

/**
 * @class ConstLapseRateAlgorithm
 * @brief Constant filling with elevation lapse rate interpolation algorithm.
 * Assuming that average values occured at the average of the elevations, the grid is filled with average values
 * reprojected to real grid elevation according to a lapse rate. The lapse rate is either calculated from the data
 * (if no extra argument is provided), or given by the user-provided the optional argument <i>"cst_lapse"</i>.
 * If followed by <i>"soft"</i>, then an attempt to calculate the lapse rate from the data is made, any only if
 * unsuccessful, then user provided lapse rate is used as a fallback. If the optional user given lapse rate is
 * followed by <i>"frac"</i>, then the lapse rate is understood as a fractional lapse rate, that is a relative change
 * of the value as a function of the elevation (for example, +0.05% per meters given as 0.0005). In this case, no attempt to calculate
 * the fractional lapse from the data is made.
 */
class ConstLapseRateAlgorithm : public InterpolationAlgorithm {
	public:
		ConstLapseRateAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, IOManager& iom)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, iom) {}
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
					const std::string& i_algo, IOManager& iom)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, iom) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
};

/**
 * @class IDWLapseAlgorithm
 * @brief Inverse Distance Weighting interpolation algorithm with elevation detrending/reprojection.
 * The input data is projected to a reference elevation and spatially interpolated using an Inverse Distance
 * Weighting interpolation algorithm (see IDWAlgorithm). Then, each value is reprojected to the real
 * elevation of the relative cell. The lapse rate is either calculated from the data
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
					const std::string& i_algo, IOManager& iom)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, iom) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
};


/**
 * @class LocalIDWLapseAlgorithm
 * @brief Inverse Distance Weighting interpolation algorithm with elevation detrending/reprojection.
 * The closest n stations (n being given as an extra argument of <i>"idw_lapse"</i>) to each pixel are
 * used to compute the local lapse rate, allowing to project the contributions of these n stations to the
 * local pixel with an inverse distance weight. We also assume a two segments regression for altitude detrending with
 * a fixed 1200m above sea level inflection point.
 */
class LocalIDWLapseAlgorithm : public InterpolationAlgorithm {
	public:
		LocalIDWLapseAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, IOManager& iom)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, iom) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
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
					const std::string& i_algo, IOManager& iom)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, iom), vecDataTA(), vecDataRH() {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		std::vector<double> vecDataTA, vecDataRH; ///<vectors of extracted TA and RH
};

/**
 * @class ILWRAlgorithm
 * @brief Incoming Long Wave Radiation interpolation algorithm.
 * Each ILWR is converted to an emissivity (using the local air temperature), interpolated using CST_LAPSE or IDW_LAPSE with
 * a fixed lapse rate and reconverted to ILWR.
 *
 * As a side effect, the user must have defined algorithms to be used for air temperature (since this is needed for
 * emissivity to ILWR conversion)
 */
class ILWRAlgorithm : public InterpolationAlgorithm {
	public:
		ILWRAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, IOManager& iom)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, iom), vecDataEA() {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		std::vector<double> vecDataEA; ///<vectors of extracted emissivities
};

/**
 * @class SimpleWindInterpolationAlgorithm
 * @brief Curvature/slope influenced  wind interpolation algorithm.
 * This is an implementation of the method described in (Liston & Elder, 2006): the wind speed and direction are
 * spatially interpolated using IDWLapseAlgorithm for the wind speed and using the user defined wind direction
 * interpolation algorithm. Then, the wind speed and direction fields are altered by wind weighting factors
 * and wind diverting factors (respectively) calculated from the local curvature and slope
 * (as taken from the DEM, see DEMObject).
 */
class SimpleWindInterpolationAlgorithm : public InterpolationAlgorithm {
	public:
		SimpleWindInterpolationAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, IOManager& iom)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, iom), vecDataVW(), vecDataDW() {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		std::vector<double> vecDataVW, vecDataDW; ///<vectors of extracted VW and DW
};

/**
 * @class USERInterpolation
 * @brief Reads user provided gridded data on the disk.
 * The grids are all in a directory that is given as the algorithm's argument. The files must be named
 * according to the following schema:
 * - {numeric date}_{capitalized meteo parameter}.asc, for example 200812011500_TA.asc
 * - Default_{capitalized meteo parameter}.asc for the grid to use when no measurements exist (which prevents
 * retrieving the date for the interpolation)
 * The meteo parameters can be found in \ref meteoparam "MeteoData". Example of use:
 * @code
 * TA::algorithms = USER
 * TA::user = ./meteo_grids
 * @endcode
 *
 */
class USERInterpolation : public InterpolationAlgorithm {
	public:
		USERInterpolation(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, IOManager& iom)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, iom), filename() {nrOfMeasurments=0;}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		std::string getGridFileName() const;
		std::string filename;
};

/**
 * @class SnowHNWInterpolation
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
 * HNW::algorithms = HNW_SNOW
 * HNW::hnw_snow = cst_lapse
 * HNW::cst_lapse = 0.0005 frac
 * @endcode
 *
 * @author Florian Kobierska, Jan Magnusson and Mathias Bavay
 */
class SnowHNWInterpolation : public InterpolationAlgorithm {
	public:
		SnowHNWInterpolation(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, IOManager& iom)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, iom) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
};

/**
 * @class OrdinaryKrigingAlgorithm
 * @brief Ordinary kriging.
 * This implements ordinary krigging (see https://secure.wikimedia.org/wikipedia/en/wiki/Kriging)
 * with user-selectable variogram model (see https://secure.wikimedia.org/wikipedia/en/wiki/Variogram).
 * The variogram and krigging
 * coefficients are re-computed fresh for each new grid (or time step).
 *
 * The variogram is currently computed with the current data (as 1/2*(X1-X2)^2), which makes it quite
 * uninteresting... The next improvement will consist in calculating the covariances (used to build the
 * variogram) from time series (thus reflecting the time-correlation between stations).
 *
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
					const std::string& i_algo, IOManager& iom)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, iom), variogram() {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	protected:
		void getDataForVariogram(std::vector<double> &distData, std::vector<double> &variData);
		bool computeVariogram();
		Fit1D variogram;
};


/**
 * @class LapseOrdinaryKrigingAlgorithm
 * @brief Ordinary kriging with detrending.
 * This is very similar to OrdinaryKrigingAlgorithm but performs detrending on the data.
 * @code
 * TA::algorithms = ODKRIG_LAPSE
 * TA::odkrig = SPHERICVARIO
 * @endcode
 *
 * @author Mathias Bavay
 */
class LapseOrdinaryKrigingAlgorithm : public OrdinaryKrigingAlgorithm {
	public:
		LapseOrdinaryKrigingAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, IOManager& iom)
			: OrdinaryKrigingAlgorithm(i_mi, i_vecArgs, i_algo, iom) {}
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
};

} //end namespace mio

#endif
