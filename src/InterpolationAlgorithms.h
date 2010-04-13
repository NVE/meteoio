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

#include "MeteoIO.h"
#include "libinterpol1D.h"

#include <vector>
#include <string>
#include <set>

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
 * @section practical Practical use
 * Practically, the user
 * has to specify in his configuration file (typically io.ini), for each parameter to be interpolated, which
 * spatial interpolations algorithms should be considered. This is provided as a space separated list of keywords
 * (one per interpolation algorithm). Please notice that some algorithms may require extra arguments.
 * Then, each algorithm will be evaluated (through the use of its rating method) and receive a grade (that might
 * depend on the number of available data, the quality of the data, etc). The algorithm that receives the higher
 * score within the user list, will be used for interpolating the selected variable.
 *
 * @section keywords Available algorithms
 * The keywords defining the algorithms are the following:
 * - STD_PRESS: standard atmospheric pressure as a function of the elevation of each cell (see StandardPressureAlgorithm)
 * - CST: constant value in each cell (see ConstAlgorithm)
 * - CST_LAPSE: constant value reprojected to the elevation of the cell (see ConstLapseRateAlgorithm)
 * - IDW: Inverse Distance Weighting averaging (see IDWAlgorithm)
 * - IDW_LAPSE: Inverse Distance Weighting averaging with reprojection to the elevation of the cell (see IDWLapseAlgorithm)
 * - RH: the dew point temperatures are interpolated using IDW_LAPSE, then reconverted locally to relative humidity (see RHAlgorithm)
 * - WIND_CURV: the wind field (VW and DW) is interpolated using IDW_LAPSE and then altered depending on the local curvature and slope (taken from the DEM, see SimpleWindInterpolationAlgorithm)
 *
 * @section example Example of configuration file
 * Here is an example of the interpolation section of an configuration file (io.ini):
 * @code
 * [Interpolations2D]
 * TA::algorithms = IDW_LAPSE CST_LAPSE
 * TA::cst_lapse = -0.008
 * 
 * RH::algorithms = RH IDW_LAPSE CST_LAPSE CST
 * 
 * HNW::algorithms = IDW_LAPSE CST_LAPSE CST
 * 
 * VW::algorithms = IDW_LAPSE CST_LAPSE
 * 
 * P::algorithms = STD_PRESS
 * @endcode
 *
 * @section lapse Lapse rates
 * Several algorithms use elevation trends, currently modelled as a linear relation. The slope of this linear relation can
 * sometimes be provided by the end user (through his io.ini configuration file), otherwise it is computed from the data.
 * In order to bring slightly more robustness, if the correlation between the input data and the computed linear regression
 * is not good enought (below 0.7, as defined in Interpol2D::LinRegression), the same regression will get re-calculated
 * with one point less (cycling throught all the points). The best result (ie: highest correlation coefficient) will be
 * kept. If the final correlation coefficient is less than 0.7, a warning is displayed.
 *
 * @section dev_use Developer usage
 * From the developer's point of view, all that has to be done is instantiate a Meteo2DInterpolator object and call its 
 * Meteo2DInterpolator::interpolate method.
 * @code
 * 	std::vector<MeteoData> vecMeteo;
 * 	std::vector<StationData> vecStation;
 * 	ConfigReader cfg("io.ini");
 * 	DEMObject dem;
 * 
 * 	[...]
 * 
 * 	//performing spatial interpolations
 * 	Meteo2DInterpolator mi(cfg, dem, vecMeteo, vecStation);
 * 	Grid2DObject param;
 * 	mi.interpolate(MeteoData::RH, param);
 * @endcode
 *
 * @section biblio Bibliography
 * The interpolation algorithms have been inspired by the following papers:
 * - "A Meteorological Distribution System for High-Resolution Terrestrial Modeling (MicroMet)", Liston and Elder, Journal of Hydrometeorology 7 (2006), 217-234.
 * - "Simulating wind ﬁelds and snow redistribution using terrain-based parameters to model snow accumulation and melt over a semi-arid mountain catchment", Adam Winstral and Danny Marks, Hydrol. Process. 16 (2002), 3585– 3603. DOI: 10.1002/hyp.1238 [NOT YET IMPLEMENTED]
 * 
 * @author Mathias Bavay
 * @date   2010-04-12
 */

/**
 * @class InterpolationAlgorithm
 * @brief A class to perform 2D spatial interpolations. For more, see \ref interpol2d
 * @author Thomas Egger
 * @date   2010-04-01
*/
class InterpolationAlgorithm {

	public:
		InterpolationAlgorithm(const Meteo2DInterpolator& _mi, 
		                       const DEMObject& _dem,
		                       const std::vector<MeteoData>& _vecMeteo,
		                       const std::vector<StationData>& _vecStation,
		                       const std::vector<std::string>& _vecArgs,
		                       const std::string _algo)
			: mi(_mi), dem(_dem), vecMeteo(_vecMeteo), vecStation(_vecStation), vecArgs(_vecArgs), algo(_algo) 
		{
			if (vecMeteo.size() != vecStation.size())
				throw InvalidArgumentException("The two data and metadata vectors don't match in size!", AT);
		}
		virtual ~InterpolationAlgorithm() {}
		virtual double getQualityRating(const MeteoData::Parameters& param) = 0;
		virtual void calculate(const MeteoData::Parameters& param, Grid2DObject& grid) = 0;
 	protected:
		const Meteo2DInterpolator& mi;
		const DEMObject& dem;
		const std::vector<MeteoData>& vecMeteo;
		const std::vector<StationData>& vecStation;
		const std::vector<std::string>& vecArgs;
		const std::string algo;

		unsigned int getData(const MeteoData::Parameters& param, std::vector<double>& vecData) const;
		unsigned int getData(const MeteoData::Parameters& param, 
		                     std::vector<double>& vecData, std::vector<StationData>& vecMeta) const;
		unsigned int getStationAltitudes(const std::vector<StationData>& vecMeta, std::vector<double>& vecData) const;
		void printInfo(const MeteoData::Parameters& param, const unsigned int& stations_used) const;
};

class AlgorithmFactory {
	public:
		static InterpolationAlgorithm* getAlgorithm(const std::string& _algoname, 
                                                            const Meteo2DInterpolator& _mi,
		                                            const DEMObject& _dem,
		                                            const std::vector<MeteoData>& _vecMeteo,
		                                            const std::vector<StationData>& _vecStation,
		                                            const std::vector<std::string>& _vecArgs);

		static std::set<std::string> setAlgorithms; ///<all algorithms that are configured
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static setAlgorithms
};

/**
 * @class ConstAlgorithm
 * @brief Constant filling interpolation algorithm. 
 * Fill the grid with the average of the inputs for this parameter.
 */
class ConstAlgorithm : public InterpolationAlgorithm {
	public:
		ConstAlgorithm(const Meteo2DInterpolator& _mi, 
		               const DEMObject& _dem,
		               const std::vector<MeteoData>& _vecMeteo,
		               const std::vector<StationData>& _vecStation,
		               const std::vector<std::string>& _vecArgs,
		               const std::string _algo)
			: InterpolationAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs, _algo) {}
		virtual double getQualityRating(const MeteoData::Parameters& param);
		virtual void calculate(const MeteoData::Parameters& param, Grid2DObject& grid);
};

/**
 * @class StandardPressureAlgorithm
 * @brief Standard atmospheric pressure interpolation algorithm.
 * Fill the grid with the standard atmosphere's pressure, depending on the local elevation.
 */
class StandardPressureAlgorithm : public InterpolationAlgorithm {
	public:
		StandardPressureAlgorithm(const Meteo2DInterpolator& _mi, 
		                          const DEMObject& _dem,
		                          const std::vector<MeteoData>& _vecMeteo,
		                          const std::vector<StationData>& _vecStation,
		                          const std::vector<std::string>& _vecArgs,
		                          const std::string _algo)
			: InterpolationAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs, _algo) {}
		virtual double getQualityRating(const MeteoData::Parameters& param);
		virtual void calculate(const MeteoData::Parameters& param, Grid2DObject& grid);
};

/**
 * @class ConstLapseRateAlgorithm
 * @brief Constant filling with elevation lapse rate interpolation algorithm.
 * Assuming that average values occured at the average of the elevations, the grid is filled with average values
 * reprojected to real grid elevation according to a user specified lapse rate. The lapse rate has
 * to be provided as an extra argument, otherwise the standard -0.0065 K/m is used (which only makes sense for temperatures)
 */
class ConstLapseRateAlgorithm : public InterpolationAlgorithm {
	public:
		ConstLapseRateAlgorithm(const Meteo2DInterpolator& _mi, 
		                        const DEMObject& _dem,
		                        const std::vector<MeteoData>& _vecMeteo,
		                        const std::vector<StationData>& _vecStation,
		                        const std::vector<std::string>& _vecArgs,
		                        const std::string _algo)
			: InterpolationAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs, _algo) {}
		virtual double getQualityRating(const MeteoData::Parameters& param);
		virtual void calculate(const MeteoData::Parameters& param, Grid2DObject& grid);
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
		IDWAlgorithm(const Meteo2DInterpolator& _mi, 
		             const DEMObject& _dem,
		             const std::vector<MeteoData>& _vecMeteo,
		             const std::vector<StationData>& _vecStation,
		             const std::vector<std::string>& _vecArgs,
		             const std::string _algo)
			: InterpolationAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs, _algo) {}
		virtual double getQualityRating(const MeteoData::Parameters& param);
		virtual void calculate(const MeteoData::Parameters& param, Grid2DObject& grid);
};

/**
 * @class IDWLapseAlgorithm
 * @brief Inverse Distance Weighting interpolation algorithm with elevation detrending/reprojection.
 * The input data is projected to a reference elevation and spatially interpolated using an Inverse Distance
 * Weighting interpolation algorithm (see IDWAlgorithm). Then, each value is reprojected to the real
 * elevation of the relative cell.
 */
class IDWLapseAlgorithm : public InterpolationAlgorithm {
	public:
		IDWLapseAlgorithm(const Meteo2DInterpolator& _mi, 
		                  const DEMObject& _dem,
		                  const std::vector<MeteoData>& _vecMeteo,
		                  const std::vector<StationData>& _vecStation,
		                  const std::vector<std::string>& _vecArgs,
		                  const std::string _algo)
			: InterpolationAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs, _algo) {}
		virtual double getQualityRating(const MeteoData::Parameters& param);
		virtual void calculate(const MeteoData::Parameters& param, Grid2DObject& grid);
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
		RHAlgorithm(const Meteo2DInterpolator& _mi, 
		            const DEMObject& _dem,
		            const std::vector<MeteoData>& _vecMeteo,
		            const std::vector<StationData>& _vecStation,
		            const std::vector<std::string>& _vecArgs,
		            const std::string _algo)
			: InterpolationAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs, _algo) {}
		virtual double getQualityRating(const MeteoData::Parameters& param);
		virtual void calculate(const MeteoData::Parameters& param, Grid2DObject& grid);
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
		SimpleWindInterpolationAlgorithm(const Meteo2DInterpolator& _mi,
		                                 const DEMObject& _dem,
		                                 const std::vector<MeteoData>& _vecMeteo,
		                                 const std::vector<StationData>& _vecStation,
		                                 const std::vector<std::string>& _vecArgs,
		                                 const std::string _algo)
			: InterpolationAlgorithm(_mi, _dem, _vecMeteo, _vecStation, _vecArgs, _algo) {}
		virtual double getQualityRating(const MeteoData::Parameters& param);
		virtual void calculate(const MeteoData::Parameters& param, Grid2DObject& grid);
};

#endif
