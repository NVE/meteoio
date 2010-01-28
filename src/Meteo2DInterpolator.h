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

#ifndef __METEO2DINTERPOLATOR_H__
#define __METEO2DINTERPOLATOR_H__

#include "Grid2DObject.h"
#include "libinterpol2D.h"
#include "MeteoData.h"
#include "StationData.h"
#include "IOUtils.h"
#include "IOExceptions.h"
#include "DEMObject.h"

#include <vector>
#include <map>

/**
 * @class Meteo2DInterpolator
 * @brief A class to spatially interpolate meteo parameters
 *
 * @author Mathias Bavay and Thomas Egger
 * @date   2010-01-14
 */
#ifdef _POPC_
class Meteo2DInterpolator : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class Meteo2DInterpolator {
#endif
 	public:
		/**
		* @brief Constructor.
		*/
		Meteo2DInterpolator(const ConfigReader& _cfg, const DEMObject& _dem, 
						const std::vector<MeteoData>& _vecData, const std::vector<StationData>& _vecMeta);


		/**
		 * @brief A generic function that can interpolate for any given MeteoData member variable
		 * 
		 * @param meteoparam Any MeteoData member variable as specified in the 
		 * 				 enum MeteoData::Parameters (e.g. MeteoData::TA)
		 * @param result A Grid2DObject that will be filled with the interpolated data
		 */
		void interpolate(const MeteoData::Parameters& meteoparam, Grid2DObject& result);

	private:
		const ConfigReader& cfg; ///< Reference to ConfigReader object, initialized during construction
		const DEMObject& dem;    ///< Reference to DEMObject object, initialized during construction
		const std::vector<MeteoData>& vecData;  ///< Reference to a vec of MeteoData, initialized during construction
		const std::vector<StationData>& vecMeta;///< Reference to a vec of StationData, initialized during construction
		
		std::map< std::string, std::vector<std::string> > mapAlgorithms; //per parameter interpolation algorithms

		unsigned int getAlgorithmsForParameter(const std::string& parname, std::vector<std::string>& vecAlgorithms);
		unsigned int getArgumentsForAlgorithm(const std::string& keyname, std::vector<std::string>& vecArguments);

		/*LEGACY*/
		//These methods are exposed until we offer a better API for requesting spatial interpolations
		void interpolateP(Grid2DObject& p);
		void interpolateHNW(Grid2DObject& hnw);
		void interpolateTA(Grid2DObject& ta);
		void interpolateRH(Grid2DObject& rh, Grid2DObject& ta);
		void interpolateVW(Grid2DObject& vw);
		void interpolateDW(Grid2DObject& dw);
		void interpolateISWR(Grid2DObject& iswr);
		void interpolateLWR(Grid2DObject& lwr);
};

#endif
