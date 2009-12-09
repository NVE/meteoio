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

/**
 * @class Meteo2DInterpolator
 * @brief A class to spatially interpolate Alpine3D's meteo parameters
 * @author Mathias Bavay
 * @date   2009-01-20
 */
class Meteo2DInterpolator {
	public:
		/**
		* @brief Constructor. It builds a vector of input data and metadata, merging meteo1D and meteo2D input files
		*/
		Meteo2DInterpolator(const DEMObject& _dem, const std::vector<MeteoData>& vecData, const std::vector<StationData>& vecMeta);

		/**
		* @brief This function calls the interpolation class for each individual meteo parameter. 
		*        It also builds a list of valid data sources for the given parameter.
		*/
		void interpolate(Grid2DObject& hnw, Grid2DObject& rh, Grid2DObject& ta, Grid2DObject& vw, Grid2DObject& p);
		void interpolate(Grid2DObject& hnw, Grid2DObject& rh, Grid2DObject& ta,
				 Grid2DObject& vw, Grid2DObject& p, Grid2DObject& iswr/*, Grid2DObject& ea*/);

	private:
		void interpolateP(Grid2DObject& p);
		void interpolateHNW(Grid2DObject& hnw);
		void interpolateTA(Grid2DObject& ta);
		void interpolateRH(Grid2DObject& rh, Grid2DObject& ta);
		void interpolateVW(Grid2DObject& vw);
		void interpolateDW(Grid2DObject& dw);
		void interpolateISWR(Grid2DObject& iswr);
		void interpolateLWR(Grid2DObject& lwr);

		const DEMObject& dem;

		std::vector<MeteoData> SourcesData;
		std::vector<StationData> SourcesMeta;
};

#endif
