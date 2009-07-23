#ifndef __METEO2DINTERPOLATOR_H__
#define __METEO2DINTERPOLATOR_H__

#include "Array2D.h"
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
		Meteo2DInterpolator(const DEMObject& dem, const vector<MeteoData>& vecData, const vector<StationData>& vecMeta);

		/**
		* @brief This function calls the interpolation class for each individual meteo parameter. 
		*        It also builds a list of valid data sources for the given parameter.
		*/
		void interpolate(Grid2DObject& nswc, Grid2DObject& rh, Grid2DObject& ta, Grid2DObject& vw, Grid2DObject& p);
		void interpolate(Grid2DObject& nswc, Grid2DObject& rh, Grid2DObject& ta,
				 Grid2DObject& vw, Grid2DObject& p, Grid2DObject& iswr/*, Grid2DObject& ea*/);
		void interpolate(Array2D<double>& nswc, Array2D<double>& rh, Array2D<double>& ta, 
				      Array2D<double>& vw, Array2D<double>& p);
		void interpolate(Array2D<double>& nswc, Array2D<double>& rh, Array2D<double>& ta, 
				      Array2D<double>& vw, Array2D<double>& p, Array2D<double>& iswr/*, Array2D<double>& ea*/);

	private:
		void interpolateP(Array2D<double>& p);
		void interpolateNSWC(Array2D<double>& nswc);
		void interpolateTA(Array2D<double>& ta);
		void interpolateRH(Array2D<double>& rh, Array2D<double>& ta);
		void interpolateVW(Array2D<double>& vw);
		void interpolateDW(Array2D<double>& dw);
		void interpolateISWR(Array2D<double>& iswr);
		//void interpolateEA(Array2D<double>& ea);

		const DEMObject& dem;

		vector<MeteoData> SourcesData;
		vector<StationData> SourcesMeta;
};

#endif
