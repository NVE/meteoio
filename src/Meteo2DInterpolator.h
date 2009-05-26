#ifndef __METEO2DINTERPOLATOR_H__
#define __METEO2DINTERPOLATOR_H__

#include "Array2D.h"
#include "libinterpol2D.h"
#include "MeteoData.h"
#include "StationData.h"
#include "slfexceptions.h"

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
  Meteo2DInterpolator(const Grid2DObject& dem, 
		      const vector<MeteoData>& vecData, 
		      const vector<StationData>& vecMeta);
  
  /**
   * @brief This function calls the interpolation class for each individual meteo parameter. 
   *        It also builds a list of valid data sources for the given parameter.
   */
  void interpolate(CArray2D<double>& nswc, CArray2D<double>& rh, CArray2D<double>& ta, CArray2D<double>& vw, CArray2D<double>& p);
  void interpolate(CArray2D<double>& nswc, CArray2D<double>& rh, CArray2D<double>& ta, CArray2D<double>& vw, CArray2D<double>& p, CArray2D<double>& iswr/*, CArray2D<double>& ea*/);

 private:
  void interpolateP(CArray2D<double>& p);
  void interpolateNSWC(CArray2D<double>& nswc);
  void interpolateTA(CArray2D<double>& ta);
  void interpolateRH(CArray2D<double>& rh, CArray2D<double>& ta);
  void interpolateVW(CArray2D<double>& vw);
  void interpolateISWR(CArray2D<double>& iswr);
  //void interpolateEA(CArray2D<double>& ea);

  const Grid2DObject& dem; //HACK: This doesn't hold for PAROC

  vector<MeteoData> SourcesData;
  vector<StationData> SourcesMeta;
};

#endif
