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

#ifndef __METEOIO_H__
#define __METEOIO_H__

//list in alphabetical order
#include <meteoio/Array2D.h>
#include <meteoio/Array3D.h>
#include <meteoio/Array.h>
#include <meteoio/BufferedIOHandler.h>
#include <meteoio/Config.h>
#include <meteoio/Coords.h>
#include <meteoio/Date.h>
#include <meteoio/Timer.h>
#include <meteoio/DEMObject.h>
#include <meteoio/DynamicLibrary.h>
#include <meteoio/FilterAlgorithms.h>
#include <meteoio/FilterProperties.h>
#include <meteoio/Grid2DObject.h>
#include <meteoio/Grid3DObject.h>
#include <meteoio/InterpolationAlgorithms.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOPlugin.h>
#include <meteoio/IOUtils.h>
#include <meteoio/libfit1D.h>
#include <meteoio/libinterpol1D.h>
#include <meteoio/libinterpol2D.h>
#include <meteoio/Matrix.h>
#include <meteoio/Meteo1DInterpolator.h>
#include <meteoio/Meteo2DInterpolator.h>
#include <meteoio/MeteoData.h>
//#include <meteoio/MeteoFilter.h> //HACK: is it obsolete?
#include <meteoio/MeteoProcessor.h> //HACK: is it obsolete?
#include <meteoio/meteofilters/FilterBlock.h>
#include <meteoio/meteofilters/ProcessingBlock.h>
#include <meteoio/meteofilters/ProcessingStack.h>
#include <meteoio/ResamplingAlgorithms.h> //HACK: is it obsolete?
#include <meteoio/StationData.h>
#include <meteoio/IOManager.h>

#include <meteoio/meteolaws/Atmosphere.h>
#include <meteoio/meteolaws/Suntrajectory.h>
#include <meteoio/meteolaws/Sun.h>
#include <meteoio/meteolaws/Meteoconst.h>

#ifdef _POPC_
#include <meteoio/IOHandler.ph>
#include <meteoio/marshal_meteoio.h>
#else
#include <meteoio/IOHandler.h>
#endif

#endif
