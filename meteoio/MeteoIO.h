/***********************************************************************************/
/*  Copyright 2009-2011 WSL Institute for Snow and Avalanche Research    SLF-DAVOS */
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

#ifdef _MSC_VER
//VC++ complains that it can not generate an assignment operator
//for some classes (those having CONST members)
	#pragma warning (disable:4512)
#endif

//list in alphabetical order
//find meteoio -name "*.h" | sort
#include <meteoio/Array1D.h>
#include <meteoio/Array2D.h>
#include <meteoio/Array3D.h>
#include <meteoio/Array4D.h>
#include <meteoio/BufferedIOHandler.h>
#include <meteoio/Config.h>
#include <meteoio/Coords.h>
#include <meteoio/DataGenerator.h>
#include <meteoio/Date.h>
#include <meteoio/DEMObject.h>
#include <meteoio/exports.h>
#include <meteoio/GeneratorAlgorithms.h>
#include <meteoio/Graphics.h>
#include <meteoio/Grid2DObject.h>
#include <meteoio/Grid3DObject.h>
#include <meteoio/InterpolationAlgorithms.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/IOHandler.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOManager.h>
#include <meteoio/IOPlugin.h>
#include <meteoio/IOUtils.h>
//#include <meteoio/MainPage.h> //only for doxygen
#include <meteoio/MathOptim.h>
#include <meteoio/Matrix.h>
//#include <meteoio/MessageBoxX11.h>
#include <meteoio/Meteo1DInterpolator.h>
#include <meteoio/Meteo2DInterpolator.h>
#include <meteoio/MeteoData.h>

#include <meteoio/meteofilters/FilterBlock.h>
//skip all the filters' implementations header files
#include <meteoio/meteofilters/ProcessingBlock.h>
#include <meteoio/meteofilters/ProcessingStack.h>
#include <meteoio/meteofilters/WindowedFilter.h>

//#include <meteoio/MeteoIO.h>
#include <meteoio/meteolaws/Atmosphere.h>
#include <meteoio/meteolaws/Meteoconst.h>
#include <meteoio/meteolaws/Sun.h>
#include <meteoio/meteolaws/Suntrajectory.h>

#include <meteoio/MeteoProcessor.h>
//#include <meteoio/meteostats/libfit1DCore.h>
#include <meteoio/meteostats/libfit1D.h>
#include <meteoio/meteostats/libinterpol1D.h>
#include <meteoio/meteostats/libinterpol2D.h>

//skip all plugins' implementations header files
#include <meteoio/plugins/libsmet.h>

#include <meteoio/ResamplingAlgorithms.h>
#include <meteoio/ResamplingAlgorithms2D.h>
#include <meteoio/StationData.h>
#include <meteoio/Timer.h>

#endif
