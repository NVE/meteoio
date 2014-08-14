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
//find meteoio -name "*.h" | sort | xargs -i echo "#include <{}>"
#include <meteoio/BufferedIOHandler.h>
#include <meteoio/Config.h>

#include <meteoio/dataClasses/Array1D.h>
#include <meteoio/dataClasses/Array2D.h>
#include <meteoio/dataClasses/Array3D.h>
#include <meteoio/dataClasses/Array4D.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/dataClasses/Date.h>
#include <meteoio/dataClasses/DEMObject.h>
#include <meteoio/dataClasses/Grid2DObject.h>
#include <meteoio/dataClasses/Grid3DObject.h>
#include <meteoio/dataClasses/Matrix.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/dataClasses/StationData.h>

#include <meteoio/DataGenerator.h>
#include <meteoio/exports.h>
#include <meteoio/FileUtils.h>
#include <meteoio/GeneratorAlgorithms.h>
#include <meteoio/Graphics.h>
#include <meteoio/InterpolationAlgorithms.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/IOHandler.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOManager.h>
#include <meteoio/IOUtils.h>
//#include <meteoio/MainPage.h> //only for doxygen
#include <meteoio/MathOptim.h>
//#include <meteoio/MessageBoxX11.h>
#include <meteoio/Meteo1DInterpolator.h>
#include <meteoio/Meteo2DInterpolator.h>

#include <meteoio/meteoFilters/FilterBlock.h>
//skip all the filters' implementations header files
#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <meteoio/meteoFilters/ProcessingStack.h>
//#include <meteoio/meteoFilters/template.h>
#include <meteoio/meteoFilters/WindowedFilter.h>

//#include <meteoio/MeteoIO.h>

#include <meteoio/meteoLaws/Atmosphere.h>
#include <meteoio/meteoLaws/Meteoconst.h>
#include <meteoio/meteoLaws/Sun.h>
#include <meteoio/meteoLaws/Suntrajectory.h>

#include <meteoio/MeteoProcessor.h>
//#include <meteoio/meteoStats/libfit1DCore.h>
#include <meteoio/meteoStats/libfit1D.h>
#include <meteoio/meteoStats/libinterpol1D.h>
#include <meteoio/meteoStats/libinterpol2D.h>

//skip all plugins' implementations header files
#include <meteoio/plugins/libncpp.h>
#include <meteoio/plugins/libsmet.h>

#include <meteoio/ResamplingAlgorithms2D.h>
#include <meteoio/ResamplingAlgorithms.h>
#include <meteoio/Timer.h>

#endif
