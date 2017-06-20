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
#include <meteoio/TimeSeriesManager.h>
#include <meteoio/GridsManager.h>
#include <meteoio/Meteo2DInterpolator.h>
#include <meteoio/meteoStats/libfit1D.h>

#include <vector>
#include <string>

namespace mio {

class Meteo2DInterpolator; // forward declaration, cyclic header include

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
		                       const std::vector< std::pair<std::string, std::string> >& /*vecArgs*/,
		                       const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager, const std::string& i_param) :
		                      algo(i_algo), mi(i_mi), tsmanager(i_tsmanager), gridsmanager(i_gridsmanager), date(0., 0), vecMeteo(), vecData(),
		                      vecMeta(), info(), param(i_param), nrOfMeasurments(0), user_lapse(IOUtils::nodata), trend_min_alt(-1e4), trend_max_alt(1e4), is_frac(false), is_soft(false) {}
		virtual ~InterpolationAlgorithm() {}
		
		//if anything is not ok (wrong parameter for this algo, insufficient data, etc) -> return zero
		virtual double getQualityRating(const Date& i_date) = 0;
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid) = 0;
		std::string getInfo() const;
		const std::string algo;

 	protected:
		std::vector<double> getData(const Date& i_date, const std::string& i_param);
		size_t getData(const Date& i_date, const std::string& i_param,
		               std::vector<double>& o_vecData, std::vector<StationData>& o_vecMeta);
		void simpleWindInterpolate(const DEMObject& dem, const std::vector<double>& vecDataVW, const std::vector<double>& vecDataDW, Grid2DObject &VW, Grid2DObject &DW, const double& scale, const double& alpha);
		Fit1D getTrend(const std::vector<double>& vecAltitudes, const std::vector<double>& vecDat) const;
		void setTrendParams(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		static std::vector<double> getStationAltitudes(const std::vector<StationData>& i_vecMeta);
		void detrend(const Fit1D& trend, const std::vector<double>& vecAltitudes, std::vector<double> &vecDat) const;
		void retrend(const DEMObject& dem, const Fit1D& trend, Grid2DObject &grid) const;

		template <class T> void parseArg(const std::pair< std::string, std::string>& arg, T& val) const {
			if (!IOUtils::convertString(val, arg.second))
				throw InvalidArgumentException("Can not parse argument "+arg.first+"::"+arg.second+"' for algorithm " + algo, AT);
		}

		Meteo2DInterpolator& mi;
		TimeSeriesManager& tsmanager;
		GridsManager& gridsmanager;
		Date date;

		std::vector<MeteoData> vecMeteo;
		std::vector<double> vecData; ///<store the measurement for the given parameter
		std::vector<StationData> vecMeta; ///<store the station data for the given parameter
		std::ostringstream info; ///<to store some extra information about the interplation process
		const std::string param; ///<the parameter that we will interpolate
		size_t nrOfMeasurments; ///<Number of stations that have been used, so this can be reported to the user
		double user_lapse; ///<when detrending the data, this is a user provided lapse_rate
		double trend_min_alt, trend_max_alt; ///<if these are provided, the detrending/retrending will be bound by a minimum and/or maximum altitude
		bool is_frac, is_soft; ///<is the lapse rate to be interpreted as fractional? Should it be used as fallback (is_soft) or it is mandatory?
};

class AlgorithmFactory {
	public:
		static InterpolationAlgorithm* getAlgorithm(const std::string& i_algoname,
		                                            Meteo2DInterpolator& i_mi,
		                                            const std::vector< std::pair<std::string, std::string> >& vecArgs, TimeSeriesManager& tsm, GridsManager& gdm, const std::string& i_param);
};

} //end namespace mio

#endif
