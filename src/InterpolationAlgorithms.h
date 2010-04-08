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
