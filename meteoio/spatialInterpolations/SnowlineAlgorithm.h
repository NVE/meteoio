/***********************************************************************************/
/*  Copyright 2020 ALPsolut.eu                                                     */
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
#ifndef SNOWLINE_ALGORITHM_H
#define SNOWLINE_ALGORITHM_H

#include <string>

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

class SnowlineAlgorithm : public InterpolationAlgorithm {
	public:
		SnowlineAlgorithm(const std::vector< std::pair<std::string, std::string> >& vecArgs,
		    const std::string& i_algo, const std::string& i_param, TimeSeriesManager& i_tsm,
		    GridsManager& i_gdm, Meteo2DInterpolator& i_mi);
		virtual double getQualityRating(const Date& i_date);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);

	private:
		typedef enum METHOD //assimilation method
		{
			BANDS, //!< elevation bands multiplied by from 0 (lowewst band) to 1 (highest band)
			CUTOFF //!< set all below snowline to zero
		} assimilation_method;

		void baseInterpol(const DEMObject& dem, Grid2DObject& grid);
		void assimilateCutoff(const DEMObject& dem, Grid2DObject& grid);
		void assimilateBands(const DEMObject& dem, Grid2DObject& grid);
		double readSnowlineFile();
		void getSnowline();
		void msg(const std::string& message);

		GridsManager& gdm;
		Meteo2DInterpolator& mi;
		std::string base_alg;
		double snowline;
		assimilation_method assim_method;
		std::string snowline_file;
		std::string where;
		double cutoff_val;
		double band_height;
		unsigned int band_no;
		bool quiet; //suppress warnings
};

} //end namespace mio

#endif //SNOWLINE_ALGORITHM_H
