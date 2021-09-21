#ifndef GRIDTIMESERIESRESAMPLING_H
#define GRIDTIMESERIESRESAMPLING_H

#include <string>
#include <utility>
#include <vector>

#include <meteoio/gridResampling/GridResamplingAlgorithms.h>

namespace mio {

class GridTimeseriesResampling : public GridResamplingAlgorithm {
	public:
		GridTimeseriesResampling(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs);

		void resample(const Date& date, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid);
		std::string toString() const;

	private:
		std::vector< std::pair<std::string, std::string> > vecArgs_;
		std::string base_algorithm_;
};

} //end namespace mio

#endif
