#ifndef GRIDLINEARRESAMPLING_H
#define GRIDLINEARRESAMPLING_H

#include <meteoio/gridResampling/GridResamplingAlgorithms.h>

namespace mio {

class GridLinearResampling : public GridResamplingAlgorithm {
	public:
		GridLinearResampling(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector< std::pair<std::string, std::string> >& vecArgs);

		void resample(const Date& date, const std::map<Date, Grid2DObject>& all_grids, Grid2DObject& resampled_grid);
		std::string toString() const;
};

} //end namespace mio

#endif
