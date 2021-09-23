// SPDX-License-Identifier: LGPL-3.0-or-later
#include <iostream>
#include <stdexcept>
#include <vector>

#include <meteoio/MeteoIO.h>

const double grid_epsilon = 1e-3; //1e-4 is still too tight because of truncated results when writing data out for the ref.

int main(int argc, char** argv) {

	if (argc < 3)
		throw std::invalid_argument("Unexpected number of input arguments: A start and end date are required.");

	mio::Config cfg("./io.ini");
	mio::IOManager io(cfg);

	const double TZ = cfg.get("TIME_ZONE", "Input"); //get dates according to time zone
	mio::Date sdate, edate;
	mio::IOUtils::convertString(sdate, argv[1], TZ); //start date
	mio::IOUtils::convertString(edate, argv[2], TZ); //end date

	mio::DEMObject dem;
	io.readDEM(dem);

	//generate a couple of mockup grids to test grid resampling with:
	const bool gen_grids = false;
	if (gen_grids) {
		mio::Grid2DObject mock_grid(dem); //base grid is the provided DEM
		mio::Date dt(sdate);
		while (dt < edate) //grids for each day between cmdline args
		{
			int year, month, day, hour, min, sec;
			dt.getDate(year, month, day, hour, min, sec);
			mio::Date dt_hours(dt);
			static const size_t hours[3] = {0, 12, 18};
			for (int ii = 0; ii < 3; ++ii) { //a couple of grids per day
				dt_hours.setDate(year, month, day, hours[ii], 0, 0);
				mock_grid = ii; //set all coordinates to some dummy value
				io.write2DGrid(mock_grid, mio::MeteoGrids::TA, dt_hours);
			}
			dt = dt + 1.;
		}
		return 0;
	} //endif gen_grids


//	/* VSTATIONS */
//	std::vector< std::vector<mio::MeteoData> > mvec;
//	io.getMeteoData(sdate, edate, mvec); //get meteo data from a VStation according to INI
//	io.writeMeteoData(mvec);

	mio::Grid2DObject grid_ta;

//	/* RAW GRID READING  */
//	std::string info;
//	io.getMeteoData(sdate, dem, mio::MeteoData::TA, grid_ta, info); //assume 1st cmd line date exists as grid
//	io.write2DGrid(grid_ta, mio::MeteoGrids::TA, sdate);

	/* TEMPORAL GRID RESAMPLING */
	mio::Date inbetween(sdate);
	inbetween = inbetween + 1.7; //pick some date that is not there as raw grid
	std::cout << "Calling MeteoIO for date " << inbetween.toString(mio::Date::ISO) << std::endl;
	io.getMeteoData(inbetween, dem, mio::MeteoData::TA, grid_ta);
	io.write2DGrid(grid_ta, mio::MeteoGrids::TA, inbetween);

	return 0;
}
