#include <iostream>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio
using namespace std;

const double epsilon = 1e-3; //1e-4 is still too tight because of truncated results when writing data out for the ref.

int main(int /*argc*/, char** argv) {
	Date d1;

	//initializing the io handlers according to the config file
	Config cfg("io.ini");
	IOManager io(cfg);

	//reading the dem (necessary for several spatial interpolations algoritms)
	DEMObject dem;
	io.readDEM(dem);

	//we assume that the time given on the command line is in TZ=+1
	IOUtils::convertString(d1,argv[1], 1.);

	std::string date_str = d1.toString(Date::ISO);
	std::replace( date_str.begin(), date_str.end(), ':', '.');

	//performing spatial interpolations
	Grid2DObject param, ref;
	io.interpolate(d1, dem, MeteoData::TA, param);
	io.read2DGrid(ref, date_str+"_TA_ref.asc");
	if(ref.grid2D.checkEpsilonEquality(param.grid2D, epsilon)==false) {
		cout << "TA grids don't match!\n"; return EXIT_FAILURE;
	}

	io.interpolate(d1, dem, MeteoData::HNW, param);
	io.read2DGrid(ref, date_str+"_HNW_ref.asc");
	if(ref.grid2D.checkEpsilonEquality(param.grid2D, epsilon)==false) {
		cout << "HNW grids don't match!\n"; return EXIT_FAILURE;
	}

	io.interpolate(d1, dem, MeteoData::RH, param);
	io.read2DGrid(ref, date_str+"_RH_ref.asc");
	if(ref.grid2D.checkEpsilonEquality(param.grid2D, epsilon)==false) {
		cout << "RH grids don't match!\n"; return EXIT_FAILURE;
	}
	//io.write2DGrid(param, MeteoGrids::RH, d1);
	io.interpolate(d1, dem, MeteoData::RSWR, param);
	io.read2DGrid(ref, date_str+"_RSWR_ref.asc");
	if(ref.grid2D.checkEpsilonEquality(param.grid2D, epsilon)==false) {
		cout << "RSWR grids don't match!\n"; return EXIT_FAILURE;
	}

	io.write2DGrid(param, MeteoGrids::RSWR, d1); //trying to write one grid out

	return 0;
}
