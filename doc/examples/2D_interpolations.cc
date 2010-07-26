#include <iostream>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio

//This is the a basic example of spatial interpolations. It does not check any exceptions, it only tries to be as c-like as possible
//provide date as ISO formatted, for example 2008-12-01T15:35:00 and 
//it will retrieve and interpolate the data for this date according to the io.ini configuration file
int main(int /*argc*/, char** argv) {
	Date d1;
	std::vector<MeteoData> vecMeteo;
	std::vector<StationData> vecStation;

	//initializing the io handlers according to the config file
	Config cfg("io.ini");
	IOHandler raw_io(cfg);
	BufferedIOHandler io(raw_io, cfg);

	//reading the dem (necessary for several spatial interpolations algoritms)
	DEMObject dem;
	io.readDEM(dem);

	//we assume that the time given on the command line is in TZ=+1
	d1.setTimeZone(1.);
	IOUtils::convertString(d1,argv[1]);
	io.readMeteoData(d1, vecMeteo, vecStation);

	//performing spatial interpolations
	Meteo2DInterpolator mi(cfg, dem, vecMeteo, vecStation);
	Grid2DObject param;
	mi.interpolate(MeteoData::RH, param);

	//writing out the spatial interpolation
	io.write2DGrid(param,"rh.asc");

	return 0;
}
