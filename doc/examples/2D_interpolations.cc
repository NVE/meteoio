#include <iostream>
#include "MeteoIO.h"

//This is the a basic example of spatial interpolations. It does not check any exceptions, it only tries to be as c-like as possible
//provide date as ISO formatted, for example 2008-12-01T15:35:00 and 
//it will retrieve and interpolate the data for this date according to the io.ini configuration file
int main(int argc, char** argv) {
	(void)argc;

	Date_IO d1;
	std::vector<MeteoData> vecMeteo;
	std::vector<StationData> vecStation;

	//initializing the io handlers according to the config file
	ConfigReader cfg("io.ini");
	IOHandler *raw_io = NULL;
	BufferedIOHandler *io = NULL;
	raw_io = new IOHandler(cfg);
	io = new BufferedIOHandler(*raw_io, cfg);

	//reading the dem (necessary for several spatial interpolations algoritms)
	DEMObject dem;
	io->readDEM(dem);

	//reading the meteorological data for the requested time step
	IOUtils::convertString(d1,argv[1]);
	io->readMeteoData(d1, vecMeteo, vecStation);

	//performing spatial interpolations
	Meteo2DInterpolator mi(cfg, dem, vecMeteo, vecStation);
	Grid2DObject param;
	mi.interpolate(MeteoData::RH, param);

	//writing out the spatial interpolation
	io->write2DGrid(param,"rh.asc");

	delete io;
	delete raw_io;
	return 0;
}
