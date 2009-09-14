#include <iostream>
#include "MeteoIO.h"

using namespace std;

int main(int argc, char** argv) {
	(void)argc;
	Date_IO d1;
	
	if(argc<2){
		printf("Error: not enough arguments !\n");
		exit(-1);
	}
	convertString(d1,argv[1]);
	
	vector<MeteoData> vecMeteo;
	vector<StationData> vecStation;
	
	IOHandler *raw_io = NULL;
	BufferedIOHandler *io = NULL;

	try {
		ConfigReader cfg("io.ini");
		raw_io = new IOHandler(cfg);
		io = new BufferedIOHandler(*raw_io, cfg);
	} catch (IOException& e){
		cout << "Problem with IOHandler creation, cause: " << e.what() << endl;
	}
	
	try {
		io->readMeteoData(d1, vecMeteo, vecStation);
	} catch (IOException& e){
		cout << "Problem when reading data, cause: " << e.what() << endl;
	}
	
	//writing some data out in order to prove that it really worked!
	for (unsigned int ii=0; ii<vecMeteo.size(); ii++) {
		cout << "---------- Station: " << (ii+1) << " / " << vecStation.size() << endl;
		cout << vecStation[ii].toString() << endl;
		cout << vecMeteo[ii].toString() << endl;
		//cout << "  Name: " << vecStation[ii].getStationName() << endl;
		//cout << "  RH: " << vecMeteo[ii].rh << endl;
	}

	//And now, doing spatial interpolations
	DEMObject dem;
	io->readDEM(dem);
	dem.update(DEMObject::CORR);
	
	//convert to local grid coordinates, an elegant way of dealing with multiple coordinates systems inputs
	dem.xllcorner = 0.;
	dem.yllcorner = 0.;
	for (unsigned i = 0; i < vecMeteo.size() ; ++i) {
		//setup reference station, convert coordinates to local grid
		IOUtils::WGS84_to_local(dem.latitude, dem.longitude, vecStation[i].latitude, vecStation[i].longitude, vecStation[i].eastCoordinate, vecStation[i].northCoordinate);
	}

	Grid2DObject    p(dem.ncols, dem.nrows, dem.xllcorner, dem.yllcorner, dem.latitude, dem.longitude, dem.cellsize);
	Grid2DObject nswc(dem.ncols, dem.nrows, dem.xllcorner, dem.yllcorner, dem.latitude, dem.longitude, dem.cellsize);
	Grid2DObject   vw(dem.ncols, dem.nrows, dem.xllcorner, dem.yllcorner, dem.latitude, dem.longitude, dem.cellsize);
	Grid2DObject   rh(dem.ncols, dem.nrows, dem.xllcorner, dem.yllcorner, dem.latitude, dem.longitude, dem.cellsize);
	Grid2DObject   ta(dem.ncols, dem.nrows, dem.xllcorner, dem.yllcorner, dem.latitude, dem.longitude, dem.cellsize);
	Meteo2DInterpolator mi(dem, vecMeteo, vecStation);

	mi.interpolate(nswc, rh, ta, vw, p);
	
	cout << "Writing the Grids to *.2d files" << endl;
	io->write2DGrid(ta, "output/ta.2d");
	io->write2DGrid(p, "output/p.2d");
	io->write2DGrid(vw, "output/vw.2d");
	io->write2DGrid(nswc, "output/nswc.2d");
	io->write2DGrid(rh, "output/rh.2d");
	
	cout << "Writing the Grids was successful" << endl;

	delete io;
	delete raw_io;

	return 0;
}
