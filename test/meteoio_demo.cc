#include <iostream>
#include "MeteoIO.h"

using namespace std;

int main(int argc, char** argv) {
	(void)argc;
	Date_IO d1;
	convertString(d1,argv[1]);
	
	vector<MeteoData> vecMeteo;
	vector<StationData> vecStation;
	
	IOHandler *io=NULL; //Initialization vital!
	
	try {
		io = new IOHandler("io.ini");
	} catch (exception& e){
		cout << "Problem with IOHandler creation, cause: " << e.what() << endl;
	}
	
	try {
		io->readMeteoData(d1, vecMeteo, vecStation);
	} catch (exception& e){
		cout << "Problem when reading data, cause: " << e.what() << endl;
	}
	
	//writing some data out in order to prove that it really worked!
	for (unsigned int ii=0; ii<vecMeteo.size(); ii++) {
		cout << "---------- Station: " << (ii+1) << " / " << vecStation.size() << endl;
		cout << "  Name: " << vecStation[ii].getStationName() << endl;
		cout << "  RH: " << vecMeteo[ii].rh << endl;
	}
	
	//And now, doing spatial interpolations
	DEMObject dem;
	io->readDEM(dem);

	int nx=dem.ncols, ny=dem.nrows;
	dem.update();
	
	//convert to local grid coordinates
	/*dem.xllcorner = 0.;
	dem.yllcorner = 0.;
	for (unsigned i = 0; i < vecMeteo.size() ; ++i) {
		//setup reference station, convert coordinates to local grid
		IOUtils::WGS84_to_local(dem.latitude, dem.longitude, vecStation[i].latitude, vecStation[i].longitude, vecStation[i].eastCoordinate, vecStation[i].northCoordinate);
	}*/

	Array2D<double> p(nx,ny), nswc(nx,ny), vw(nx,ny), rh(nx,ny), ta(nx,ny);
	Meteo2DInterpolator mi(dem, vecMeteo, vecStation);

	mi.interpolate(nswc, rh, ta, vw, p);
	
	cout << "Writing the Grids to *.2d files" << endl;
	io->write2DGrid(ta, dem.xllcorner, dem.yllcorner, dem.cellsize, "output/ta.2d");
	io->write2DGrid(p, dem.xllcorner, dem.yllcorner, dem.cellsize, "output/p.2d");
	io->write2DGrid(vw, dem.xllcorner, dem.yllcorner, dem.cellsize, "output/vw.2d");
	io->write2DGrid(nswc, dem.xllcorner, dem.yllcorner, dem.cellsize, "output/nswc.2d");
	io->write2DGrid(rh, dem.xllcorner, dem.yllcorner, dem.cellsize, "output/rh.2d");
	
	cout << "Writing the Grids was successful" << endl;
	return 0;
}
