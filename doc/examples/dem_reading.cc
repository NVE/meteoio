#include "MeteoIO.h"

//This is a basic example of using as dem: the dem is read, the grid coordinates of a point given by its (lat,long) are retrieved
//and a sub-dem is extracted starting at these coordinates and extending dist_x and dist_y and writen out.
int main(void) {
	const double lat1=46.1592, lon1=8.12993;
	const double dist_x=700, dist_y=1200;
	DEMObject dem;
	IOHandler *io = NULL;
	unsigned int i,j;
	
	try {
		ConfigReader cfg("io.ini");
		io = new IOHandler(cfg);
	} catch (IOException& e){
		std::cout << "Problem with IOHandler creation, cause: " << e.what() << std::endl;
	}
	//reading dem
	io->readDEM(dem);
	//retrieving grid coordinates of a real world point
	dem.WGS84_to_grid(lat1, lon1, i,j);

	//computing grid distances from real world distances
	const unsigned int ncols = (unsigned int)ceil(dist_x/dem.cellsize);
	const unsigned int nrows = (unsigned int)ceil(dist_y/dem.cellsize);

	//extracting a sub-dem starting at the given coordinates and extending a given distance along x and along y
	DEMObject sub_dem(dem, i, j, ncols, nrows);
	//writing the sub-dem out
	io->write2DGrid(sub_dem,"sub_dem.dem");

	return 0;
}
