#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio

//This is a basic example of using as dem: the dem is read, the grid coordinates of a point given by its (lat,long) are retrieved
//and a sub-dem is extracted starting at these coordinates and extending dist_x and dist_y and written out.
int main(void) {
	DEMObject dem;
	Config cfg("io.ini");
	IOManager io(cfg);

	//reading dem
	dem.setUpdatePpt(DEMObject::SLOPE);
	io.readDEM(dem);

	//writing some statistics about this dem
	//dem.grid2D.getMin() scans the DEM grid to get the min, while dem.min_altitude is cached and therefore very cheap
	//The raw content of the 2D grids can also be accessed, for example dem.grid2D.getMin(IOUtils::RAW_NODATA). In this case, there would be no interpretation of some values as nodata.
	std::cout << "DEM information: \n";
	std::cout << "\tmin=" << dem.grid2D.getMin() << " max=" << dem.grid2D.getMax() << " mean=" << dem.grid2D.getMean() << "\n";
	std::cout << "\tmin slope=" << dem.min_slope << " max slope=" << dem.max_slope << std::endl;


	Grid2DObject slope(dem.ncols, dem.nrows, dem.cellsize, dem.llcorner, dem.slope);
	io.write2DGrid(slope, MeteoGrids::SLOPE, Date(0.));
	Grid2DObject azi(dem.ncols, dem.nrows, dem.cellsize, dem.llcorner, dem.azi);
	io.write2DGrid(azi, MeteoGrids::AZI, Date(0.));

	//retrieving grid coordinates of a real world point
	Coords point;
	point.copyProj(dem.llcorner); //we use the same projection parameters as the DEM
	point.setLatLon(46.232103, 7.362185, IOUtils::nodata);
	dem.gridify(point);

	//computing grid distances from real world distances
	const double dist_x=70000., dist_y=120000.;
	const unsigned int ncols = (unsigned int)ceil(dist_x/dem.cellsize);
	const unsigned int nrows = (unsigned int)ceil(dist_y/dem.cellsize);

	//extracting a sub-dem starting at the given coordinates and extending a given distance along x and along y
	DEMObject sub_dem(dem, point.getGridI(), point.getGridJ(), ncols, nrows);

	//writing the sub-dem out
	io.write2DGrid(sub_dem, MeteoGrids::DEM, Date(0.));

	return 0;
}
