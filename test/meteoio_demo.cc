#include <iostream>
#include "MeteoIO.h"

int main() {
  Date_IO d1(2009,01,01,18,00);

  vector<MeteoData> vecMeteo;
  vector<StationData> vecStation;

  IOInterface *ioTest=NULL; //Initialization vital!

  try {
    ioTest = new IOHandler("io.ini");
  } catch (exception& e){
    cout << "Problem with IOHandler creation, cause: " << e.what() << endl;
  }

 try {
    ioTest->readMeteoData(d1, vecMeteo, vecStation);
 } catch (exception& e){
    cout << "Problem when reading data, cause: " << e.what() << endl;
  }

  //writing some data out in order to prove that it really worked!
  for (unsigned int ii=0; ii<vecMeteo.size(); ii++) {
      cout << "---------- Station: " << (ii+1) << " / " << vecStation.size() << endl;
      cout << "  Name: " << vecStation[ii].getStationName() << endl;
      cout << "  Air Temperature: " << vecMeteo[ii].ta << endl;
    }

  Grid2DObject dem;
  ioTest->readDEM(dem);
  int nx,ny;
  ioTest->get2DGridSize(nx,ny);

  CArray2D<double> p, nswc, vw, rh, ta;
  ta.Create(nx,ny);
  p.Create(nx,ny);
  nswc.Create(nx,ny);
  rh.Create(nx,ny);
  vw.Create(nx,ny);
  
  Meteo2DInterpolator mi(dem, vecMeteo, vecStation);
  mi.interpolate(nswc, rh, ta, vw, p);

  Grid2DObject p2, nswc2, vw2, rh2, ta2;
  cout << "Convert CArray2D to Grid2DObject" << endl;
  p2.set(dem.ncols, dem.nrows, dem.xllcorner, dem.yllcorner, dem.cellsize, dem.nodata, p);
  nswc2.set(dem.ncols, dem.nrows, dem.xllcorner, dem.yllcorner, dem.cellsize, dem.nodata, nswc);
  ta2.set(dem.ncols, dem.nrows, dem.xllcorner, dem.yllcorner, dem.cellsize, dem.nodata, ta);
  rh2.set(dem.ncols, dem.nrows, dem.xllcorner, dem.yllcorner, dem.cellsize, dem.nodata, rh);
  vw2.set(dem.ncols, dem.nrows, dem.xllcorner, dem.yllcorner, dem.cellsize, dem.nodata, vw);
  cout << "conversion was successful" << endl;

  cout << "Writing the Grids to *.2d files" << endl;
  ioTest->write2DGrid(ta2, "output/ta.2d");
  ioTest->write2DGrid(p2, "output/p.2d");
  ioTest->write2DGrid(vw2, "output/vw.2d");
  ioTest->write2DGrid(nswc2, "output/nswc.2d");
  ioTest->write2DGrid(rh2, "output/rh.2d");

  cout << "Writing the Grids was successful" << endl;
}

//compile with: g++ meteoio_demo.cc ../src/Laws.c -I ../src/ -I ../src/filter/ -L ../lib/ -lmeteoIO -ldl -lm -o meteoio_demo
