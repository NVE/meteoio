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

}


//compile with: g++ test/meteoio_demo.cc -I ./src/ -I ./src/filter/ -L ./lib/ -lmeteoIO -ldl
