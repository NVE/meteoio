#include <iostream>
#include "MeteoIO.h"

int main() {
  Date d1(2009,01,01,18,00);

  vector<MeteoData> vecMeteo;
  vector<StationData> vecStation;

  IOHandler *ioTest=NULL; //Initialization vital!

  try {
    ioTest = new A3DIO("io.ini");
  } catch (exception& e){
    cout << "Problem with A3DIO handler cration, cause: " << e.what() << endl;
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



//compile with: g++ slfio_demo.cc -I ../snowpack/snowpack_core/ -I ../common/ -I ../snowpack/ -I filter -I . Date.o A3DIO.o ASCIIFileIO.o slfexceptions.o  slfutils.o MeteoBuffer.o MeteoData.o StationData.o ConfigReader.o DynamicLibrary.o -ldl IOHandler.o Meteo1DResampler.o libinterpol1D.o Grid2DObject.o
