#include <iostream>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio

//This is the most basic example. It does not check any exceptions, it only tries to be as c-like as possible
//provide date as ISO formatted, for example 2008-12-01T15:35:00 and 
//it will retrieve the data for this date according to the io.ini configuration file
int main(int /*argc*/, char** argv) {
	Date d1;
	std::vector<MeteoData> vecMeteo;
	std::vector<StationData> vecStation;

	ConfigReader cfg("io.ini");
	IOHandler raw_io(cfg);
	BufferedIOHandler io(raw_io, cfg);

	IOUtils::convertString(d1,argv[1]);
	io.readMeteoData(d1, vecMeteo, vecStation);

	//writing some data out in order to prove that it really worked!
	for (unsigned int ii=0; ii < vecMeteo.size(); ii++) {
		std::cout << "---------- Station: " << (ii+1) << " / " << vecStation.size() << std::endl;
		std::cout << vecStation[ii] << std::endl;
		std::cout << vecMeteo[ii] << std::endl;
	}

	return 0;
}
