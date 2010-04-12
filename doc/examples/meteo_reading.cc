#include <iostream>
#include "MeteoIO.h"

//This is the most basic example. It does not check any exceptions, it only tries to be as c-like as possible
//provide date as ISO formatted, for example 2008-12-01T15:35:00 and 
//it will retrieve the data for this date according to the io.ini configuration file
int main(int argc, char** argv) {
	(void)argc;

	Date_IO d1;
	std::vector<MeteoData> vecMeteo;
	std::vector<StationData> vecStation;

	IOHandler *raw_io = NULL;
	BufferedIOHandler *io = NULL;

	ConfigReader cfg("io.ini");
	raw_io = new IOHandler(cfg);
	io = new BufferedIOHandler(*raw_io, cfg);
	
	IOUtils::convertString(d1,argv[1]);
	io->readMeteoData(d1, vecMeteo, vecStation);

	//writing some data out in order to prove that it really worked!
	for (unsigned int ii=0; ii < vecMeteo.size(); ii++) {
		std::cout << "---------- Station: " << (ii+1) << " / " << vecStation.size() << std::endl;
		std::cout << vecStation[ii] << std::endl;
		std::cout << vecMeteo[ii] << std::endl;
	}

	delete io;
	delete raw_io;

	return 0;
}
