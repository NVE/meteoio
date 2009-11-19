#include <iostream>
#include "MeteoIO.h"

using namespace std;

int main(int argc, char** argv) {
	//provide date as ISO formatted, for example 2008-12-01T15:35:00
	//please look in the online documentation for more code examples!!
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
		ConfigReader cfg("io_demo.ini");
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
	}

	delete io;
	delete raw_io;

	return 0;
}
