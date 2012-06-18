#include <iostream>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio

//This example takes two ISO-formatted dates on the command line
//for example ./data_converter 2008-12-01T00:00:00 2008-12-31T23:00
//It will retrieve the data for this time interval and write it out as specified
//in the io.ini configuration
int main(int argc, char** argv) {
	if(argc!=2) {
		std::cout << "Invalid number of arguments! Please provide a date range!\n";
		exit(0);
	}

	Config cfg("io.ini");
	Date d1, d2;
	const double TZ = cfg.get("Input", "TIME_ZONE");
	IOUtils::convertString(d1,argv[1], TZ);
	IOUtils::convertString(d2,argv[2], TZ);

	std::vector< std::vector<MeteoData> > vecMeteo;
	IOManager io(cfg);
	//io.setProcessingLevel(IOManager::raw);
	std::cout << "Reading input data" << std::endl;

	//Very basic conversion: get the whole data set at once, with its original sampling rate
	//io.getMeteoData(d1, d2, vecMeteo);

	//More elaborate conversion: sample the data to a specific rate
	//by looping over the time and calling readMeteoData for each timestep
	std::vector<MeteoData> Meteo; //we need some intermediate storage, for storing data sets for 1 timestep
	io.getMeteoData(d1, Meteo); //we need to know how many stations will be available
	vecMeteo.insert(vecMeteo.begin(), Meteo.size(), std::vector<MeteoData>()); //allocation for the vectors
	for(; d1<=d2; d1+=1./24.) { //time loop, sampling rate = 1/24 day = 1 hour
		io.getMeteoData(d1, Meteo); //read 1 timestep at once, forcing resampling to the timestep
		for(unsigned int ii=0; ii<Meteo.size(); ii++) {
			vecMeteo.at(ii).push_back(Meteo[ii]); //fill the data manually into the vector of vectors
		}
	}

	//In both case, we write the data out
	std::cout << "Writing output data" << std::endl;
	io.writeMeteoData(vecMeteo);

	std::cout << "Done!!" << std::endl;
	return 0;
}
