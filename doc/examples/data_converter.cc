#include <iostream>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio

//This example takes two ISO-formatted dates and a sampling rate (in h) on the command line
//for example ./data_converter 2008-12-01T00:00:00 2008-12-31T23:00 1
//It will retrieve the data for this time interval and write it out once per 1 hour as specified
//in the io.ini configuration

void real_main(int argc, char** argv) {
	if(argc!=4) {
		std::cout << "Invalid number of arguments! Please provide a date range and a sampling rate (in hours)\n";
		exit(0);
	}

	Config cfg("io.ini");
	const double TZ = cfg.get("TIME_ZONE", "Input");
	Date d1, d2;
	double Tstep;
	IOUtils::convertString(d1, argv[1], TZ);
	IOUtils::convertString(d2, argv[2], TZ);
	IOUtils::convertString(Tstep, argv[3]);
	Tstep /= 24.; //convert to sampling rate in days

	IOManager io(cfg);
	std::cout << "Reading input data" << std::endl;

	Timer timer;
	timer.start();
	
	std::map<std::string, size_t> mapIDs; //over a large time range, the number of stations might change... this is the way to make it work
	std::vector<MeteoData> Meteo; //we need some intermediate storage, for storing data sets for 1 timestep
	std::vector< std::vector<MeteoData> > vecMeteo; //so we can keep and output the data that has been read
	
	for(; d1<=d2; d1+=Tstep) { //time loop
		io.getMeteoData(d1, Meteo); //read 1 timestep at once, forcing resampling to the timestep
		size_t insert_position = 0;
		for(size_t ii=0; ii<Meteo.size(); ii++) {
			const std::string stationID( Meteo[ii].meta.stationID );
			if (mapIDs.count( stationID )==0) { //if this is the first time we encounter this station, save where it should be inserted
				mapIDs[ stationID ] = insert_position++;
				vecMeteo.push_back( std::vector<MeteoData>() ); //allocating the new station
			}
			vecMeteo[ mapIDs[stationID] ].push_back(Meteo[ii]); //fill the data manually into the vector of vectors
		}
	}

	//io.getMeteoData(d1, d2, vecMeteo); //This would be the call that does NOT resample the data, instead of the above "for" loop

	timer.stop();
	//In both case, we write the data out
	std::cout << "Writing output data" << std::endl;
	io.writeMeteoData(vecMeteo);

	std::cout << "Done!! in " << timer.getElapsed() << " s" << std::endl;
}

int main(int argc, char** argv) {
	try {
		real_main(argc, argv);
	} catch(const std::exception &e) {
		std::cerr << e.what();
		exit(1);
	}
	return 0;
}
