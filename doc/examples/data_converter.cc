#include <iostream>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio

//This is the most basic example. It does not check any exceptions, it only tries to be as c-like as possible
//provide date as ISO formatted, for example 2008-12-01T15:35:00 and 
//it will retrieve the data for this date according to the io.ini configuration file
int main(int /*argc*/, char** /*argv*/) {
	std::ofstream fout;

	Date d1(2008,12,01,00,00);
	Date d2(2009,01,31,23,00);
	std::vector< std::vector<MeteoData> > vecMeteo;
	std::vector< std::vector<StationData> > vecStation;

	Config cfg("io.ini");
	IOHandler raw_io(cfg);
	BufferedIOHandler io(raw_io, cfg);

	std::cout << "Reading input data" << std::endl;

	//Very basic conversion: get the whole data set at once, with its original sampling rate
	/*io.readMeteoData(d1, d2, vecMeteo, vecStation);*/

	//More elaborate conversion: sample the data to a specific rate
	//by looping over the time and calling readMeteoData for each timestep
	std::vector<MeteoData> Meteo; //we need some intermediate storage, for storing data sets for 1 timestep
	std::vector<StationData> Station;
	io.readMeteoData(d1, Meteo, Station); //we need to know how many stations will be available
	vecMeteo.insert(vecMeteo.begin(), Station.size(), std::vector<MeteoData>()); //allocation for the vectors
	vecStation.insert(vecStation.begin(), Station.size(), std::vector<StationData>());
	for(; d1<=d2; d1+=1./24.) { //time loop, sampling rate = 1/24 day = 1 hour
		io.readMeteoData(d1, Meteo, Station); //read 1 timestep at once, forcing resampling to the timestep
		for(unsigned int ii=0; ii<Station.size(); ii++) {
			vecMeteo.at(ii).push_back(Meteo[ii]); //fill the data manually into the vector of vectors
			vecStation.at(ii).push_back(Station[ii]);
		}
	}

	//In both case, we write the data out
	std::cout << "Writing output data" << std::endl;
	io.writeMeteoData(vecMeteo, vecStation);

	std::cout << "Done!!" << std::endl;
	return 0;
}
