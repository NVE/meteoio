#include <iostream>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio

//This is the most basic example. It does not check any exceptions, it only tries to be as c-like as possible
//provide date as ISO formatted, for example 2008-12-01T15:35:00 and 
//it will retrieve the data for this date according to the io.ini configuration file
int main(int /*argc*/, char** /*argv*/) {
	std::ofstream fout;

	Date d1(2010,07,1,00,00);
	Date d2(2010,10,1,23,00);
	std::vector< std::vector<MeteoData> > vecMeteo;

	Config cfg("io.ini");
	IOManager io(cfg);

	std::cout << "Reading input data" << std::endl;

	//Very basic conversion: get the whole data set at once, with its original sampling rate
	/*io.readMeteoData(d1, d2, vecMeteo, vecStation);*/

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
