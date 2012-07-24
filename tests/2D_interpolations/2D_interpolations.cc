#include <iostream>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio
using namespace std;


// HACK SPEACK WITH Mathias if use ndiff-2 !!! (Problem no-data -999, relativ error, do it here ???)

//This is the a basic example of spatial interpolations. It does not check any exceptions, it only tries to be as c-like as possible
//provide date as ISO formatted, for example 2008-12-01T15:35:00 and
//it will retrieve and interpolate the data for this date according to the io.ini configuration file
int main(int /*argc*/, char** argv) {
	Date d1;

	//initializing the io handlers according to the config file
	Config cfg("io.ini");
	IOManager io(cfg);

	//reading the dem (necessary for several spatial interpolations algoritms)
	DEMObject dem;
	io.readDEM(dem);

	//we assume that the time given on the command line is in TZ=+1
	IOUtils::convertString(d1,argv[1], 1.);

	//performing spatial interpolations
	Grid2DObject param;
	io.interpolate(d1, dem, MeteoData::TA, param);
	io.write2DGrid(param, MeteoGrids::TA, d1);
	io.interpolate(d1, dem, MeteoData::HNW, param);
	io.write2DGrid(param, MeteoGrids::HNW, d1);
	io.interpolate(d1, dem, MeteoData::RH, param);
	io.write2DGrid(param, MeteoGrids::RH, d1);
	io.interpolate(d1, dem, MeteoData::RSWR, param);
	io.write2DGrid(param, MeteoGrids::RSWR, d1);

	string search_string = "-999.000";
	string replace_string = "0";
	string inbuf;
	ifstream input_file("2009-01-19T12.00_HNW.asc");
	ofstream output_file("first.txt");

	while (!input_file.eof()){

		size_t spot=0;
		size_t hiho = 0;

		getline(input_file, inbuf);

		do {

			spot = inbuf.find(search_string);
			if(spot >= 0){
				string tmpstring = inbuf.substr(0,spot);
				tmpstring += replace_string;
				tmpstring += inbuf.substr(spot+search_string.length(), inbuf.length());
				inbuf = tmpstring;
			}

			cout << hiho++<< endl;

		} while (spot != string::npos);

		output_file << inbuf << endl;
	}


	//ifstream ifout("");
	return 0;
}
