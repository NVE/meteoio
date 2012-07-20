#include <stdio.h>
#include <stdlib.h>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio
using namespace std;

// ----- Constants ----
// HACK create  this durations... better way to construct them ??
const Duration minute = mio::Date(2008, 12, 01, 3, 36, 00, 0) - mio::Date(2008, 12, 01, 3, 35, 00, 0);
const Duration day = mio::Date(2008, 12, 02, 0, 00, 00, 0) - mio::Date(2008, 12, 01, 0, 00, 00, 0);
const int minutes_per_day = 60 * 24;

mio::SunObject Sun(46.77181, 9.86820, 2192.); //Stillberg station
const double slope_azi=38., slope_elev=35.;   //Stillberg station
const double TZ=1.;
const double TA = 273.15+11., RH = 0.5, mean_albedo = 0.5;

string f_reference("reverence_results.txt");
string f_output("generated_results.txt");

// ---- INPUT Data to controll----
// if you add an date here, don't forget to add it to list with push_pack in main !!!!
const mio::Date d1(2000, 12, 16, 01, 00, 00, 1); // Winter, good Weather
const mio::Date d2(2000, 12, 26, 01, 00, 00, 1); // Winter, "bad" Weather
const mio::Date d3(2001, 05, 29, 01, 00, 00, 1); // Summer, good Weather
const mio::Date d4(2001, 05, 26, 01, 00, 00, 1); // Summer, "bad" Weather
const mio::Date d5(2001, 04, 14, 01, 00, 00, 1); // Higher income then possible
const mio::Date d6(2001, 05, 22, 01, 00, 00, 1); // double peaks

const double iswr_ref []= {-1000., -100., -10., -1., -0.1, -0.01, 0, 0.01, 0.1, 1., 10., 100., 1000};

// write out ref file
bool writeSun24h(ofstream& os, const mio::Date start_date, const double iswr_ref) {

	for(mio::Date date(start_date);date <= (start_date+day); date=date+minute){
		os << date.toString(Date::ISO) << "\t";
		os << std::setprecision(10);
		
		Sun.setDate(date.getJulianDate(), date.getTimeZone()); //local julian date and timezone
		
		Sun.calculateRadiation(TA, RH, mean_albedo);

		//radiation in the beam'
		double b_toa, b_direct, b_diffuse, md_beam;
		Sun.getBeamRadiation(b_toa, b_direct, b_diffuse);
		md_beam=Sun.getSplitting(b_toa,iswr_ref);
		
		os << b_toa << "\t" << b_direct << "\t" << b_diffuse << "\t";
		os << md_beam << "\t";

		//radiation on the horizontal
		double h_toa, h_direct, h_diffuse, md_horizontal;
		Sun.getHorizontalRadiation(h_toa, h_direct, h_diffuse);
		md_horizontal=Sun.getSplitting(h_toa,iswr_ref);
		
		os << h_toa << "\t" << h_direct << "\t" << h_diffuse << "\t";
		os << md_horizontal << "\t";

		//radiation on the slope
		double s_toa, s_direct, s_diffuse, md_slope;
		Sun.getSlopeRadiation(slope_azi, slope_elev, s_toa, s_direct, s_diffuse);
		md_slope=Sun.getSplitting(s_toa,iswr_ref);
		
		os << s_toa << "\t" << s_direct << "\t" << s_diffuse << "\t";
		os << md_slope << "\t";

		//other sun stuff
		double solar_azi, solar_elev, eccentricity;
		double sunrise, sunset, daylight;
		Sun.position.getHorizontalCoordinates(solar_azi, solar_elev, eccentricity);
		Sun.position.getDaylight(sunrise, sunset, daylight, date.getTimeZone());
		
		os << Sun.getElevationThresh() << "\t";
		os << Sun.position.getAngleOfIncidence(slope_azi,slope_elev) << "\t";
		os << eccentricity << "\t";
		//os << Hour Angle << "\t";
		os << solar_azi << "\t" << solar_elev << "\t";
		os << IOUtils::printFractionalDay(sunrise)<< "\t";
		os << IOUtils::printFractionalDay(sunset)<< "\t";
		os << IOUtils::printFractionalDay(daylight/minutes_per_day)<< "\t";
		
		// end line
		os << endl;
	}
	return true;
}


// print out header to know which line is which
void printHeader(ostream& os){
	os << "# date" << "\t" ;
	os << "r. beam total" << "\t"<< "r. beam direct" << "\t"<< "r. beam diffues" << "\t"<< "r. beam splitting" << "\t";
	os << "r. horizontal total" << "\t"<< "r. horizontal direct" << "\t"<< "r. horizontal diffuse" << "\t"<< "r. horizontal splitting" << "\t";
	os << "r. slope total" << "\t"<< "r. slope direct" << "\t"<< "r. slope diffuse" << "\t"<< "r. slope splitting" << "\t";
	os << "Elevation Thresh" << "\t"<< "Eccentricity " << "\t"<< "solar Azimut" << "\t"<< "solar elevation" << "\t"<< "Sunrise Time" << "\t"<< "Sunset Time" << "\t"<< "Sunglight" << "\t";
	os << endl;
}

//Test if sun simulation at Stilberg station at different dates and different iswr_ref
int main() {
	
	// ----- Cenerate list for loops --------
	cout << " --- Init Variables \n";
	
	list<mio::Date> date;
	date.push_back(d1);
	date.push_back(d2);
	date.push_back(d3);
	date.push_back(d4);
	date.push_back(d5);
	date.push_back(d6);

	list<double> iswr(iswr_ref, iswr_ref + sizeof(iswr_ref) / sizeof(double));
	
	// ----- Write reference file ------
	cout << " --- Start writing Output file : \n";
	
	ofstream ofs(f_output.c_str(), ofstream::out | ofstream::trunc);
	
	printHeader(ofs);
	for (list<mio::Date>::iterator it_date = date.begin(); it_date != date.end(); it_date++) {
		cout << " -- read information for date : " << (*it_date).toString(Date::ISO) << endl;
		for (list<double>::iterator it_iswr = iswr.begin(); it_iswr != iswr.end(); it_iswr++) {
			cout << " - read information for reference iswr : " << *it_iswr << endl;
			if(!writeSun24h(ofs, *it_date, *it_iswr)){
				exit(1);
			}
		}
	}
	ofs.close();
	
	// ------ Compare reference file with generated results ---------
	cout << " --- Start comparing reference file with output file generated before \n";	
	ifstream ifref(f_reference.c_str());
	ifstream ifout(f_output.c_str());
	string l_ref, l_out;
	
	while (!ifref.eof())
	{
		if(ifout.eof()){
			cerr << "Not enough lines generated as result!!!" << endl;
			exit(1);
		}
		
		getline(ifref,l_ref);
		getline(ifout,l_out);
		if (l_ref!=l_out) {
			cerr << " ERROR, Sun generatet error at following point an error " << endl;
			printHeader(cerr);
			cerr << "ref : \n " << l_ref << endl;
			cerr << "out : \n " << l_out << endl;
			exit(1);
		}
	}
	
	if(!ifout.eof()){
		cerr << "To much lines generated as result!!!" << endl;
		exit(1);
	}
	ifout.close(); 
	ifref.close();
	
	return 0;
}
