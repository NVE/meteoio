#include <iostream>
#include <meteoio/MeteoIO.h>

using namespace mio;

//This is the a basic example of time manipulation
//provide date as ISO formatted, for example 2008-12-01T15:35:00
int main(int argc, char** argv) {
	(void)argc;
	const double TZ=-4.;

	Date tmp;
	tmp.setTimeZone(TZ); //HACK
	tmp.setFromSys();
	std::cout << "now=" << tmp;
	double julian=tmp.getJulianDate();
	julian = floor(julian*24.*2. + 0.5) / (24.*2.);
	tmp.setDate(julian, TZ);
	std::cout << "Rounded now=" << tmp;
	return 0;

	Date d1;
	IOUtils::convertString(d1,argv[1], 0);
	std::cout << "In timezone GMT+0:\n";
	std::cout << d1 << std::endl;

	std::cout << "In timezone GMT " << TZ << ":\n";
	d1.setTimeZone(TZ,false);
	std::cout << d1 << std::endl;

	std::cout << "Same, directly read in timezone GMT " << TZ << ":\n";
	d1.setTimeZone(TZ,false);
	IOUtils::convertString(d1,argv[1], TZ);
	std::cout << d1 << std::endl;

	std::cout << "And swapped back to timezone GMT+0:\n";
	d1.setTimeZone(0.,false);
	std::cout << d1 << std::endl;

	std::cout << "Dates comparison in same TZ:\n";
	d1.setTimeZone(TZ,false);
	Date d2(d1);
	d2-=1./24.;
	d2.setTimeZone(0.,false);
	std::cout << (d1<=d2) << std::endl;
	
	return 0;
}
