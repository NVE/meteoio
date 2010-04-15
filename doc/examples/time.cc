#include <iostream>
#include "MeteoIO.h"

//This is the a basic example of time manipulation
//provide date as ISO formatted, for example 2008-12-01T15:35:00
int main(int argc, char** argv) {
	(void)argc;

	Date_IO d1;
	IOUtils::convertString(d1,argv[1]);
	std::cout << "In timezone GMT+0:\n";
	std::cout << d1 << std::endl;

	std::cout << "In timezone GMT+1:\n";
	d1.setTimeZone(+1.,false);
	std::cout << d1 << std::endl;

	std::cout << "Same, directly read in timezone GMT+1:\n";
	d1.setTimeZone(+1.,false);
	IOUtils::convertString(d1,argv[1]);
	std::cout << d1 << std::endl;

	std::cout << "And swapped back to timezone GMT+0:\n";
	d1.setTimeZone(0.,false);
	std::cout << d1 << std::endl;
	
	return 0;
}
