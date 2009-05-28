#define __IOEXCEPTION_CC__
#include "IOExceptions.h"
using namespace std;

IOException::IOException(const std::string& message, const std::string& position)
{
	if (position=="") {
		msg = "At unknown position: " + message;
	} else {
		msg = position + ": " + message;
	}
}

IOException::~IOException() throw(){
}


const char* IOException::what() const throw()
{
	return msg.c_str();
}
#undef __IOException_CC__
