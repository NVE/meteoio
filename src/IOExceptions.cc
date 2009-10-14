#include "IOExceptions.h"
using namespace std;

#ifdef _POPC_
IOException::IOException(const std::string& message, const std::string& position):POPException(STD_EXCEPTION)
#else
IOException::IOException(const std::string& message, const std::string& position)
#endif
{
	if (position=="") {
		msg = "At unknown position: " + message;
	} else {
		msg = position + ": " + message;
	}
#ifdef _POPC_
//	printf("IOException(%d): %s\n",Code(),msg.c_str());
	SetExtra(msg.c_str());
#endif
}

IOException::~IOException() throw(){
}


const char* IOException::what() const throw()
{
	return msg.c_str();
}
