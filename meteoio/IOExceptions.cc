/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <meteoio/IOExceptions.h>

#include <string.h>
#if defined(LINUX) && !defined(ANDROID) && !defined(CYGWIN)
	#include <execinfo.h> //needed for the backtracing of the stack
	#if defined(__GNUC__)
		#include <sstream>
		#include <cxxabi.h>
	#endif
#endif
#if defined(WIN32)
    #include <windows.h>
#endif

using namespace std;

namespace mio {

#ifdef _POPC_
IOException::IOException(const std::string& message, const std::string& position) : POPException(STD_EXCEPTION)
#else
IOException::IOException(const std::string& message, const std::string& position) : msg()
#endif
{
	if (position=="") {
		msg = "At unknown position: " + message;
	} else {
		msg = "[" + (strrchr(position.c_str(), '/') ? strrchr(position.c_str(), '/') + 1 : position) + "] " + message;
		//msg = position + ": " + message;
	}
#if defined(LINUX) && !defined(ANDROID) && !defined(CYGWIN)
	void* tracearray[25]; //maximal size for backtrace: 25 pointers
	size_t tracesize = backtrace(tracearray, 25); //obtains backtrace for current thread
	char** symbols = backtrace_symbols(tracearray, tracesize); //translate pointers to strings
	msg += "\n\n\033[01;30m**** backtrace ****\n"; //we use ASCII color codes to make the backtrace less visible/aggressive
	for (unsigned int ii=1; ii<tracesize; ii++) {
	#ifdef __GNUC__
		std::stringstream ss;
		char *mangled_name = 0, *offset_begin = 0, *offset_end = 0;
		for (char *p = symbols[ii]; *p; ++p) {
			// find parantheses and +address offset surrounding mangled name
			if (*p == '(') mangled_name = p;
			else if (*p == '+') offset_begin = p;
			else if (*p == ')') offset_end = p;
		}
		if (mangled_name && offset_begin && offset_end && mangled_name < offset_begin) {
			//the line could be processed, attempt to demangle the symbol
			*mangled_name++ = '\0'; *offset_begin++ = '\0'; *offset_end++ = '\0';
			int status;
			char *real_name = abi::__cxa_demangle(mangled_name, 0, 0, &status);
			// if demangling is successful, output the demangled function name
			if (status == 0) {
				ss << "\t(" << ii << ") in " << real_name << " from " << symbols[ii] << " " << offset_end << "+" << offset_begin;
			} else { // otherwise, output the mangled function name
				ss << "\t(" << ii << ") in " << mangled_name << " from " << symbols[ii] << " " << offset_end << "+" << offset_begin;
			}
			free(real_name);
		} else { // otherwise, print the whole line
			ss << "\t(" << ii << ") at " << symbols[ii];
		}
		msg += ss.str()+"\n";
	#else
		msg += "\tat " + string(symbols[ii]) + "\n";
	#endif
	}
	msg += "\033[0m"; //back to normal color
	free(symbols);
#endif
#ifdef _POPC_
	//printf("IOException(%d): %s\n",Code(),msg.c_str());
	SetExtra(msg.c_str());
#endif
}

IOException::~IOException() throw(){
}

const char* IOException::what() const throw()
{
#if defined(WIN32)
    const string tmp = msg + "\n\nPlease check the terminal for more information!";
    MessageBox ( NULL, tmp.c_str(), TEXT("Oops, something went wrong!"), MB_OK | MB_ICONERROR );
#endif
	return msg.c_str();
}

} //namespace
