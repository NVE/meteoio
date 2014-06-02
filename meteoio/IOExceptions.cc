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
	#if defined(MSG_BOX)
		#include <meteoio/MessageBoxX11.h>
	#endif
#endif
#if defined(MSWIN)
	#include <windows.h>
	#undef max
	#undef min
#endif
#if defined(APPLE) && defined(MSG_BOX)
	#include <CoreFoundation/CoreFoundation.h>
#endif

using namespace std;

namespace mio {

#if defined(MSG_BOX)
	void messageBox(const std::string& msg) {
	#if defined(LINUX) && !defined(ANDROID) && !defined(CYGWIN)
		const string box_msg = msg + "\n\nPlease check the terminal for more information!";
		MessageBoxX11("Oops, something went wrong!", box_msg.c_str());
	#endif
	#if defined(MSWIN)
		const string box_msg = msg + "\n\nPlease check the terminal for more information!";
		MessageBox( NULL, box_msg.c_str(), TEXT("Oops, something went wrong!"), MB_OK | MB_ICONERROR );
	#endif
	#if defined(APPLE)
		const string box_msg = msg + "\n\nPlease check the terminal for more information!";
		const void* keys[] = { kCFUserNotificationAlertHeaderKey,
				kCFUserNotificationAlertMessageKey };
		const void* values[] = { CFSTR("Oops, something went wrong!"),
					CFStringCreateWithCString(NULL, box_msg.c_str(), kCFStringEncodingMacRoman) };
		CFDictionaryRef dict = CFDictionaryCreate(0, keys, values,
				sizeof(keys)/sizeof(*keys), &kCFTypeDictionaryKeyCallBacks, &kCFTypeDictionaryValueCallBacks);
		SInt32 error = 0;
		CFUserNotificationCreate(NULL, 0, kCFUserNotificationStopAlertLevel, &error, dict);
	#endif
	}
#endif

#ifdef _POPC_
IOException::IOException(const std::string& message, const std::string& position) : POPException(STD_EXCEPTION)
#else
IOException::IOException(const std::string& message, const std::string& position) : msg(), full_output()
#endif
{
#if defined(MSWIN)
	const char *delim = strrchr(position.c_str(), '\\');
#else
	const char *delim = strrchr(position.c_str(), '/');
#endif
	const std::string where = (position.empty())? "unknown position" : ((delim)? delim+1 : position);
	msg = "[" + where + "] " + message;

#if defined(LINUX) && !defined(ANDROID) && !defined(CYGWIN)
	void* tracearray[25]; //maximal size for backtrace: 25 pointers
	const int tracesize = backtrace(tracearray, 25); //obtains backtrace for current thread
	char** symbols = backtrace_symbols(tracearray, tracesize); //translate pointers to strings
	std::string backtrace_info = "\n\033[01;30m**** backtrace ****\n"; //we use ASCII color codes to make the backtrace less visible/aggressive
	for (unsigned int ii=1; ii<(unsigned)tracesize; ii++) {
	#ifdef __GNUC__
		std::ostringstream ss;
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
				const std::string tmp(real_name);
				const size_t pos = tmp.find_first_of("(");
				const std::string func_name = tmp.substr(0, pos);
				const std::string func_args = tmp.substr(pos);
				ss << "\t(" << ii << ") in \033[4m" << func_name << "\033[0m\033[01;30m" << func_args << " from " << symbols[ii] << " " << offset_end << "+" << offset_begin;
			} else { // otherwise, output the mangled function name
				ss << "\t(" << ii << ") in " << mangled_name << " from " << symbols[ii] << " " << offset_end << "+" << offset_begin;
			}
			free(real_name);
		} else { // otherwise, print the whole line
			ss << "\t(" << ii << ") at " << symbols[ii];
		}
		backtrace_info += ss.str()+"\n";
	#else
		backtrace_info += "\tat " + string(symbols[ii]) + "\n";
	#endif
	}
	backtrace_info += "\033[0m"; //back to normal color
	full_output = backtrace_info + "[" + where + "] \033[31;1m" + message + "\033[0m\n";
	free(symbols);
#else
	full_output = msg+"\n";
#endif
#ifdef _POPC_
	const string tmp = backtrace_info + "\n\n" + msg;
	SetExtra(tmp.c_str());
#endif
}

IOException::~IOException() throw() {}

const char* IOException::what() const throw()
{
#if defined(MSG_BOX)
	messageBox(msg);
#endif
	return full_output.c_str();
}

} //namespace
