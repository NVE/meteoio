#ifndef __IOPLUGIN_H__
#define __IOPLUGIN_H__

#include "DynamicLibrary.h"
#include "IOInterface.h"

class IOPlugin {
	public:
		std::string libname;
		std::string classname;
		IOInterface *io;
		DynamicLibrary *dynLibrary;
		
		IOPlugin(std::string _s1, std::string _s2, IOInterface *p1, DynamicLibrary *p2) : libname(_s1), classname(_s2), io(p1), dynLibrary(p2){}
		IOPlugin() : libname(""), classname(""), io(NULL), dynLibrary(NULL){}
};

#endif
