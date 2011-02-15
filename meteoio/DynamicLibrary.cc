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
#include <meteoio/DynamicLibrary.h>

namespace mio {

#ifdef WIN32
DynamicLibrary::DynamicLibrary(HINSTANCE objFile) : _objFile(objFile){}
#else
DynamicLibrary::DynamicLibrary(void* objFile) : _objFile(objFile){}
#endif

DynamicLibrary::~DynamicLibrary(void)
{
#ifdef WIN32
    FreeLibrary(_objFile);
#else
    dlclose(_objFile);
#endif
}

PluginObject* DynamicLibrary::newObject(const std::string& name, const Config& cfg)
{
  // If there is no valid library, return null
	if(_objFile == NULL) {
		return NULL;
	}

	// Get the loadObject() function.  If it doesn't exist, return NULL.
#ifdef WIN32
	const void (*loadSym)(const std::string&, const Config&) = (const void (*)(const std::string&, const Config&))GetProcAddress(_objFile, "loadObject");
#else
	const void* loadSym = dlsym(_objFile, "loadObject");
#endif

	if(loadSym == NULL) {
		return NULL;
	}

//HACK: this has to stay until c++ standard handles this case...
#ifdef __GNUC__
__extension__
#endif
	// Load a new instance of the requested class, and return it
	void* obj = ((void* (*)(const std::string&, const Config&))(loadSym))(name, cfg);
	return reinterpret_cast<PluginObject*>(obj);
}

DynamicLibrary* DynamicLoader::loadObjectFile(const std::string& file)
{
#ifdef WIN32
	HINSTANCE objFile = LoadLibrary(TEXT(file.c_str()));
#else
	void* objFile = dlopen(file.c_str(), RTLD_NOW);
#endif

	if(objFile == NULL) {
		return NULL;
	}

	return new DynamicLibrary(objFile);
}

std::string DynamicLoader::getErrorMessage(){
#ifdef WIN32
	std::stringstream ss;
	ss << GetLastError();
	return ss.str();
#else
	return std::string(dlerror());
#endif
}

PluginObject::PluginObject(void (*delObj)(void*)) : _deleteObject(delObj){}

PluginObject::~PluginObject(void){}

void PluginObject::deleteSelf(void)
{
	(*_deleteObject)(reinterpret_cast<void*>(this));
}

} //namespace
