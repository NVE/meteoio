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

#ifdef _WIN32
DynamicLibrary::DynamicLibrary(HINSTANCE i_objFile) : objFile(i_objFile){}
#else
DynamicLibrary::DynamicLibrary(void* i_objFile) : objFile(i_objFile){}
#endif

DynamicLibrary::~DynamicLibrary(void)
{
#ifdef _WIN32
    FreeLibrary(objFile);
#else
    dlclose(objFile);
#endif
}

PluginObject* DynamicLibrary::newObject(const std::string& name, const Config& cfg)
{
  // If there is no valid library, return null
	if(objFile == NULL) {
		return NULL;
	}

	// Get the loadObject() function.  If it doesn't exist, return NULL.
#ifdef _WIN32
	#pragma warning(disable:4191) //GetProcAddress does NOT return a FARPROC, the warning misses it...
	const void (*loadSym)(const std::string&, const Config&) = (const void (*)(const std::string&, const Config&))GetProcAddress(objFile, "loadObject");
#else
	const void* loadSym = dlsym(objFile, "loadObject");
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
#ifdef _WIN32
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
#ifdef _WIN32
	std::stringstream ss;
	ss << GetLastError();
	return ss.str();
#else
	return std::string(dlerror());
#endif
}

PluginObject::PluginObject(void (*i_delObj)(void*)) : deleteObject(i_delObj){}

PluginObject::~PluginObject(void){}

void PluginObject::deleteSelf(void)
{
	(*deleteObject)(reinterpret_cast<void*>(this));
}

} //namespace
