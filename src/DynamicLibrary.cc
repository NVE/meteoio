#include "DynamicLibrary.h"

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

PluginObject* DynamicLibrary::newObject(const std::string& name, const std::string& filename)
{
  // If there is no valid library, return null
	if(_objFile == NULL) {
		return NULL;
	}

	// Get the loadObject() function.  If it doesn't exist, return NULL.
#ifdef WIN32
	void (*loadSym)(const std::string&, const std::string&) = (void (*)(const std::string&, const std::string&))GetProcAddress(_objFile, "loadObject");
#else
	void* loadSym = dlsym(_objFile, "loadObject");
#endif

	if(loadSym == NULL) {
		return NULL;
	}

	// Load a new instance of the requested class, and return it
	void* obj = ((void* (*)(const std::string&, const std::string&))(loadSym))(name, filename);
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
