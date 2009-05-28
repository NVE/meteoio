#include "DynamicLibrary.h"

DynamicLibrary::DynamicLibrary(void* objFile) : _objFile(objFile){}

DynamicLibrary::~DynamicLibrary(void)
{
	dlclose(_objFile);
}

PluginObject* DynamicLibrary::newObject(const std::string& name, const std::string& filename)
{
  // If there is no valid library, return null
	if(_objFile == NULL) {
		return NULL;
	}

	// Get the loadObject() function.  If it doesn't exist, return NULL.
	void* loadSym = dlsym(_objFile, "loadObject");
	if(loadSym == NULL) {
		return NULL;
	}

	// Load a new instance of the requested class, and return it
	void* obj = ((void* (*)(const std::string&, const std::string&))(loadSym))(name, filename);
	return reinterpret_cast<PluginObject*>(obj);
}

DynamicLibrary* DynamicLoader::loadObjectFile(const std::string& file, int flags)
{
	void* objFile = dlopen(file.c_str(), flags);
	if(objFile == NULL) {
		return NULL;
	}

	return new DynamicLibrary(objFile);
}

PluginObject::PluginObject(void (*delObj)(void*)) : _deleteObject(delObj){}

PluginObject::~PluginObject(void){}

void PluginObject::deleteSelf(void)
{
	(*_deleteObject)(reinterpret_cast<void*>(this));
}
