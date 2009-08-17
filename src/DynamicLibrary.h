#ifndef __DYNAMICLIB_H__
#define __DYNAMICLIB_H__


#ifndef USE_PRECOMPILED_HEADERS
#ifdef WIN32
#include <direct.h>
#include <windows.h>
#else
#include <dlfcn.h>
#endif
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#endif

/**
 * @class PluginObject
 * @brief The PluginObject is an interface for all dynamically loadable Objects, its main task is to register a callback destructor function
 *        for the object dynamically allocated
 *
 * @author Thomas Egger
 * @date   2009-03-10
 */
class PluginObject {
	private:
		// Callback function that should be called to delete dynamic object
		void (*_deleteObject)(void*);
	public:
		// The constructor sets the callback function to use
		PluginObject(void (*delObj)(void*));

		// The destructor
		virtual ~PluginObject(void);

		// Sends "this" to the callback destructor function.
		void deleteSelf(void);
};

/**
 * @class DynamicLibrary
 * @brief Manages lifetime of an open dynamic library
 *
 * @author Thomas Egger
 * @date   2009-03-10
 */
class DynamicLibrary {
	protected:
		// The handle to the shared library that was opened
#ifdef WIN32
		HINSTANCE _objFile;
#else
		void *_objFile;
#endif

		// Since an instance of DynamicLibrary manages lifetime of an open 
		// library, it is important to make sure that the object isn't 
		// copied.
		DynamicLibrary(const DynamicLibrary&) {}
		DynamicLibrary& operator=(const DynamicLibrary&) {return *this;}

		// Creates a new library, with the object file handle that is passed 
		// in. Protected so that only the DynamicLoader can create an 
		// instance (since it is declared friend.
#ifdef WIN32
		DynamicLibrary(HINSTANCE objFile);
#else
		DynamicLibrary(void* objFile);
#endif

	public:
		// Destructor, closes the open shared library
		~DynamicLibrary(void);

		// Creates a new instance of the named class, or returns NULL is the 
		// class isn't found. 
		PluginObject* newObject(const std::string& name, const std::string& filename);

		friend class DynamicLoader; ///< The friend class DynamicLoader can solely instantiate the DynamicLibrary class (protected constructor)
};

/**
 * @class DynamicLoader
 * @brief The dynamic loader class, used for loading DynamicLibraries.
 *
 * @author Thomas Egger
 * @date   2009-03-10
 */
class DynamicLoader {
	public:
		static DynamicLibrary* loadObjectFile(const std::string& file);
		// Loads a DynamicLibrary, given the shared library file
		// "file", with the dlopen flags supplied.

		//Return last error message of dynamic link libraries
		static std::string getErrorMessage();
};

#endif
