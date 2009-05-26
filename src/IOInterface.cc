#include "IOInterface.h"

IOInterface::IOInterface(void (*delObj)(void*)) : PluginObject(delObj){}

