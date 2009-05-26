#include "IOHandler.h"

IOHandler::IOHandler(void (*delObj)(void*)) : PluginObject(delObj){}

