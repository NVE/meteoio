#define __SLFEXCEPTION_CC__
#include "slfexceptions.h"
using namespace std;

SLFException::SLFException(const std::string& message, const std::string& position){
  if (position=="")
    msg = "At unknown position: " + message;
  else 
    msg = position + ": " + message;
}

SLFException::~SLFException() throw(){
}


const char* SLFException::what() const throw(){
  return msg.c_str();
}
#undef __SLFEXCEPTION_CC__
