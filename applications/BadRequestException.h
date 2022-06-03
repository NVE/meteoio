#ifndef BADREQUESTEXCEPTION_H
#define BADREQUESTEXCEPTION_H

#include <exception>
#include <string>

using namespace std;

class BadRequestException : public exception
{
public:
    BadRequestException(const string &msg) : message(msg) {}

    virtual const char *what() const throw()
    {
        return message.c_str();
    }
private:
    const string message;
};

#endif // BADREQUESTEXCEPTION_H