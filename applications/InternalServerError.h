#ifndef INTERNALSERVERERROR_H
#define INTERNALSERVERERROR_H

#include <exception>
#include <string>

using namespace std;

class InternalServerError : public exception
{
public:
    InternalServerError(const string &msg) : message(msg) {}

    virtual const char *what() const throw()
    {
        return message.c_str();
    }
private:
    const string message;
};

#endif // INTERNALSERVERERROR_H