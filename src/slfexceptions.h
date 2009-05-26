#ifndef __SLFEXCEPTIONS_H__
#define __SLFEXCEPTIONS_H__

#include <exception>
#include <string>
#include <iostream>

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)

/**
 * @class SLFException
 * @brief The basic exception class adjusted for the needs of SLF software
 *
 * @author Thomas Egger
 */
class SLFException : public std::exception {
 public:
  SLFException(const std::string& message="SLFException occured", const std::string& position="");
  ~SLFException() throw();
  const char* what() const throw();

 protected:
  std::string msg;
};

/**
 * @class FileNotFoundException
 * @brief thrown when a there is an unsuccessful attempt to locate a file
 *
 * @author Thomas Egger
 */
class FileNotFoundException : public SLFException {
 public:
 FileNotFoundException(const std::string& filename="",
		       const std::string& position="") : SLFException("FileNotFoundException: " + filename,position){}
};

/**
 * @class FileAccessException
 * @brief thrown when a there are insufficient rights to access the file in a certain way (e.g. read, write)
 *
 * @author Thomas Egger
 */
class FileAccessException : public SLFException {
 public:
 FileAccessException(const std::string& filename="",
		     const std::string& position="") : SLFException("FileAccessException: " + filename,position){}
};

/**
 * @class InvalidFileNameException
 * @brief thrown when a filename given is not valid (e.g. "..", "." or empty)
 *
 * @author Thomas Egger
 */
class InvalidFileNameException : public SLFException {
 public:
  InvalidFileNameException(const std::string& filename="",
			   const std::string& position="") : SLFException("InvalidFileNameException: " + filename, position){}
};

/**
 * @class InvalidFormatException
 * @brief thrown when parsed data does not reflect an expected format (e.g. premature end of a line, file)
 *
 * @author Thomas Egger
 */
class InvalidFormatException : public SLFException {
 public:
  InvalidFormatException(const std::string& message="",
			 const std::string& position="") : SLFException("InvalidFormatException: " + message, position){}
};

/**
 * @class IndexOutOfBoundsException
 * @brief thrown when an index is out of bounds
 *
 * @author Thomas Egger
 */
class IndexOutOfBoundsException : public SLFException {
 public:
  IndexOutOfBoundsException(const std::string& message="",
			    const std::string& position="") : SLFException("IndexOutOfBoundsException: " + message, position){}
};

/**
 * @class ConversionFailedException
 * @brief thrown when an unsuccessful to convert data types/classes is made (e.g. attempt to convert a literal into a number)
 *
 * @author Thomas Egger
 */
class ConversionFailedException : public SLFException {
 public:
  ConversionFailedException(const std::string& message="",
			    const std::string& position="") : SLFException("ConversionFailedException: " + message, position){}
};

/**
 * @class InvalidArgumentException
 * @brief thrown when encountered an unexpected function's argument (e.g. bad index, bad or missing parameter name, etc.)
 *
 * @author Florian Hof
 */
class InvalidArgumentException : public SLFException {
 public:
  InvalidArgumentException(const std::string& message="",
			 const std::string& position="") : SLFException("InvalidArgumentException: " + message, position){}
};

/**
 * @class UnknownValueException
 * @brief thrown when encountered an unexpected value (e.g. unknown name or key)
 *
 * @author Florian Hof
 */
class UnknownValueException : public SLFException {
 public:
  UnknownValueException(const std::string& message="",
			 const std::string& position="") : SLFException("UnknownValueException: " + message, position){}
};

/**
 * @class NoAvailableDataException
 * @brief thrown when no data is available 
 *
 * @author Florian Hof
 */
class NoAvailableDataException : public SLFException {
 public:
  NoAvailableDataException(const std::string& message="",
			 const std::string& position="") : SLFException("NoAvailableDataException: " + message, position){}
};


/// HACK for POPC by Laurent Winkler : Print error message instead of throw
// lwk debug : a THROW macro is defined
#ifdef _PAROC_
#ifdef THROW
#undef THROW
#endif
#define THROW
#ifndef __SLFEXCEPTION_CC__
#define SLFException(a,b) std::cout<<"SLFException ("<<(a)<<", "<<(b)<<") at line "<<__LINE__<<" of file "<<__FILE__<<"\n"
#define FileNotFoundException(a,b) std::cout<<"FileNotFoundException ("<<(a)<<", "<<(b)<<") at line "<<__LINE__<<" of file "<<__FILE__<<"\n"
#define FileAccessException(a,b) std::cout<<"FileAccessException ("<<(a)<<", "<<(b)<<") at line "<<__LINE__<<" of file "<<__FILE__<<"\n"
#define InvalidFileNameException(a,b) std::cout<<"InvalidFileNameException ("<<(a)<<", "<<(b)<<") at line "<<__LINE__<<" of file "<<__FILE__<<"\n"
#define InvalidFormatException(a,b) std::cout<<"InvalidFormatException ("<<(a)<<", "<<(b)<<") at line "<<__LINE__<<" of file "<<__FILE__<<"\n"
#define IndexOutOfBoundsException(a,b) std::cout<<"IndexOutOfBoundsException ("<<(a)<<", "<<(b)<<") at line "<<__LINE__<<" of file "<<__FILE__<<"\n"
#define ConversionFailedException(a,b) std::cout<<"ConversionFailedException ("<<(a)<<", "<<(b)<<") at line "<<__LINE__<<" of file "<<__FILE__<<"\n"
#define InvalidArgumentException(a,b) std::cout<<"InvalidArgumentException ("<<(a)<<", "<<(b)<<") at line "<<__LINE__<<" of file "<<__FILE__<<"\n"
#define UnknownValueException(a,b) std::cout<<"UnknownValueException ("<<(a)<<", "<<(b)<<") at line "<<__LINE__<<" of file "<<__FILE__<<"\n"
#define NoAvailableDataException(a,b) std::cout<<"NoAvailableDataException ("<<(a)<<", "<<(b)<<") at line "<<__LINE__<<" of file "<<__FILE__<<"\n"
#endif /*__SLFEXCEPTION_CC__*/
#else
#define THROW throw 
#endif /*_PAROC_*/

#endif /*__SLFEXCEPTIONS_H__*/
