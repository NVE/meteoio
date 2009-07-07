#ifndef FILTERBASE_H_INCLUDED
#define FILTERBASE_H_INCLUDED

#include "IOUtils.h"
#include "MeteoData.h"
#include "StationData.h"
#include "MeteoBuffer.h"
#include "IOExceptions.h"

#include <string>
#include <vector>
#include <set>
#include <map>

using namespace std;

/**
 * @class FilterBase
 * @brief Base class (interface) for data filters. 
 *        Defines the interface for concrete filters, except the method doing the check itself, 
 *        which is defined in FilterBase1Stn and FilterBaseMStn. 
 *        This class also handle the parameters, the reporting and other things. 
 * @author Florian Hof
 * @date   2009-02-26
 */
class FilterBase {

	public:

	// base definition

	/**
	* Constructs a base filter.
	* To be called by any subclasses' constructor.
	*/
	FilterBase();

	/**
	* Destructor.
	*/
	virtual ~FilterBase(); //HACK: the destructors should be implemented in the derived classes!!

	/**
	* Returns the name of the type of the filter, used as the filter's identifier.
	* Has to be defined by the concrete filter implementation, with one unique name per type of filter.
	* Note that to be accessible, the filter also has to be registered somewhere, for example in the FilterFacade::registerFilters() method.
	*/
	virtual const string getName() const = 0;

	/**
	* Returns the minimal size of the window that the filter requires.
	* Both conditions has to be fullfilled for the filter to make the proper decision (otherwise, a warning is printed during doCheck). 
	* Before calling getMinimalWindow, the parameters have to be filled and the prepareCheck method called!
	* @param minNbPoints   [out] The minimal number of measures. Returns 1 for single-value filters (for ex. min-max).
	* @param minDeltaTime   [out] The minimal delta of the time frame. Returns an empty delta for single-value filters. 
	*/
	virtual void getMinimalWindow(unsigned int& minNbPoints, Date_IO& minDeltaTime) = 0;

	// parameters handling

	/**
	* Returns the name of the allowed filter's parameters.
	*/
	const set<string>& getParamsName() const;

	/**
	* Returns the value of the given filter's parameters.
	*/
	const string getParamValue(const string& name);

	/**
	* Sets the value of the given filter's parameters.
	* The parameter's name has to be declared (in getParamsName).
	*/
	void setParamValue(const string& name, const string& value);

	// base parameters accessors

	/**
	* Is the filter soft (warn and force the value) or hard (warn and set to no_data).
	* Filled by the preCheck method.
	*/
	bool isSoft() const;

	/**
	* Gets the name of the filter instance, usually the section name in the .ini config file.
	* @return The name of the filter instance, or empty if undefined. 
	*/
	const string getSectionName() const;

	// check handling

	/**
	* Make some initialisation and parameters' validation before the real check.
	* Has to be called before any doCheck and getMinimalWindow methods and after any parameters modification!
	*/
	virtual void prepareCheck();

	/**
	* Report a warning about a non-fulfilled condition. 
	* Used by filter implementation to report warnings.
	* @param message   The (explicit please) warning message
	*/
	void reportNF(const char* message);

	/**
	* Report a warning about a non-fulfilled condition. 
	* Used by filter implementation to report warnings.
	* @param message   The (explicit please) warning message
	*/
	void reportNF(const string& message);

	/**
	* Report a warning about a filtering problem (for example not enough data to decide). 
	* Used by filter implementation to report warnings.
	* @param message   The (explicit please) warning message
	*/
	void reportP(const char* message);

	/**
	* Report a warning about a filtering problem (for example not enough data to decide). 
	* Used by filter implementation to report warnings.
	* @param message   The (explicit please) warning message
	*/
	void reportP(const string& message);

	protected:

	// raw parameters

	/**
	* Names of the allowed filter's parameters.
	* Has to be filled in the constructor.
	*/
	set<string> m_paramsName;

	/**
	* Values of the filter's parameters.
	* The values' list is managed by this base class.
	*/
	map<string,string> m_paramsValue;

	private:

	// interpreted parameters

	/**
	* Is the filter soft (warn and force the value) or hard (warn and set to no_data).
	* Filled by the preCheck method.
	*/
	bool m_isSoft;

	// friends

	friend class FilterValue;

};

#endif // FILTERBASE_H_INCLUDED
