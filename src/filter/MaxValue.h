#ifndef MAXVALUE_H_INCLUDED
#define MAXVALUE_H_INCLUDED

#include "FilterBase1Stn.h"
#include "FilterValue.h"
#include "IOUtils.h"

using namespace IOUtils;
using namespace std;

/**
 * @class MaxValue
 * @brief Filtering of values bigger than a maxima.
 * @author Florian Hof
 * @date   2009-03-16
 */
class MaxValue : public FilterValue1Stn {

	public:

		// base definition

		MaxValue();

		const string getName() const;

		void getMinimalWindow(unsigned int& minNbPoints, Date_IO& minDeltaTime);

		/** 
		* Registers the filter's definition to the FilterFacade's filter list. 
		* It is usually called by FilterFacade::registerFilters, where all filter classes are listed. 
		*/
		static void registerFilter();

		// check handling

		void prepareCheck();

		void doCheck(MeteoBuffer& unfilteredMeteoBuffer, MeteoBuffer& filteredMeteoBuffer, unsigned int iFilteredElement);

	private:

		// interpreted parameters

		/** The limit of the value (here the maximal value) */
		double m_limitValue;

};

#endif
