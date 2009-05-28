#ifndef MAXCHANGERATE_H_INCLUDED
#define MAXCHANGERATE_H_INCLUDED

#include "FilterBase1Stn.h"
#include "FilterValue.h"
#include "IOUtils.h"

using namespace IOUtils;
using namespace std;

/**
 * @class MaxChangeRate
 * @brief Filtering of values with change bigger to a maxima (in unit/hour).
 * @author Florian Hof
 * @date   2009-03-16
 */
class MaxChangeRate : public FilterValue1Stn {

	public:

		// base definition

		MaxChangeRate();

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

		/** The maximal incremental rate, in unit/h, always positive or nodata */
		double m_maxIncrRate;

		/** The maximal decremental rate, in unit/h, always positive or nodata */
		double m_maxDecrRate;

};

#endif
