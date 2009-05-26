#ifndef MINMAXVALUE_H_INCLUDED
#define MINMAXVALUE_H_INCLUDED

#include "FilterBase1Stn.h"
#include "FilterValue.h"

using namespace std;

/**
 * @class MinMaxValue
 * @brief Filtering of values smaller than a minima and bigger than a maxima.
 * @author Florian Hof
 * @date   2009-03-19
 */
class MinMaxValue : public FilterValue1Stn {

 public:

  // base definition

  MinMaxValue();
  
  const string getName() const;

  void getMinimalWindow(unsigned int& minNbPoints, Date& minDeltaTime);

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

  /** The minimal value */
  double m_minValue;

  /** The maximal value */
  double m_maxValue;

};

#endif
