#ifndef NOOBSERVEDCHANGE_H_INCLUDED
#define NOOBSERVEDCHANGE_H_INCLUDED

#include "FilterBase1Stn.h"
#include "FilterValue.h"

using namespace std;

/**
 * @class NoObservedChange
 * @brief Filtering of values with no change over time (as frozen).
 *        Both the minimals number of points and time frame have to 
 *        be reached before declaring a NOC (No Observed Change). 
 * @author Florian Hof
 * @date   2009-03-16
 */
class NoObservedChange : public FilterValue1Stn {

 public:

  // base definition

  NoObservedChange();
  
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

  /** The minimal number of points (aka measures) before declaring a NOC  (default is 2) */
  unsigned int m_minNbPoints;

  /** The minimal time frame (delta time) before declaring a NOC (default is 1 minute) */
  Date_IO m_minDeltaTime;

};

#endif
