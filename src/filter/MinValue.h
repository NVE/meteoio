#ifndef MINVALUE_H_INCLUDED
#define MINVALUE_H_INCLUDED

#include "FilterBase1Stn.h"
#include "FilterValue.h"

using namespace std;

/**
 * @class MinValue
 * @brief Filtering of values smaller than a minima.
 * @author Florian Hof
 * @date   2009-02-26
 */
class MinValue : public FilterValue1Stn {

 public:

  // base definition

  MinValue();
  
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
  //void doCheck(MeteoBuffer& unfilteredMeteoBuffer, MeteoBuffer& filteredMeteoBuffer, unsigned int iUnfilteredElement, unsigned int iFilteredElement, unsigned int iStation);
  //void doCheckOne(vector<StationData>& stations, vector<vector<MeteoData> >& data, int iStation, int iData);

 private:

  // interpreted parameters

  /** The limit of the value (here the minimal value) */
  double m_limitValue;

};

#endif 
