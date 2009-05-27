#include "NoObservedChange.h"
#include "FilterFacade.h"

// filter name
static const string c_name = string("NoObservedChange");

// parameters name
static const string c_minNbPoints = string("minNbPoints");
static const string c_minDeltaTime = string("minDeltaTime");

NoObservedChange::NoObservedChange() : 
  FilterValue1Stn()
{
  // filling parameters' list
  m_paramsName.insert(c_minNbPoints);
  m_paramsName.insert(c_minDeltaTime);

  // initialization of interpreted parameters to default values
  m_minNbPoints = 2;
  m_minDeltaTime = Date_IO(2000, 1, 1, 0, 1) - Date_IO(2000, 1, 1, 0, 0); // 1 min
}

FilterBase1Stn* createNoObservedChange() {
  return (FilterBase1Stn*)new NoObservedChange();
}

// to be called in FilterFacade::registerFilters(), otherwise the filter won't be visible!
void NoObservedChange::registerFilter() {
  FilterFacade::registerFilter(c_name, &createNoObservedChange);
}

const string NoObservedChange::getName() const {
  return c_name;
}

void NoObservedChange::getMinimalWindow(unsigned int& minNbPoints, Date_IO& minDeltaTime) {
  minNbPoints = m_minNbPoints;
  minDeltaTime = m_minDeltaTime;
}

void NoObservedChange::prepareCheck() {
  // read the ancestor's parameters
  FilterValue1Stn::prepareCheck();

  // read the parameters
  if (m_paramsValue.find(c_minNbPoints) != m_paramsValue.end()) {
    if (!IOUtils::convertString<unsigned int>(m_minNbPoints, m_paramsValue[c_minNbPoints])) 
      THROW InvalidArgumentException("parameter '"+c_minNbPoints+"' has to be an unsigned integer", AT);
  }
  if (m_paramsValue.find(c_minDeltaTime) != m_paramsValue.end()) {
    if (!IOUtils::convertString<Date_IO>(m_minDeltaTime, m_paramsValue[c_minDeltaTime])) 
      THROW InvalidArgumentException("parameter '"+c_minDeltaTime+"' has to be a time (format HH:MM)", AT);
  }
  if (m_paramsValue.find(c_minNbPoints) != m_paramsValue.end() && m_paramsValue.find(c_minDeltaTime) != m_paramsValue.end()) {
    THROW InvalidArgumentException("at least one of the 2 parameters "+c_minNbPoints+" or "+c_minDeltaTime+" has to be set", AT);
  }
}

void NoObservedChange::doCheck(MeteoBuffer& unfilteredMeteoBuffer, MeteoBuffer& filteredMeteoBuffer, unsigned int iFilteredElement) {
  stringstream tmpStringStream;
  // get the values (with a reference to it)
  MeteoData& currData = filteredMeteoBuffer.getMeteoData(iFilteredElement);
  double& currValue = getMeasureValue(currData);
  if (currValue == nodata) return;
  MeteoBufferIterator meteoIter(unfilteredMeteoBuffer, filteredMeteoBuffer, iFilteredElement);
  for (unsigned int i = 2 ; ; i++) {
    try {
      MeteoData& prevData = meteoIter.getPreviousUnfiltered();
      double& prevValue = getMeasureValue(prevData);
      if (currValue != prevValue && prevValue != nodata) break; // found change, end successfully
      if (i >= m_minNbPoints && currData.date - prevData.date >= m_minDeltaTime) {
        tmpStringStream << "measure of "<<getMeasureName()<<" at "<<currData.date.toString()<<
          " has no change over an extended period";
        if (isSoft()) {
          tmpStringStream << ", kept unchanged";
        } else {
          currValue = nodata;
          tmpStringStream << ", forced to nodata";
        }
        reportNF(tmpStringStream.str());
        break;
      }
    } catch (NoAvailableDataException& ex) {
      reportP("cannot find enough previous measure for the NoObservedChange filter, filter not applied on measure at "+currData.date.toString());
      return;
    }
  }
}

