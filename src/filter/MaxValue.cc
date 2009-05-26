#include "MaxValue.h"
#include "FilterFacade.h"

// filter name
static const string c_name = string("MaxValue");

// parameters name
static const string c_limitValue = string("limitValue");

MaxValue::MaxValue() : 
  FilterValue1Stn()
{
  // filling parameters' list
  m_paramsName.insert(c_limitValue);

  // initialization of interpreted parameters
  m_limitValue = MeteoData::nodata;
}

FilterBase1Stn* createMaxValue() {
  return (FilterBase1Stn*)new MaxValue();
}

// to be called in FilterFacade::registerFilters(), otherwise the filter won't be visible!
void MaxValue::registerFilter() {
  FilterFacade::registerFilter(c_name, &createMaxValue);
}

const string MaxValue::getName() const {
  return c_name;
}

void MaxValue::getMinimalWindow(unsigned int& minNbPoints, Date_IO& minDeltaTime) {
  minNbPoints = 1;
  minDeltaTime = Date_IO(2000, 1, 1, 0, 0) - Date_IO(2000, 1, 1, 0, 0);
}

void MaxValue::prepareCheck() {
  // read the ancestor's parameters
  FilterValue1Stn::prepareCheck();

  // read the "limitValue" parameter
  if (m_paramsValue.find(c_limitValue) != m_paramsValue.end()) {
    if (!IOUtils::convertString<double>(m_limitValue, m_paramsValue[c_limitValue]))
      THROW InvalidArgumentException("parameter '"+c_limitValue+"' has to be a float (or double)", AT);
  } else {
    THROW InvalidArgumentException("mandatory parameter '"+c_limitValue+"' not found", AT);
  }
}

void MaxValue::doCheck(MeteoBuffer& unfilteredMeteoBuffer, MeteoBuffer& filteredMeteoBuffer, unsigned int iFilteredElement) {
  stringstream tmpStringStream;
  (void)unfilteredMeteoBuffer; //the compiler will know we intend not to use this parameter
  // get the value (with a reference to it)
  MeteoData& currData = filteredMeteoBuffer.getMeteoData(iFilteredElement);
  double& currValue = getMeasureValue(currData);
  // check and handle if beyond the limit
  if ((currValue != MeteoData::nodata) && (currValue > m_limitValue)) {
    tmpStringStream << "measure of "<<getMeasureName()<<" at "<<currData.date.toString()<<
      " value "<<currValue<<" is beyond the max "<<m_limitValue;
    if (isSoft()) {
      currValue = m_limitValue;
      tmpStringStream << ", forced to max";
    } else {
      currValue = MeteoData::nodata;
      tmpStringStream << ", forced to nodata";
    }
    reportNF(tmpStringStream.str());
  }
}

