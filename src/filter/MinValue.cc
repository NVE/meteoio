#include "MinValue.h"
#include "FilterFacade.h"

// filter name
static const string c_name = string("MinValue");

// parameters name
static const string c_limitValue = string("limitValue");

MinValue::MinValue() : 
  FilterValue1Stn()
{
  // filling parameters' list
  m_paramsName.insert(c_limitValue);

  // initialization of interpreted parameters
  m_limitValue = MeteoData::nodata;
}

FilterBase1Stn* createMinValue() {
  return (FilterBase1Stn*)new MinValue();
}

// to be called in FilterFacade::registerFilters(), otherwise the filter won't be visible!
void MinValue::registerFilter() {
  FilterFacade::registerFilter(c_name, &createMinValue);
}

const string MinValue::getName() const {
  return c_name;
}

void MinValue::getMinimalWindow(unsigned int& minNbPoints, Date_IO& minDeltaTime) {
  minNbPoints = 1;
  minDeltaTime = Date_IO(2000, 1, 1, 0, 0) - Date_IO(2000, 1, 1, 0, 0);
}

void MinValue::prepareCheck() {
  // read the ancestor's parameters
  FilterValue1Stn::prepareCheck();

  // read the "limitValue" parameter
  if (m_paramsValue.find(c_limitValue) != m_paramsValue.end()) {
    //m_limitValue = atof(m_paramsValue[c_limitValue].c_str());
    if (!IOUtils::convertString<double>(m_limitValue, m_paramsValue[c_limitValue]))
      THROW InvalidArgumentException("parameter '"+c_limitValue+"' has to be a float (or double)", AT);
  } else {
    THROW InvalidArgumentException("mandatory parameter '"+c_limitValue+"' not found", AT);
  }
}

void MinValue::doCheck(MeteoBuffer& unfilteredMeteoBuffer, MeteoBuffer& filteredMeteoBuffer, unsigned int iFilteredElement) {
  stringstream tmpStringStream;
  (void)unfilteredMeteoBuffer; //the compiler will know we intend not to use this parameter
  // get the value (with a reference to it)
  MeteoData& currData = filteredMeteoBuffer.getMeteoData(iFilteredElement);
  double& currValue = getMeasureValue(currData);
  // check and handle if beyond the limit
  if ((currValue != MeteoData::nodata) && (currValue < m_limitValue)) {
    tmpStringStream << "measure of "<<getMeasureName()<<" at "<<currData.date.toString()<<
      " value "<<currValue<<" is beyond the min "<<m_limitValue;
    if (isSoft()) {
      currValue = m_limitValue;
      tmpStringStream << ", forced to min";
    } else {
      currValue = MeteoData::nodata;
      tmpStringStream << ", forced to nodata";
    }
    reportNF(tmpStringStream.str());
  }
}

/*void MinValue::doCheckOne(vector<StationData>& stations, vector<vector<MeteoData> >& data, int iStation, int iData) {
  stringstream tmpStringStream;
  // get the value (with a reference to it)
  double& value = getMeasureValue(data[iStation][iData]);
  cout<<"checking measure "<<getMeasureName()<<" value "<<value<<endl;
  // check and handle if beyond the limit
  if ((value != MeteoData::nodata) && (value < m_limitValue)) {
    tmpStringStream << "measure of "<<getMeasureName()<<" at "<<value<<" is beyond the min value "<<m_limitValue;
    if (isSoft()) {
      value = m_limitValue;
      tmpStringStream << ", forced to min";
    } else {
      value = MeteoData::nodata;
      tmpStringStream << ", forced to nodata";
    }
    report(tmpStringStream.str());
  }
  cout<<"checked measure "<<getMeasureName()<<" value "<<value<<endl;
}*/

