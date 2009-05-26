#include "FilterValue.h"

// parameters name
static const string c_measureName = string("measureName");

FilterValue::FilterValue(FilterBase* filter) : 
  m_filter(filter)
{
  // filling parameters' list
  m_filter->m_paramsName.insert(c_measureName);

  // initialization of interpreted parameters
  m_measureName = string();
  m_measurePtr = NULL;
}

void FilterValue::prepareCheck() {
  // read the ancestor's parameters
  //m_filter->prepareCheck(); // to be done where it is used

  // read the "measureName" parameter
  if (m_filter->m_paramsValue.find(c_measureName) != m_filter->m_paramsValue.end()) {
    m_measurePtr = getMeasurePtr(m_filter->m_paramsValue[c_measureName]);
    m_measureName = m_filter->m_paramsValue[c_measureName]; // assign if previous call has not failled
  } else {
    THROW InvalidArgumentException("mandatory parameter '"+c_measureName+"' not found", AT);
  }
}

const string& FilterValue::getMeasureName() const {
  return m_measureName;
}

double& FilterValue::getMeasureValue(MeteoData& data) const {
  if (! m_measurePtr) {
    THROW IOException("parameter '"+c_measureName+"' not well defined (perhaps prepareCheck has not been called?)", AT);
  }
  return data.*m_measurePtr;
}

const double& FilterValue::getMeasureValue(const MeteoData& data) const {
  if (! m_measurePtr) {
    THROW IOException("parameter '"+c_measureName+"' not well defined (perhaps prepareCheck has not been called?)", AT);
  }
  return data.*m_measurePtr;
}


double& FilterValue::getMeasureValue(MeteoData& data, double MeteoData::* measurePtr) {
  return data.*measurePtr;
}

const double& FilterValue::getMeasureValue(const MeteoData& data, const double MeteoData::* measurePtr) {
  return data.*measurePtr;
}

double MeteoData::* FilterValue::getMeasurePtr(const string& measureName) {
  if        (measureName == "ta") {
    return &MeteoData::ta;
  } else if (measureName == "iswr") {
    return &MeteoData::iswr;
  } else if (measureName == "vw") {
    return &MeteoData::vw;
  } else if (measureName == "rh") {
    return &MeteoData::rh;
  } else if (measureName == "lwr") {
    return &MeteoData::lwr;
  } else if (measureName == "nswc") {
    return &MeteoData::nswc;
  } else if (measureName == "ts0") {
    return &MeteoData::ts0;
  } else {
    THROW InvalidArgumentException("parameter '"+c_measureName+"' has unexpected value "+measureName+
        ", expected are ta, iswr, vw, rh, lwr, nswc, ts0", AT);
  }
}


FilterValue1Stn::FilterValue1Stn() : 
  FilterBase1Stn()
{
  m_filterValue = new FilterValue(this);
}

FilterValue1Stn::~FilterValue1Stn() {
  delete m_filterValue;
  m_filterValue = NULL;
}

void FilterValue1Stn::prepareCheck() {
  // read the ancestor's parameters
  FilterBase1Stn::prepareCheck();
  // read the FilterValue's parameters
  m_filterValue->prepareCheck();
}

const string& FilterValue1Stn::getMeasureName() const {
  return m_filterValue->getMeasureName();
}

double& FilterValue1Stn::getMeasureValue(MeteoData& data) const {
  return m_filterValue->getMeasureValue(data);
}

const double& FilterValue1Stn::getMeasureValue(const MeteoData& data) const {
  return m_filterValue->getMeasureValue(data);
}

FilterValueMStn::FilterValueMStn() : 
  FilterBaseMStn()
{
  m_filterValue = new FilterValue(this);
}

FilterValueMStn::~FilterValueMStn() {
  delete m_filterValue;
  m_filterValue = NULL;
}

void FilterValueMStn::prepareCheck() {
  // read the ancestor's parameters
  FilterBaseMStn::prepareCheck();
  // read the FilterValue's parameters
  m_filterValue->prepareCheck();
}

const string& FilterValueMStn::getMeasureName() const {
  return m_filterValue->getMeasureName();
}

double& FilterValueMStn::getMeasureValue(MeteoData& data) const {
  return m_filterValue->getMeasureValue(data);
}

const double& FilterValueMStn::getMeasureValue(const MeteoData& data) const {
  return m_filterValue->getMeasureValue(data);
}

