#include "MaxChangeRate.h"
#include "FilterFacade.h"

// filter name
static const string c_name = string("MaxChangeRate");

// parameters name
static const string c_maxRate = string("maxRate");
static const string c_maxIncrRate = string("maxIncrRate");
static const string c_maxDecrRate = string("maxDecrRate");

MaxChangeRate::MaxChangeRate() : 
  FilterValue1Stn()
{
  // filling parameters' list
  m_paramsName.insert(c_maxRate);
  m_paramsName.insert(c_maxIncrRate);
  m_paramsName.insert(c_maxDecrRate);

  // initialization of interpreted parameters
  m_maxIncrRate = MeteoData::nodata;
  m_maxDecrRate = MeteoData::nodata;
}

FilterBase1Stn* createMaxChangeRate() {
  return (FilterBase1Stn*)new MaxChangeRate();
}

// to be called in FilterFacade::registerFilters(), otherwise the filter won't be visible!
void MaxChangeRate::registerFilter() {
  FilterFacade::registerFilter(c_name, &createMaxChangeRate);
}

const string MaxChangeRate::getName() const {
  return c_name;
}

void MaxChangeRate::getMinimalWindow(unsigned int& minNbPoints, Date& minDeltaTime) {
  minNbPoints = 2;
  minDeltaTime = Date(2000, 1, 1, 0, 0) - Date(2000, 1, 1, 0, 0);
}

void MaxChangeRate::prepareCheck() {
  // read the ancestor's parameters
  FilterValue1Stn::prepareCheck();

  // read the parameters
  if (m_paramsValue.find(c_maxRate) != m_paramsValue.end()) {
    if (m_paramsValue.find(c_maxIncrRate) == m_paramsValue.end() && m_paramsValue.find(c_maxDecrRate) == m_paramsValue.end()) {
      if (!slfutils::convertString<double>(m_maxIncrRate, m_paramsValue[c_maxRate])) 
        THROW InvalidArgumentException("parameter '"+c_maxRate+"' has to be a float (or double)", AT);
      m_maxDecrRate = m_maxIncrRate;
    } else {
      THROW InvalidArgumentException("when parameter '"+c_maxRate+"' is set, neither "+c_maxIncrRate+" nor "+c_maxDecrRate+" should be set", AT);
    }
  } else {
    if (m_paramsValue.find(c_maxIncrRate) != m_paramsValue.end()) {
      if (!slfutils::convertString<double>(m_maxIncrRate, m_paramsValue[c_maxIncrRate])) 
        THROW InvalidArgumentException("parameter '"+c_maxIncrRate+"' has to be a float (or double)", AT);
    }
    if (m_paramsValue.find(c_maxDecrRate) != m_paramsValue.end()) {
      if (!slfutils::convertString<double>(m_maxDecrRate, m_paramsValue[c_maxDecrRate])) 
        THROW InvalidArgumentException("parameter '"+c_maxDecrRate+"' has to be a float (or double)", AT);
    }
    if (m_maxIncrRate == MeteoData::nodata && m_maxDecrRate == MeteoData::nodata) {
      THROW InvalidArgumentException("at least one of the 3 parameters "+c_maxRate+", "+c_maxIncrRate+" or "+c_maxDecrRate+" has to be set", AT);
    }
  }
}

void MaxChangeRate::doCheck(MeteoBuffer& unfilteredMeteoBuffer, MeteoBuffer& filteredMeteoBuffer, unsigned int iFilteredElement) {
  stringstream tmpStringStream;
  // get the values (with a reference to it)
    MeteoData& currData = filteredMeteoBuffer.getMeteoData(iFilteredElement);
    double& currValue = getMeasureValue(currData);
    MeteoBufferIterator meteoIter(unfilteredMeteoBuffer, filteredMeteoBuffer, iFilteredElement);
  try {
    MeteoData& prevData = meteoIter.getPrevious();
    double& prevValue = getMeasureValue(prevData);
    // check and handle if beyond the limit
    if ((currValue != MeteoData::nodata) && (prevValue != MeteoData::nodata)) {
      if (currValue > prevValue && m_maxIncrRate != MeteoData::nodata) {
        double deltaHour = (currData.date.getJulian() - prevData.date.getJulian()) * 24.;
        if ( (currValue - prevValue) / deltaHour > m_maxIncrRate ) {
          tmpStringStream << "measure of "<<getMeasureName()<<" at "<<currData.date.toString()<<
            " changing to fast (from "<<prevValue<<" up to "<<currValue<<" in "<<deltaHour<<" hours)";
          if (isSoft()) {
            tmpStringStream << ", kept unchanged";
          } else {
            currValue = MeteoData::nodata;
            tmpStringStream << ", forced to nodata";
          }
          reportNF(tmpStringStream.str());
        } else {}
      } else if (currValue < prevValue && m_maxDecrRate != MeteoData::nodata) {
        double deltaHour = (currData.date.getJulian() - prevData.date.getJulian()) * 24.;
        if ( (prevValue - currValue) / deltaHour > m_maxDecrRate ) {
          tmpStringStream << "measure of "<<getMeasureName()<<" at "<<currData.date.toString()<<
            " changing to fast (from "<<prevValue<<" down to "<<currValue<<" in "<<deltaHour<<" hours)";
          if (isSoft()) {
            tmpStringStream << ", kept unchanged";
          } else {
            currValue = MeteoData::nodata;
            tmpStringStream << ", forced to nodata";
          }
          reportNF(tmpStringStream.str());
        } else {}
      }
    }
  } catch (NoAvailableDataException& ex) {
    tmpStringStream<<"cannot find 1 previous measure for the MaxChangeRate filter, filter not applied on measure at "<<currData.date.toString();
    reportP(tmpStringStream.str());
    return;
  }
}

