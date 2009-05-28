#include "MinMaxValue.h"
#include "FilterFacade.h"

// filter name
static const string c_name = string("MinMaxValue");

// parameters name
static const string c_minValue = string("minValue");
static const string c_maxValue = string("maxValue");

MinMaxValue::MinMaxValue() : FilterValue1Stn()
{
	// filling parameters' list
	m_paramsName.insert(c_minValue);
	m_paramsName.insert(c_maxValue);

	// initialization of interpreted parameters
	m_minValue = nodata;
	m_maxValue = nodata;
}

FilterBase1Stn* createMinMaxValue()
{
	return (FilterBase1Stn*)new MinMaxValue();
}

// to be called in FilterFacade::registerFilters(), otherwise the filter won't be visible!
void MinMaxValue::registerFilter()
{
	FilterFacade::registerFilter(c_name, &createMinMaxValue);
}

const string MinMaxValue::getName() const
{
	return c_name;
}

void MinMaxValue::getMinimalWindow(unsigned int& minNbPoints, Date_IO& minDeltaTime)
{
	minNbPoints = 1;
	minDeltaTime = Date_IO(2000, 1, 1, 0, 0) - Date_IO(2000, 1, 1, 0, 0);
}

void MinMaxValue::prepareCheck()
{
	// read the ancestor's parameters
	FilterValue1Stn::prepareCheck();

	// read the "minValue" parameter
	if (m_paramsValue.find(c_minValue) != m_paramsValue.end()) {
		if (!IOUtils::convertString<double>(m_minValue, m_paramsValue[c_minValue])) {
			THROW InvalidArgumentException("parameter '"+c_minValue+"' has to be a float (or double)", AT);
		}
	}
	if (m_paramsValue.find(c_maxValue) != m_paramsValue.end()) {
		if (!IOUtils::convertString<double>(m_maxValue, m_paramsValue[c_maxValue])) {
			THROW InvalidArgumentException("parameter '"+c_maxValue+"' has to be a float (or double)", AT);
		}
	}
	if (m_minValue == nodata && m_maxValue == nodata) {
		THROW InvalidArgumentException("at least one of the 2 parameters "+c_minValue+" or "+c_maxValue+" has to be set", AT);
	}
}

void MinMaxValue::doCheck(MeteoBuffer& unfilteredMeteoBuffer, MeteoBuffer& filteredMeteoBuffer, unsigned int iFilteredElement)
{
	stringstream tmpStringStream;
	(void)unfilteredMeteoBuffer; //the compiler will know we intend not to use this parameter
	// get the value (with a reference to it)
	MeteoData& currData = filteredMeteoBuffer.getMeteoData(iFilteredElement);
	double& currValue = getMeasureValue(currData);
	// check and handle if beyond the limit
	if ((m_minValue != nodata) && (currValue != nodata) && (currValue < m_minValue)) {
		tmpStringStream << "measure of "<<getMeasureName()<<" at "<<currData.date.toString()<<
		" value "<<currValue<<" is beyond the min "<<m_minValue;
		if (isSoft()) {
			currValue = m_minValue;
			tmpStringStream << ", forced to min";
		} else {
			currValue = nodata;
			tmpStringStream << ", forced to nodata";
		}
		reportNF(tmpStringStream.str());
	}
	if ((m_maxValue != nodata) && (currValue != nodata) && (currValue > m_maxValue)) {
		tmpStringStream << "measure of "<<getMeasureName()<<" at "<<currData.date.toString()<<
		" value "<<currValue<<" is beyond the max "<<m_maxValue;
		if (isSoft()) {
			currValue = m_maxValue;
			tmpStringStream << ", forced to max";
		} else {
			currValue = nodata;
			tmpStringStream << ", forced to nodata";
		}
		reportNF(tmpStringStream.str());
	}
}

