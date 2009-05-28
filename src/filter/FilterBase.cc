#include "FilterBase.h"

// parameters name
static const string c_isSoft = "isSoft";
static const string c_sectionName = "sectionName";

FilterBase::FilterBase()
{
	// filling parameters' list
	m_paramsName.insert(c_isSoft);
	m_paramsName.insert(c_sectionName);

	// initialization of interpreted parameters
	m_isSoft = false;
	m_paramsValue[c_sectionName] = string();
}

FilterBase::~FilterBase() {
}

const set<string>& FilterBase::getParamsName() const
{
	return m_paramsName;
}

const string FilterBase::getParamValue(const string& name)
{
	return m_paramsValue[name];
}

void FilterBase::setParamValue(const string& name, const string& value)
{
	if (m_paramsName.find(name) == m_paramsName.end()) {
		THROW InvalidArgumentException("invalid parameter named " + name, AT);
	} else {
		m_paramsValue[name] = value;
	}
}

void FilterBase::prepareCheck()
{
	// read the "isSoft" parameter
	if (m_paramsValue.find(c_isSoft) != m_paramsValue.end()) {
		if (!IOUtils::convertString<bool>(m_isSoft, m_paramsValue[c_isSoft])) {
			THROW InvalidArgumentException("parameter '"+c_isSoft+"' has to be a boolean",AT);
		}
	}
}

void FilterBase::reportNF(const char* message)
{
	printf("[W] NON-FULFILLED FILTER %s: %s\n", getSectionName().c_str(), message);
}

void FilterBase::reportNF(const string& message)
{
	reportNF(message.c_str());
}

void FilterBase::reportP(const char* message)
{
	printf("[W] PROBLEM ON FILTER %s: %s\n", getSectionName().c_str(), message);
}

void FilterBase::reportP(const string& message)
{
	reportP(message.c_str());
}

bool FilterBase::isSoft() const
{
	return m_isSoft;
}

const string FilterBase::getSectionName() const
{
	return (const_cast<FilterBase*>(this))->m_paramsValue[c_sectionName];
}
