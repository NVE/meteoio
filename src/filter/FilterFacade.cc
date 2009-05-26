#include "FilterFacade.h"
#include "IOHandler.h"
#include "ConfigReader.h"

// references to each filters
#include "MinValue.h"
#include "MaxValue.h"
#include "MinMaxValue.h"
#include "MaxChangeRate.h"
#include "NoObservedChange.h"

using namespace std;

static const string c_name = "name";
static const string c_sectionName = "sectionName";

map<string, Filter1StnCreator> FilterFacade::s_registeredFilters1StnByName;
map<string, FilterMStnCreator> FilterFacade::s_registeredFiltersMStnByName;
static bool initialized = false;

FilterFacade::FilterFacade() {
  if (!initialized) registerAllFilters();
  m_isEasyMode = false;
}

FilterFacade::FilterFacade(const string& filename) {
  if (!initialized) registerAllFilters();
  m_isEasyMode = true;
  addActiveFilters(filename);
  prepareCheck();
}

FilterFacade::~FilterFacade() {
  if (m_isEasyMode) {
    clearActiveFilters();
  }
}

void FilterFacade::registerAllFilters() {
  if (!initialized) {
    FilterFacade::s_registeredFilters1StnByName = map<string, Filter1StnCreator>();
    FilterFacade::s_registeredFiltersMStnByName = map<string, FilterMStnCreator>();
    initialized = true;
    
    // this is the list of known filters; add your filter to the list below  **********
    MinValue::registerFilter();
    MaxValue::registerFilter();
    MinMaxValue::registerFilter();
    MaxChangeRate::registerFilter();
    NoObservedChange::registerFilter();
  }
}

void FilterFacade::registerFilter(const string& name, const Filter1StnCreator creator) {
  if (!initialized) registerAllFilters();
  stringstream tmpStringStream;
  if (name.empty()) {
    tmpStringStream<<"Invalid empty name for filter.";
    THROW InvalidArgumentException(tmpStringStream.str(), AT);
  }
  if (!creator) {
    tmpStringStream<<"Invalid null creator for filter '"<<name<<"'.";
    THROW InvalidArgumentException(tmpStringStream.str(), AT);
  }
  if (FilterFacade::s_registeredFilters1StnByName.find(name) != FilterFacade::s_registeredFilters1StnByName.end()) {
    tmpStringStream<<"A filter already exists with name '"<<name<<"'.";
    THROW InvalidArgumentException(tmpStringStream.str(), AT);
  }
  if (FilterFacade::s_registeredFiltersMStnByName.find(name) != FilterFacade::s_registeredFiltersMStnByName.end()) {
    tmpStringStream<<"A multi-station filter already exists with name '"<<name<<"'.";
    THROW InvalidArgumentException(tmpStringStream.str(), AT);
  }
  FilterFacade::s_registeredFilters1StnByName[name] = creator;
}

void FilterFacade::registerFilter(const string& name, const FilterMStnCreator creator) {
  if (!initialized) registerAllFilters();
  stringstream tmpStringStream;
  if (name.empty()) {
    tmpStringStream<<"Invalid empty name for filter.";
    THROW InvalidArgumentException(tmpStringStream.str(), AT);
  }
  if (!creator) {
    tmpStringStream<<"Invalid null creator for filter '"<<name<<"'.";
    THROW InvalidArgumentException(tmpStringStream.str(), AT);
  }
  if (FilterFacade::s_registeredFilters1StnByName.find(name) != FilterFacade::s_registeredFilters1StnByName.end()) {
    tmpStringStream<<"A single-station filter already exists with name '"<<name<<"'.";
    THROW InvalidArgumentException(tmpStringStream.str(), AT);
  }
  if (FilterFacade::s_registeredFiltersMStnByName.find(name) != FilterFacade::s_registeredFiltersMStnByName.end()) {
    tmpStringStream<<"A filter already exists with name '"<<name<<"'.";
    THROW InvalidArgumentException(tmpStringStream.str(), AT);
  }
  FilterFacade::s_registeredFiltersMStnByName[name] = creator;
}

map<string, Filter1StnCreator>* FilterFacade::getRegisteredFilters1Stn() {
  if (!initialized) registerAllFilters();
  return &s_registeredFilters1StnByName;
}

map<string, FilterMStnCreator>* FilterFacade::getRegisteredFiltersMStn() {
  if (!initialized) registerAllFilters();
  return &s_registeredFiltersMStnByName;
}

vector<FilterBase1Stn*>* FilterFacade::getActiveFilters1Stn() {
  return &m_activeFilters1Stn;
}

vector<FilterBaseMStn*>* FilterFacade::getActiveFiltersMStn() {
  return &m_activeFiltersMStn;
}

void FilterFacade::addActiveFilters(const string& filename) {
  ifstream inputFileStream;
  stringstream tmpStringStream;
  int lineType = ConfigReader::CfgLineComment /*dummy neutral initial value */;
  string str1, str2, section;
  unsigned int lineNb = 0;
  FilterBase1Stn* filter1Stn = NULL;
  FilterBaseMStn* filterMStn = NULL;

  if (!slfutils::fileExists(filename)) THROW FileNotFoundException(filename, AT);
  inputFileStream.open (filename.c_str(), ifstream::in);
  if (inputFileStream.fail()) THROW FileAccessException(filename,AT);

  while (ConfigReader::readConfigLine(inputFileStream, ++lineNb, lineType, str1, str2)) {
    if (lineType == ConfigReader::CfgLineSection) {
      // sections
      if (filter1Stn) { // save the current filter
        addActiveFilter(filter1Stn);
        filter1Stn = NULL;
      } else if (filterMStn) { // save the current filter
        addActiveFilter(filterMStn);
        filterMStn = NULL;
      }
      section = str1;
      // the next filter will be created later
    } else if (lineType == ConfigReader::CfgLineKeyValue) {
      if (str1 == c_name) {
        // "name" parameter
        if (!filter1Stn && !filterMStn) {
          if (FilterFacade::s_registeredFilters1StnByName.find(str2) != FilterFacade::s_registeredFilters1StnByName.end()) {
            Filter1StnCreator creator = FilterFacade::s_registeredFilters1StnByName[str2];
            filter1Stn = (*creator)(); // create the filter
            filter1Stn->setParamValue(c_sectionName, section);
          } else if (FilterFacade::s_registeredFiltersMStnByName.find(str2) != FilterFacade::s_registeredFiltersMStnByName.end()) {
            FilterMStnCreator creator = FilterFacade::s_registeredFiltersMStnByName[str2];
            filterMStn = (*creator)(); // create the filter
            filterMStn->setParamValue(c_sectionName, section);
          } else {
            tmpStringStream<<"Unknown filter named '"<<str2<<"'.";
            THROW UnknownValueException(tmpStringStream.str(), AT);
          }
        } else {
          tmpStringStream<<"Missplaced '"<<c_name<<"' parameter at line "<<lineNb<<" , it must be just after a section start.";
          THROW InvalidFormatException(tmpStringStream.str(), AT);
        }
      } else if (str1 == c_sectionName) {
          tmpStringStream<<"Unallowed used of reserved '"<<c_sectionName<<"' parameter at line "<<lineNb<<".";
          THROW InvalidFormatException(tmpStringStream.str(), AT);
      } else {
        // other parameters
        if (filter1Stn) {
          filter1Stn->setParamValue(str1, str2);
        } else if (filterMStn) {
          filterMStn->setParamValue(str1, str2);
        } else {
          tmpStringStream<<"Missplaced parameter at line "<<lineNb<<" , the first parameter of a filter has to be '"<<c_name<<"'.";
          THROW InvalidFormatException(tmpStringStream.str(), AT);
        }
      }
    }
  }
  if (filter1Stn) { // save the last filter
    addActiveFilter(filter1Stn);
    filter1Stn = NULL;
  } else if (filterMStn) {
    addActiveFilter(filterMStn);
    filterMStn = NULL;
  }
}

void FilterFacade::addActiveFilter(FilterBase1Stn* filter) {
  m_activeFilters1Stn.push_back(filter);
}

void FilterFacade::addActiveFilter(FilterBaseMStn* filter) {
  m_activeFiltersMStn.push_back(filter);
}

bool FilterFacade::removeActiveFilter(const FilterBase1Stn* filter) {
  vector<FilterBase1Stn*>::iterator iter;
  for( iter = m_activeFilters1Stn.begin(); iter != m_activeFilters1Stn.end(); ++iter ) {
    if (*iter == filter) {
      m_activeFilters1Stn.erase(iter);
      return true;
    }
  }
  return false;
}

bool FilterFacade::removeActiveFilter(const FilterBaseMStn* filter) {
  vector<FilterBaseMStn*>::iterator iter;
  for( iter = m_activeFiltersMStn.begin(); iter != m_activeFiltersMStn.end(); ++iter ) {
    if (*iter == filter) {
      m_activeFiltersMStn.erase(iter);
      return true;
    }
  }
  return false;
}

void FilterFacade::removeActiveFilters() {
  m_activeFilters1Stn.clear();
  m_activeFiltersMStn.clear();
}

void FilterFacade::clearActiveFilters() {
  vector<FilterBase1Stn*>::iterator iter1;
  for( iter1 = m_activeFilters1Stn.begin(); iter1 != m_activeFilters1Stn.end(); ++iter1 ) {
    delete (*iter1);
    (*iter1) = NULL;
  }
  m_activeFilters1Stn.clear();
  vector<FilterBaseMStn*>::iterator iterM;
  for( iterM = m_activeFiltersMStn.begin(); iterM != m_activeFiltersMStn.end(); ++iterM ) {
    delete (*iterM);
    (*iterM) = NULL;
  }
  m_activeFiltersMStn.clear();
}

void FilterFacade::getMinimalWindow(unsigned int& minNbPoints, Date& minDeltaTime) {
  minNbPoints = 1;
  minDeltaTime = Date(2000, 1, 1, 0, 0) - Date(2000, 1, 1, 0, 0);
  unsigned int minNbPointsTmp = minNbPoints;
  Date minDeltaTimeTmp = minDeltaTime;
  vector<FilterBase1Stn*>::iterator iter1;
  for( iter1 = m_activeFilters1Stn.begin(); iter1 != m_activeFilters1Stn.end(); ++iter1 ) {
    (*iter1)->getMinimalWindow(minNbPointsTmp, minDeltaTimeTmp);
    if (minNbPointsTmp > minNbPoints) minNbPoints = minNbPointsTmp;
    if (minDeltaTimeTmp > minDeltaTime) minDeltaTime = minDeltaTimeTmp;
  }
  vector<FilterBaseMStn*>::iterator iterM;
  for( iterM = m_activeFiltersMStn.begin(); iterM != m_activeFiltersMStn.end(); ++iterM ) {
    (*iterM)->getMinimalWindow(minNbPointsTmp, minDeltaTimeTmp);
    if (minNbPointsTmp > minNbPoints) minNbPoints = minNbPointsTmp;
    if (minDeltaTimeTmp > minDeltaTime) minDeltaTime = minDeltaTimeTmp;
  }
}

void FilterFacade::prepareCheck() {
  vector<FilterBase1Stn*>::iterator iter1;
  for( iter1 = m_activeFilters1Stn.begin(); iter1 != m_activeFilters1Stn.end(); ++iter1 ) {
    (*iter1)->prepareCheck();
  }
  vector<FilterBaseMStn*>::iterator iterM;
  for( iterM = m_activeFiltersMStn.begin(); iterM != m_activeFiltersMStn.end(); ++iterM ) {
    (*iterM)->prepareCheck();
  }
}

void FilterFacade::doCheck(MeteoBuffer& unfilteredMeteoBuffer, MeteoBuffer& filteredMeteoBuffer, unsigned int iUnfilteredElement) {
  // check single-station filters
  unsigned int iFilteredElement;
  stringstream tmpStringStream;
  MeteoData& unfilteredMeteo = unfilteredMeteoBuffer.getMeteoData(iUnfilteredElement);
  iFilteredElement = filteredMeteoBuffer.seek(unfilteredMeteo.date);
  if (iFilteredElement == MeteoBuffer::npos) {
    filteredMeteoBuffer.put(unfilteredMeteo, unfilteredMeteoBuffer.getStationData(iUnfilteredElement));
    iFilteredElement = filteredMeteoBuffer.seek(unfilteredMeteo.date);
  } else if (filteredMeteoBuffer.getMeteoData(iFilteredElement).date == unfilteredMeteo.date) {
    ;// all right
  } else {
    tmpStringStream << "Existing content in filteredMeteoBuffer prevents to insert at a previous date";
    THROW SLFException(tmpStringStream.str(), AT);
  }
  vector<FilterBase1Stn*>::iterator iter;
  for( iter = m_activeFilters1Stn.begin(); iter != m_activeFilters1Stn.end(); ++iter ) {
    (*iter)->doCheck(unfilteredMeteoBuffer, filteredMeteoBuffer, iFilteredElement);
  }
}

void FilterFacade::doCheck(vector<StationData>& stations, vector<MeteoBuffer>& filteredMeteoBuffers, vector<MeteoData>& meteoDatas) {
  // check single-station filters
  for (unsigned int i=0 ; i < stations.size() ; i++) {
    MeteoBuffer newMeteoBuffer(1);
    newMeteoBuffer.put(meteoDatas[i], stations[i]);
    vector<FilterBase1Stn*>::iterator iter1;
    for( iter1 = m_activeFilters1Stn.begin(); iter1 != m_activeFilters1Stn.end(); ++iter1 ) {
      (*iter1)->doCheck(filteredMeteoBuffers[i], newMeteoBuffer, 0);
    }
    meteoDatas[i] = newMeteoBuffer.getMeteoData(0);
  }
  // check multi-stations filters
  vector<FilterBaseMStn*>::iterator iterM;
  for( iterM = m_activeFiltersMStn.begin(); iterM != m_activeFiltersMStn.end(); ++iterM ) {
    (*iterM)->doCheck(stations, filteredMeteoBuffers, meteoDatas);
  }
}
