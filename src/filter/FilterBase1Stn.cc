#include "FilterBase1Stn.h"

using namespace std;

MeteoBufferIterator::MeteoBufferIterator(MeteoBuffer& unfilteredMeteoBuffer, MeteoBuffer& filteredMeteoBuffer, unsigned int iFilteredElement) :
  m_unfilteredMeteoBuffer(unfilteredMeteoBuffer),
  m_filteredMeteoBuffer(filteredMeteoBuffer),
  m_iUnfilteredElement(MeteoBuffer::npos),
  m_iFilteredElement(iFilteredElement)
{
}

MeteoData& MeteoBufferIterator::getCurrent()
{
	if (m_iFilteredElement != MeteoBuffer::npos) {
		return m_filteredMeteoBuffer.getMeteoData(m_iFilteredElement);
	} else if (m_iUnfilteredElement != MeteoBuffer::npos) {
		return m_unfilteredMeteoBuffer.getMeteoData(m_iUnfilteredElement);
	} else {
		THROW NoAvailableDataException("no corresponding data available in the buffers", AT);
	}
}

MeteoData& MeteoBufferIterator::getPrevious()
{
	if (m_iFilteredElement == MeteoBuffer::npos) {
		return getPreviousUnfiltered();
	} else if (m_iFilteredElement == 0) {
		// switch from filtered buffer to unfiltered buffer, which has usually a bigger history
		return getPreviousUnfiltered();
	} else {
		m_iFilteredElement--;
		return m_filteredMeteoBuffer.getMeteoData(m_iFilteredElement);
	}
}

MeteoData& MeteoBufferIterator::getPreviousUnfiltered()
{
	if (m_iFilteredElement != MeteoBuffer::npos) {
		// switch from filtered buffer to unfiltered buffer
		Date_IO& currDate = m_filteredMeteoBuffer.getMeteoData(m_iFilteredElement).date;
		m_iFilteredElement = MeteoBuffer::npos;
		m_iUnfilteredElement = m_unfilteredMeteoBuffer.seek(currDate);
	}
	
	if (m_iUnfilteredElement == MeteoBuffer::npos || m_iUnfilteredElement == 0) {
		THROW NoAvailableDataException("no previous data available in the buffers", AT);
	} else {
		m_iUnfilteredElement--;
		return m_unfilteredMeteoBuffer.getMeteoData(m_iUnfilteredElement);
	}
}

