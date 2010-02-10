/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ARRAY_H
#define ARRAY_H

#include "IOUtils.h"
#include <vector>
#include <limits>
#include "IOExceptions.h"

#define NOSAFECHECKS

/**
 * @class Array
 * @brief The template class Array is a 1D array (vector) able to hold any type of object as datatype
 *
 * @author Thomas Egger
 * @date   2009-05-02
 */
template<class T> class Array {
	public:
		Array(unsigned int asize=0);

		T& operator [](unsigned int index);
		const T operator [](unsigned int index) const;

		unsigned int size();
		void resize(unsigned int asize);
		void clear();
		void insertAt(int index, T e);
		void removeAt(unsigned int index);
		/**
		* @brief returns the minimum value contained in the grid
		* @param flag_nodata specify how to process nodata values (see NODATA_HANLDING)
		* @return minimum value
		*/
		T getMin(const IOUtils::nodata_handling flag_nodata=IOUtils::PARSE_NODATA) const;
		/**
		* @brief returns the maximum value contained in the grid
		* @param flag_nodata specify how to process nodata values (see NODATA_HANLDING)
		* @return maximum value
		*/
		T getMax(const IOUtils::nodata_handling flag_nodata=IOUtils::PARSE_NODATA) const;
		/**
		* @brief returns the mean value contained in the grid
		* @param flag_nodata specify how to process nodata values (see NODATA_HANLDING)
		* @return mean value
		*/
		T getMean(const IOUtils::nodata_handling flag_nodata=IOUtils::PARSE_NODATA) const;


	protected:
		std::vector<T> vecData; ///<the actual data structure, that holds the objects of type T
		unsigned int arraySize; ///<this is introduced to omit the costly vecData.size()
};

template<class T> Array<T>::Array(unsigned int asize) {
	resize(asize);
	arraySize = asize;
}

template<class T> unsigned int Array<T>::size() {
	return arraySize;
}

template<class T> void Array<T>::resize(unsigned int asize) {
	if (asize != vecData.size()) {
		vecData.resize(asize);
		arraySize = asize;
	}
}

template<class T> T& Array<T>::operator [](unsigned int index) {
#ifndef NOSAFECHECKS
	if (index >= arraySize) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif
	
	return vecData[index];
}

template<class T> const T Array<T>::operator [](unsigned int index) const {
#ifndef NOSAFECHECKS
	if (index >= arraySize) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif

	return vecData[index];
}

template<class T> void Array<T>::clear() {
	vecData.clear();
	arraySize = 0;
}

template<class T> void Array<T>::insertAt(int index, T e) {
	if (index < 0) {
		vecData.push_back(e);
                arraySize++;
	} else if ((index >= 0) && (index < (int)vecData.size())) {
		vecData.insert(vecData.begin() + index, e);
		arraySize++;
	} else {
		throw IndexOutOfBoundsException("", AT);
	}
}

template<class T> void Array<T>::removeAt(unsigned int index) {
	if (index < vecData.size()) {
		vecData.erase(vecData.begin()+index);
		arraySize--;
	}
}

template<class T> T Array<T>::getMin(const IOUtils::nodata_handling flag_nodata) const {

	T min = std::numeric_limits<T>::max();

	if(flag_nodata==IOUtils::RAW_NODATA) {
		for (unsigned int ii=0; ii<arraySize; ii++) {
			const T val = vecData[ii];
			if(val<min) min=val;
		}
		return min;
	} else if(flag_nodata==IOUtils::PARSE_NODATA) {
		for (unsigned int ii=0; ii<arraySize; ii++) {
			const T val = vecData[ii];
			if(val!=IOUtils::nodata && val<min) min=val;
		}
		if(min!=std::numeric_limits<T>::max()) return min;
		else return (T)IOUtils::nodata;
	} else {
		throw InvalidArgumentException("Unknown nodata_handling flag",AT);
	}
}

template<class T> T Array<T>::getMax(const IOUtils::nodata_handling flag_nodata) const {

	T max = -std::numeric_limits<T>::max();

	if(flag_nodata==IOUtils::RAW_NODATA) {
		for (unsigned int ii=0; ii<arraySize; ii++) {
			const T val = vecData[ii];
			if(val>max) max=val;
		}
		return max;
	} else if(flag_nodata==IOUtils::PARSE_NODATA) {
		for (unsigned int ii=0; ii<arraySize; ii++) {
			const T val = vecData[ii];
			if(val!=IOUtils::nodata && val>max) max=val;
		}
		if(max!=-std::numeric_limits<T>::max()) return max;
		else return (T)IOUtils::nodata;
	} else {
		throw InvalidArgumentException("Unknown nodata_handling flag",AT);
	}
}

template<class T> T Array<T>::getMean(const IOUtils::nodata_handling flag_nodata) const {

	T mean = 0;

	if(flag_nodata==IOUtils::RAW_NODATA) {
		for (unsigned int ii=0; ii<arraySize; ii++) {
			const T val = vecData[ii];
			mean += val;
		}
		const unsigned int count = arraySize;
		if(count>0) return mean/(T)(count);
		else return (T)0;
	} else if(flag_nodata==IOUtils::PARSE_NODATA) {
		unsigned int count = 0;
		for (unsigned int ii=0; ii<arraySize; ii++) {
			const T val = vecData[ii];
			if(val!=IOUtils::nodata) {
				mean += val;
				count++;
			}
		}
		if(count>0) return mean/(T)(count);
		else return (T)IOUtils::nodata;
	} else {
		throw InvalidArgumentException("Unknown nodata_handling flag",AT);
	}
}

struct POINT
{
  unsigned int ix;
  unsigned int iy;
};

typedef Array<POINT> POINTSArray;

#endif
