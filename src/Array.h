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
		T getMin();
		T getMax();

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

template<class T> T Array<T>::getMin() {

	T min = std::numeric_limits<T>::max();

	for (unsigned int ii=0; ii<arraySize; ii++) {
		const T val = vecData[ii];
		if(val<min) min=val;
	}
	
	return min;
}

template<class T> T Array<T>::getMax() {

	T max = -std::numeric_limits<T>::max();

	for (unsigned int ii=0; ii<arraySize; ii++) {
		const T val = vecData[ii];
		if(val>max) max=val;
	}

	return max;
}

struct SPECIAL_PTS
{
  unsigned int ix;
  unsigned int iy;
};

typedef Array<SPECIAL_PTS> CSpecialPTSArray;

#endif
