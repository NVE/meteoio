#ifndef ARRAY_H
#define ARRAY_H

#include <vector>
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

struct SPECIAL_PTS
{
  unsigned int ix;
  unsigned int iy;
};

typedef Array<SPECIAL_PTS> CSpecialPTSArray;

#endif
