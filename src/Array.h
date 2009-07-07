#ifndef ARRAY_H
#define ARRAY_H

#include <vector>
#include "IOExceptions.h"

/**
 * @class CArray
 * @brief The template class CArray is a 1D array (vector) able to hold any type of object as datatype
 *
 * @author Thomas Egger
 * @date   2009-05-02
 */
template<class T> class CArray {
	public:
		CArray(unsigned int asize=0);
		~CArray();

		T& operator [](unsigned int index);
		const T operator [](unsigned int index) const;
		CArray<T>& operator =(CArray<T>& val);

		unsigned int GetSize();
		void SetSize(unsigned int asize);
		void RemoveAll();
		void InsertAt(int index, T e);
		void RemoveAt(unsigned int index);

	protected:
		std::vector<T> vecData; ///<the actual data structure, that holds the objects of type T
		unsigned int size;      ///<this is introduced to omit the costly vecData.size()
};

template<class T> CArray<T>::CArray(unsigned int asize) {
	SetSize(asize);
	size = asize;
}

template<class T> CArray<T>::~CArray() {
	RemoveAll();
}

template<class T> unsigned int CArray<T>::GetSize() {
	return size;
}

template<class T> void CArray<T>::SetSize(unsigned int asize) {
	if (asize != vecData.size()) {
		vecData.resize(asize);
		size = asize;
	}
}

template<class T> T& CArray<T>::operator [](unsigned int index) {
#ifndef NOSAFECHECKS
	if (index >= size) {
		THROW IndexOutOfBoundsException("", AT);
	}
#endif
	
	return vecData[index];
}

template<class T> const T CArray<T>::operator [](unsigned int index) const {
#ifndef NOSAFECHECKS
	if (index >= size) {
		THROW IndexOutOfBoundsException("", AT);
	}
#endif

	return vecData[index];
}

template<class T> CArray<T>& CArray<T>::operator=(CArray & val) {
	vecData = val.vecData;
	size = val.size;
	return *this;
}

template<class T> void CArray<T>::RemoveAll() {
	vecData.clear();
	size = 0;
}

template<class T> void CArray<T>::InsertAt(int index, T e) {
	if (index < 0) {
		vecData.push_back(e);
                size++;
	} else if ((index >= 0) && (index < (int)vecData.size())) {
		vecData.insert(vecData.begin() + index, e);
		size++;
	} else {
		THROW IndexOutOfBoundsException("", AT);
	}
}

template<class T> void CArray<T>::RemoveAt(unsigned int index) {
	if (index < vecData.size()) {
		vecData.erase(vecData.begin()+index);
		size--;
	}
}

struct SPECIAL_PTS
{
  int ix;
  int iy;
};

typedef CArray<SPECIAL_PTS> CSpecialPTSArray;

#endif
