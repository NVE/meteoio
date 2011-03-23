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
#include <iostream>

#include <meteoio/IOUtils.h>
#include <meteoio/IOExceptions.h>

#define NOSAFECHECKS
namespace mio {

/**
 * @class Array
 * @brief The template class Array is a 1D array (vector) able to hold any type of object as datatype
 *
 * @ingroup data_str
 * @author Thomas Egger
 * @date   2009-05-02
 */

template<class T> class Array {
	public:
		Array(const unsigned int& asize=0);

		/**
		* A constructor that creates an array filled with constant values
		* @param asize size of the new array
		* @param init initial value to fill the array with
		*/
		Array(const unsigned int& asize, const T& init);

		unsigned int size();
		void resize(const unsigned int& asize);
		void resize(const unsigned int& asize, const T& init);
		void clear();
		void insertAt(const int& index, T e);
		void removeAt(const unsigned int& index);
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
		/**
		* @brief returns the sum of values contained in the grid
		* @param flag_nodata specify how to process nodata values (see NODATA_HANLDING)
		* @return sum
		*/
		T getSum(const IOUtils::nodata_handling flag_nodata=IOUtils::PARSE_NODATA) const;
		/**
		* @brief returns the grid of the absolute value of values contained in the grid
		* @param flag_nodata specify how to process nodata values (see NODATA_HANLDING)
		* @return grid of abs(grid)
		*/
		const Array<T> getAbs(const IOUtils::nodata_handling flag_nodata=IOUtils::PARSE_NODATA) const;
		void abs(const IOUtils::nodata_handling flag_nodata=IOUtils::PARSE_NODATA);


		template<class P> friend std::ostream& operator<<(std::ostream& os, const Array<P>& array);
		T& operator [](const unsigned int& index);
		const T operator [](const unsigned int& index) const;
		T& operator ()(const unsigned int& index);
		const T operator ()(const unsigned int& index) const;

		Array<T>& operator =(const Array<T>&);
		Array<T>& operator =(const T& value);
		
		Array<T>& operator+=(const T& rhs);
		const Array<T> operator+(const T& rhs);
		Array<T>& operator+=(const Array<T>& rhs);
		const Array<T> operator+(const Array<T>& rhs);

		Array<T>& operator-=(const T& rhs);
		const Array<T> operator-(const T& rhs);
		Array<T>& operator-=(const Array<T>& rhs);
		const Array<T> operator-(const Array<T>& rhs);

		Array<T>& operator*=(const T& rhs);
		const Array<T> operator*(const T& rhs);
		Array<T>& operator*=(const Array<T>& rhs);
		const Array<T> operator*(const Array<T>& rhs);

		Array<T>& operator/=(const T& rhs);
		const Array<T> operator/(const T& rhs);
		Array<T>& operator/=(const Array<T>& rhs);
		const Array<T> operator/(const Array<T>& rhs);

	protected:
		std::vector<T> vecData; ///<the actual data structure, that holds the objects of type T
		unsigned int nx; ///<this is introduced to omit the costly vecData.size()
};

template<class T> Array<T>::Array(const unsigned int& asize) {
	resize(asize);
}

template<class T> Array<T>::Array(const unsigned int& asize, const T& init) {
	resize(asize, init);
}

template<class T> unsigned int Array<T>::size() {
	return nx;
}

template<class T> void Array<T>::resize(const unsigned int& asize) {
	vecData.resize(asize);
	nx = asize;
}

template<class T> void Array<T>::resize(const unsigned int& asize, const T& init) {
	resize(asize);
	std::fill(vecData.begin(), vecData.end(), init);
}


template<class T> T& Array<T>::operator()(const unsigned int& index) {
#ifndef NOSAFECHECKS
	if (index >= nx) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif
	return vecData[index];
}

template<class T> const T Array<T>::operator()(const unsigned int& index) const {
#ifndef NOSAFECHECKS
	if (index >= nx) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif
	return vecData[index];
}

template<class T> T& Array<T>::operator [](const unsigned int& index) {
#ifndef NOSAFECHECKS
	if (index >= nx) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif
	
	return vecData[index];
}

template<class T> const T Array<T>::operator [](const unsigned int& index) const {
#ifndef NOSAFECHECKS
	if (index >= nx) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif

	return vecData[index];
}

template<class T> void Array<T>::clear() {
	vecData.clear();
	nx = 0;
}

template<class T> std::ostream& operator<<(std::ostream& os, const Array<T>& array) {
	os << "<array1d>\n";
	for(unsigned int ii=0; ii<array.nx; ii++) {
		os << array(ii) << " ";
	}
	os << "\n</array1d>\n";
	return os;
}

template<class T> void Array<T>::insertAt(const int& index, T e) {
	if (index < 0) {
		vecData.push_back(e);
                nx++;
	} else if ((index >= 0) && (index < (int)vecData.size())) {
		vecData.insert(vecData.begin() + index, e);
		nx++;
	} else {
		throw IndexOutOfBoundsException("", AT);
	}
}

template<class T> void Array<T>::removeAt(const unsigned int& index) {
	if (index < vecData.size()) {
		vecData.erase(vecData.begin()+index);
		nx--;
	}
}

template<class T> T Array<T>::getMin(const IOUtils::nodata_handling flag_nodata) const {

	T min = std::numeric_limits<T>::max();

	if(flag_nodata==IOUtils::RAW_NODATA) {
		for (unsigned int ii=0; ii<nx; ii++) {
			const T val = vecData[ii];
			if(val<min) min=val;
		}
		return min;
	} else if(flag_nodata==IOUtils::PARSE_NODATA) {
		for (unsigned int ii=0; ii<nx; ii++) {
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
		for (unsigned int ii=0; ii<nx; ii++) {
			const T val = vecData[ii];
			if(val>max) max=val;
		}
		return max;
	} else if(flag_nodata==IOUtils::PARSE_NODATA) {
		for (unsigned int ii=0; ii<nx; ii++) {
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
		for (unsigned int ii=0; ii<nx; ii++) {
			const T val = vecData[ii];
			mean += val;
		}
		const unsigned int count = nx;
		if(count>0) return mean/(T)(count);
		else return (T)0;
	} else if(flag_nodata==IOUtils::PARSE_NODATA) {
		unsigned int count = 0;
		for (unsigned int ii=0; ii<nx; ii++) {
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

template<class T> T Array<T>::getSum(const IOUtils::nodata_handling flag_nodata) const {

	T sum = 0;

	if(flag_nodata==IOUtils::RAW_NODATA) {
		for (unsigned int ii=0; ii<nx; ii++) {
			sum += vecData[ii];
		}
		return sum;
	} else if(flag_nodata==IOUtils::PARSE_NODATA) {
		bool empty=true;
		for (unsigned int ii=0; ii<nx; ii++) {
			const T val = vecData[ii];
			if(val!=IOUtils::nodata) {
				sum += val;
				empty=false;
			}
		}
		if(!empty) return sum;
		else return (T)IOUtils::nodata;
	} else {
		throw InvalidArgumentException("Unknown nodata_handling flag",AT);
	}
}

template<class T> void Array<T>::abs(const IOUtils::nodata_handling flag_nodata) {
	if(std::numeric_limits<T>::is_signed) {
		if(flag_nodata==IOUtils::RAW_NODATA) {
			for (unsigned int ii=0; ii<nx; ii++) {
				T& val = operator()(ii);
				if(val<0) val=-val;
			}
		} else if(flag_nodata==IOUtils::PARSE_NODATA) {
			for (unsigned int ii=0; ii<nx; ii++) {
				T& val = operator()(ii);
				if(val<0 && val!=IOUtils::nodata) val=-val;
			}
		} else {
			throw InvalidArgumentException("Unknown nodata_handling flag",AT);
		}
	}
}


template<class T> const Array<T> Array<T>::getAbs(const IOUtils::nodata_handling flag_nodata) const {
	Array<T> result = *this; //make a copy
	result.abs(flag_nodata); //already implemented

	return result;
}

//arithmetic operators
template<class T> Array<T>& Array<T>::operator=(const Array<T>& source) {
	if(this != &source) {
		vecData = source.vecData;
		nx = source.nx;
	}
	return *this;
}

template<class T> Array<T>& Array<T>::operator=(const T& value) {
	//reset every single member of the Array3D<T>
	std::fill(vecData.begin(), vecData.end(), value);
	return *this;
}

template<class T> Array<T>& Array<T>::operator+=(const Array<T>& rhs)
{
	//They have to have equal size
	if (rhs.nx != nx)
		throw IOException("Trying to add two Array objects with different dimensions", AT);

	//Add to every single member of the Array<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		operator()(ii) += rhs(ii);
	}

	return *this;
}

template<class T> const Array<T> Array<T>::operator+(const Array<T>& rhs)
{
	Array<T> result = *this; //make a copy
	result += rhs; //already implemented

	return result;
}

template<class T> Array<T>& Array<T>::operator+=(const T& rhs)
{
	//Add to every single member of the Array<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		operator()(ii) += rhs;
	}

	return *this;
}

template<class T> const Array<T> Array<T>::operator+(const T& rhs)
{
	Array<T> result = *this;
	result += rhs; //already implemented

	return result;
}

template<class T> Array<T>& Array<T>::operator-=(const Array<T>& rhs)
{
	//They have to have equal size
	if (rhs.nx != nx)
		throw IOException("Trying to substract two Array objects with different dimensions", AT);

	//Substract to every single member of the Array<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		operator()(ii) -= rhs(ii);
	}

	return *this;
}

template<class T> const Array<T> Array<T>::operator-(const Array<T>& rhs)
{
	Array<T> result = *this; //make a copy
	result -= rhs; //already implemented

	return result;
}

template<class T> Array<T>& Array<T>::operator-=(const T& rhs)
{
	//Substract to every single member of the Array<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		operator()(ii) -= rhs;
	}

	return *this;
}

template<class T> const Array<T> Array<T>::operator-(const T& rhs)
{
	Array<T> result = *this;
	result -= rhs; //already implemented

	return result;
}

template<class T> Array<T>& Array<T>::operator*=(const Array<T>& rhs)
{
	//They have to have equal size
	if (rhs.nx != nx)
		throw IOException("Trying to multiply two Array objects with different dimensions", AT);

	//Add to every single member of the Array<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		operator()(ii) *= rhs(ii);
	}

	return *this;
}

template<class T> const Array<T> Array<T>::operator*(const Array<T>& rhs)
{
	Array<T> result = *this; //make a copy
	result *= rhs; //already implemented

	return result;
}

template<class T> Array<T>& Array<T>::operator*=(const T& rhs)
{
	//Add to every single member of the Array<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		operator()(ii) *= rhs;
	}

	return *this;
}

template<class T> const Array<T> Array<T>::operator*(const T& rhs)
{
	Array<T> result = *this;
	result *= rhs; //already implemented

	return result;
}

template<class T> Array<T>& Array<T>::operator/=(const Array<T>& rhs)
{
	//They have to have equal size
	if (rhs.nx != nx)
		throw IOException("Trying to divide two Array objects with different dimensions", AT);

	//Divide every single member of the Array<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		operator()(ii) /= rhs(ii);
	}

	return *this;
}

template<class T> const Array<T> Array<T>::operator/(const Array<T>& rhs)
{
	Array<T> result = *this; //make a copy
	result /= rhs; //already implemented

	return result;
}

template<class T> Array<T>& Array<T>::operator/=(const T& rhs)
{
	//Divide every single member of the Array<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		operator()(ii) /= rhs;
	}

	return *this;
}

template<class T> const Array<T> Array<T>::operator/(const T& rhs)
{
	Array<T> result = *this;
	result /= rhs; //already implemented

	return result;
}

} //end namespace mio

#endif
