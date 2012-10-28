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
#ifndef ARRAY1D_H
#define ARRAY1D_H

#include <vector>
#include <limits>
#include <iostream>

#include <meteoio/IOUtils.h>
#include <meteoio/IOExceptions.h>

namespace mio {

/**
 * @class Array1D
 * @brief The template class Array1D is a 1D array (vector) able to hold any type of object as datatype.
 * If the compilation flag NOSAFECHECKS is used, bounds check is turned off (leading to increased performances).
 *
 * @ingroup data_str
 * @author Thomas Egger
 * @date   2009-05-02
 */

template<class T> class Array1D {
	public:
		Array1D(const unsigned int& asize=0);

		/**
		* A constructor that creates an array filled with constant values
		* @param asize size of the new array
		* @param init initial value to fill the array with
		*/
		Array1D(const unsigned int& asize, const T& init);

		/**
		* @brief set how to process nodata values (ie: as nodata or as normal numbers)
		* @param i_keep_nodata true means that NODATA is interpreted as NODATA, false means that it is a normal number
		By default, arrays keep nodata.
		*/
		void setKeepNodata(const bool i_keep_nodata);

		/**
		* @brief get how to process nodata values (ie: as nodata or as normal numbers)
		* @return true means that NODATA is interpreted as NODATA, false means that it is a normal number
		*/
		bool getKeepNodata();

		void size(unsigned int& nx) const;
		unsigned int getNx() const;

		void resize(const unsigned int& asize);
		void resize(const unsigned int& asize, const T& init);
		void clear();
		bool isEmpty() const;
		void insertAt(const int& index, T e);
		void removeAt(const unsigned int& index);
		/**
		* @brief returns the minimum value contained in the grid
		* @return minimum value
		*/
		T getMin() const;
		/**
		* @brief returns the maximum value contained in the grid
		* @return maximum value
		*/
		T getMax() const;
		/**
		* @brief returns the mean value contained in the grid
		* @return mean value
		*/
		T getMean() const;
		/**
		* @brief returns the number of points contained in the grid.
		* If setNodataHandling(IOUtils::RAW_NODATA), then the number of points is the size of the grid.
		* If setNodataHandling(IOUtils::PARSE_NODATA), then it is the number of non-nodata values in the grid
		* @return count
		*/
		size_t getCount() const;
		/**
		* @brief returns the grid of the absolute value of values contained in the grid
		* @return grid of abs(grid)
		*/
		const Array1D<T> getAbs() const;
		void abs();


		template<class P> friend std::ostream& operator<<(std::ostream& os, const Array1D<P>& array);

		bool checkEpsilonEquality(const Array1D<double>& rhs, const double& epsilon) const;
		static bool checkEpsilonEquality(const Array1D<double>& rhs1, const Array1D<double>& rhs2, const double& epsilon);

		T& operator [](const unsigned int& index);
		const T operator [](const unsigned int& index) const;
		T& operator ()(const unsigned int& index);
		const T operator ()(const unsigned int& index) const;

		Array1D<T>& operator =(const Array1D<T>&);
		Array1D<T>& operator =(const T& value);

		Array1D<T>& operator+=(const T& rhs);
		const Array1D<T> operator+(const T& rhs);
		Array1D<T>& operator+=(const Array1D<T>& rhs);
		const Array1D<T> operator+(const Array1D<T>& rhs);

		Array1D<T>& operator-=(const T& rhs);
		const Array1D<T> operator-(const T& rhs);
		Array1D<T>& operator-=(const Array1D<T>& rhs);
		const Array1D<T> operator-(const Array1D<T>& rhs);

		Array1D<T>& operator*=(const T& rhs);
		const Array1D<T> operator*(const T& rhs);
		Array1D<T>& operator*=(const Array1D<T>& rhs);
		const Array1D<T> operator*(const Array1D<T>& rhs);

		Array1D<T>& operator/=(const T& rhs);
		const Array1D<T> operator/(const T& rhs);
		Array1D<T>& operator/=(const Array1D<T>& rhs);
		const Array1D<T> operator/(const Array1D<T>& rhs);

	protected:
		std::vector<T> vecData; ///<the actual data structure, that holds the objects of type T
		unsigned int nx; ///<this is introduced to omit the costly vecData.size()
		bool keep_nodata;
};

template<class T> Array1D<T>::Array1D(const unsigned int& asize)
                  : vecData(asize), nx(asize), keep_nodata(true)
{
	//resize(asize);
}

template<class T> Array1D<T>::Array1D(const unsigned int& asize, const T& init)
                 : vecData(asize, init), nx(asize), keep_nodata(true)
{
	//resize(asize, init);
}

template<class T> void Array1D<T>::setKeepNodata(const bool i_keep_nodata) {
	keep_nodata = i_keep_nodata;
}

template<class T> bool Array1D<T>::getKeepNodata() {
	return keep_nodata;
}

template<class T> void Array1D<T>::size(unsigned int& o_nx) const {
	o_nx = nx;
}

template<class T> unsigned int Array1D<T>::getNx() const {
	return nx;
}

template<class T> void Array1D<T>::resize(const unsigned int& asize) {
	vecData.clear();
	vecData.resize(asize);
	nx = asize;
}

template<class T> void Array1D<T>::resize(const unsigned int& asize, const T& init) {
	vecData.clear();
	vecData.resize(asize, init);
	nx = asize;
}


template<class T> inline T& Array1D<T>::operator()(const unsigned int& index) {
#ifndef NOSAFECHECKS
	if (index >= nx) {
		std::stringstream ss;
		ss << "Trying to access array(" << index << ")";
		ss << " while array is (" << nx << ")";
		throw IndexOutOfBoundsException(ss.str(), AT);
	}
#endif
	return vecData[index];
}

template<class T> inline const T Array1D<T>::operator()(const unsigned int& index) const {
#ifndef NOSAFECHECKS
	if (index >= nx) {
		std::stringstream ss;
		ss << "Trying to access array(" << index << ")";
		ss << " while array is (" << nx << ")";
		throw IndexOutOfBoundsException(ss.str(), AT);
	}
#endif
	return vecData[index];
}

template<class T> inline T& Array1D<T>::operator [](const unsigned int& index) {
#ifndef NOSAFECHECKS
	return vecData.at(index);
#else
	return vecData[index];
#endif
}

template<class T> inline const T Array1D<T>::operator [](const unsigned int& index) const {
#ifndef NOSAFECHECKS
	return vecData.at(index);
#else
	return vecData[index];
#endif
}

template<class T> void Array1D<T>::clear() {
	vecData.clear();
	nx = 0;
}

template<class T> bool Array1D<T>::isEmpty() const {
	return (nx==0);
}

template<class T> std::ostream& operator<<(std::ostream& os, const Array1D<T>& array) {
	os << "<array1d>\n";
	for(unsigned int ii=0; ii<array.nx; ii++) {
		os << array(ii) << " ";
	}
	os << "\n</array1d>\n";
	return os;
}

template<class T> void Array1D<T>::insertAt(const int& index, T e) {
	if (index < 0) {
		vecData.push_back(e);
                nx++;
	} else if ((index >= 0) && (index < (int)vecData.size())) {
		vecData.insert(vecData.begin() + index, e);
		nx++;
	} else {
		std::stringstream ss;
		ss << "Inserting an element at (" << index << ") in an array of size [" << nx << "]";
		throw IndexOutOfBoundsException(ss.str(), AT);
	}
}

template<class T> void Array1D<T>::removeAt(const unsigned int& index) {
	if (index < vecData.size()) {
		vecData.erase(vecData.begin()+index);
		nx--;
	}
}

template<class T> T Array1D<T>::getMin() const {

	T min = std::numeric_limits<T>::max();

	if(keep_nodata==false) {
		for (unsigned int ii=0; ii<nx; ii++) {
			const T val = vecData[ii];
			if(val<min) min=val;
		}
		return min;
	} else {
		for (unsigned int ii=0; ii<nx; ii++) {
			const T val = vecData[ii];
			if(val!=IOUtils::nodata && val<min) min=val;
		}
		if(min!=std::numeric_limits<T>::max()) return min;
		else return (T)IOUtils::nodata;
	}
}

template<class T> T Array1D<T>::getMax() const {

	T max = -std::numeric_limits<T>::max();

	if(keep_nodata==false) {
		for (unsigned int ii=0; ii<nx; ii++) {
			const T val = vecData[ii];
			if(val>max) max=val;
		}
		return max;
	} else {
		for (unsigned int ii=0; ii<nx; ii++) {
			const T val = vecData[ii];
			if(val!=IOUtils::nodata && val>max) max=val;
		}
		if(max!=-std::numeric_limits<T>::max()) return max;
		else return (T)IOUtils::nodata;
	}
}

template<class T> T Array1D<T>::getMean() const {

	T mean = 0;

	if(keep_nodata==false) {
		for (unsigned int ii=0; ii<nx; ii++) {
			const T val = vecData[ii];
			mean += val;
		}
		const unsigned int count = nx;
		if(count>0) return mean/(T)(count);
		else return (T)0;
	} else {
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
	}
}

template<class T> size_t Array1D<T>::getCount() const
{
	if(keep_nodata==false) {
		return (size_t)nx;
	} else {
		size_t count = 0;
		for (unsigned int ii=0; ii<nx; ii++) {
			if(vecData[ii]!=IOUtils::nodata) count++;
		}
		return count;
	}
}

template<class T> void Array1D<T>::abs() {
	if(std::numeric_limits<T>::is_signed) {
		if(keep_nodata==false) {
			for (unsigned int ii=0; ii<nx; ii++) {
				T& val = operator()(ii);
				if(val<0) val=-val;
			}
		} else {
			for (unsigned int ii=0; ii<nx; ii++) {
				T& val = operator()(ii);
				if(val<0 && val!=IOUtils::nodata) val=-val;
			}
		}
	}
}


template<class T> const Array1D<T> Array1D<T>::getAbs() const {
	Array1D<T> result = *this; //make a copy
	result.abs(); //already implemented

	return result;
}

//arithmetic operators
template<class T> bool Array1D<T>::checkEpsilonEquality(const Array1D<double>& rhs, const double& epsilon) const {
	if(nx!=rhs.nx) return false;

	for (unsigned int jj=0; jj<nx; jj++)
		if(IOUtils::checkEpsilonEquality(vecData[jj], rhs.vecData[jj], epsilon)==false) return false;

	return true;
}

template<class T> bool Array1D<T>::checkEpsilonEquality(const Array1D<double>& rhs1, const Array1D<double>& rhs2, const double& epsilon) {
	return rhs1.checkEpsilonEquality(rhs2, epsilon);
}

template<class T> Array1D<T>& Array1D<T>::operator=(const Array1D<T>& source) {
	if(this != &source) {
		vecData = source.vecData;
		nx = source.nx;
		keep_nodata = source.keep_nodata;
	}
	return *this;
}

template<class T> Array1D<T>& Array1D<T>::operator=(const T& value) {
	//reset every single member of the Array1D<T>
	std::fill(vecData.begin(), vecData.end(), value);
	return *this;
}

template<class T> Array1D<T>& Array1D<T>::operator+=(const Array1D<T>& rhs)
{
	//They have to have equal size
	if (rhs.nx != nx) {
		std::stringstream ss;
		ss << "Trying to add two Array1D objects with different dimensions: ";
		ss << "(" << nx << ") + (" << rhs.nx << ")";
		throw IOException(ss.str(), AT);
	}

	//Add to every single member of the Array1D<T>
	if(keep_nodata==false) {
		for (unsigned int ii=0; ii<nx; ii++) {
			operator()(ii) += rhs(ii);
		}
	} else {
		for (unsigned int ii=0; ii<nx; ii++) {
			if(operator()(ii)==IOUtils::nodata || rhs(ii)==IOUtils::nodata)
				operator()(ii) = IOUtils::nodata;
			else
				operator()(ii) += rhs(ii);
		}
	}

	return *this;
}

template<class T> const Array1D<T> Array1D<T>::operator+(const Array1D<T>& rhs)
{
	Array1D<T> result = *this; //make a copy
	result += rhs; //already implemented

	return result;
}

template<class T> Array1D<T>& Array1D<T>::operator+=(const T& rhs)
{
	//Add to every single member of the Array1D<T>
	if(keep_nodata==false) {
		for (unsigned int ii=0; ii<nx; ii++) {
			operator()(ii) += rhs;
		}
	} else {
		for (unsigned int ii=0; ii<nx; ii++) {
			if(operator()(ii)!=IOUtils::nodata)
				operator()(ii) += rhs;
		}
	}

	return *this;
}

template<class T> const Array1D<T> Array1D<T>::operator+(const T& rhs)
{
	Array1D<T> result = *this;
	result += rhs; //already implemented

	return result;
}

template<class T> Array1D<T>& Array1D<T>::operator-=(const Array1D<T>& rhs)
{
	//They have to have equal size
	if (rhs.nx != nx) {
		std::stringstream ss;
		ss << "Trying to substract two Array1D objects with different dimensions: ";
		ss << "(" << nx << ") - (" << rhs.nx << ")";
		throw IOException(ss.str(), AT);
	}

	//Substract to every single member of the Array1D<T>
	if(keep_nodata==false) {
		for (unsigned int ii=0; ii<nx; ii++) {
			operator()(ii) -= rhs(ii);
		}
	} else {
		for (unsigned int ii=0; ii<nx; ii++) {
			if(operator()(ii)==IOUtils::nodata || rhs(ii)==IOUtils::nodata)
				operator()(ii) = IOUtils::nodata;
			else
				operator()(ii) -= rhs(ii);
		}
	}

	return *this;
}

template<class T> const Array1D<T> Array1D<T>::operator-(const Array1D<T>& rhs)
{
	Array1D<T> result = *this; //make a copy
	result -= rhs; //already implemented

	return result;
}

template<class T> Array1D<T>& Array1D<T>::operator-=(const T& rhs)
{
	//Substract to every single member of the Array1D<T>
	if(keep_nodata==false) {
		for (unsigned int ii=0; ii<nx; ii++) {
			operator()(ii) -= rhs;
		}
	} else {
		for (unsigned int ii=0; ii<nx; ii++) {
			if(operator()(ii)!=IOUtils::nodata)
				operator()(ii) -= rhs;
		}
	}

	return *this;
}

template<class T> const Array1D<T> Array1D<T>::operator-(const T& rhs)
{
	Array1D<T> result = *this;
	result -= rhs; //already implemented

	return result;
}

template<class T> Array1D<T>& Array1D<T>::operator*=(const Array1D<T>& rhs)
{
	//They have to have equal size
	if (rhs.nx != nx){
		std::stringstream ss;
		ss << "Trying to multiply two Array1D objects with different dimensions: ";
		ss << "(" << nx << ") * (" << rhs.nx << ")";
		throw IOException(ss.str(), AT);
	}
	//Multiply every single member of the Array1D<T>
	if(keep_nodata==false) {
		for (unsigned int ii=0; ii<nx; ii++) {
			operator()(ii) *= rhs(ii);
		}
	} else {
		for (unsigned int ii=0; ii<nx; ii++) {
			if(operator()(ii)==IOUtils::nodata || rhs(ii)==IOUtils::nodata)
				operator()(ii) = IOUtils::nodata;
			else
				operator()(ii) *= rhs(ii);
		}
	}

	return *this;
}

template<class T> const Array1D<T> Array1D<T>::operator*(const Array1D<T>& rhs)
{
	Array1D<T> result = *this; //make a copy
	result *= rhs; //already implemented

	return result;
}

template<class T> Array1D<T>& Array1D<T>::operator*=(const T& rhs)
{
	//Multiply every single member of the Array1D<T>
	if(keep_nodata==false) {
		for (unsigned int ii=0; ii<nx; ii++) {
			operator()(ii) *= rhs;
		}
	} else {
		for (unsigned int ii=0; ii<nx; ii++) {
			if(operator()(ii)!=IOUtils::nodata)
				operator()(ii) *= rhs;
		}
	}

	return *this;
}

template<class T> const Array1D<T> Array1D<T>::operator*(const T& rhs)
{
	Array1D<T> result = *this;
	result *= rhs; //already implemented

	return result;
}

template<class T> Array1D<T>& Array1D<T>::operator/=(const Array1D<T>& rhs)
{
	//They have to have equal size
	if (rhs.nx != nx){
		std::stringstream ss;
		ss << "Trying to divide two Array1D objects with different dimensions: ";
		ss << "(" << nx << ") / (" << rhs.nx << ")";
		throw IOException(ss.str(), AT);
	}
	//Divide every single member of the Array1D<T>
	if(keep_nodata==false) {
		for (unsigned int ii=0; ii<nx; ii++) {
			operator()(ii) /= rhs(ii);
		}
	} else {
		for (unsigned int ii=0; ii<nx; ii++) {
			if(operator()(ii)==IOUtils::nodata || rhs(ii)==IOUtils::nodata)
				operator()(ii) = IOUtils::nodata;
			else
				operator()(ii) /= rhs(ii);
		}
	}

	return *this;
}

template<class T> const Array1D<T> Array1D<T>::operator/(const Array1D<T>& rhs)
{
	Array1D<T> result = *this; //make a copy
	result /= rhs; //already implemented

	return result;
}

template<class T> Array1D<T>& Array1D<T>::operator/=(const T& rhs)
{
	//Divide every single member of the Array1D<T>
	if(keep_nodata==false) {
		for (unsigned int ii=0; ii<nx; ii++) {
			operator()(ii) /= rhs;
		}
	} else {
		for (unsigned int ii=0; ii<nx; ii++) {
			if(operator()(ii)!=IOUtils::nodata)
				operator()(ii) /= rhs;
		}
	}

	return *this;
}

template<class T> const Array1D<T> Array1D<T>::operator/(const T& rhs)
{
	Array1D<T> result = *this;
	result /= rhs; //already implemented

	return result;
}

} //end namespace mio

#endif
