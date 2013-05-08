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
#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <meteoio/IOExceptions.h>
#include <meteoio/IOUtils.h>

#include <vector>
#include <limits>
#include <iostream>

namespace mio {

template <class T> class Array2D;

/**
 * @class Array2DProxy
 * @brief The template class Array2DProxy is a helper class for the template class Array2D
 *        with the purpose of adding the [][] operator to Array2D
 *
 * @author Thomas Egger
 */
template <class T> class Array2DProxy {
	public:
		friend class Array2D<T>;
		T& operator[](const unsigned int& j) {
			return array2D(anx, j);
		}

	private:
		Array2DProxy(Array2D<T>& i_array2D, const unsigned int& i_anx) : array2D(i_array2D), anx(i_anx){}
		Array2D<T>& array2D;
		const unsigned int anx;
};

/**
 * @class Array2D
 * @brief The template class Array2D is a 2D Array (Matrix) able to hold any type of object as datatype.
 * It relies on the Array2DProxy class to provide the [][] operator (slower than the (i,j) call).
 * If the compilation flag NOSAFECHECKS is used, bounds check is turned off (leading to increased performances).
 *
 * @ingroup data_str
 * @author Thomas Egger
 */
template<class T> class Array2D {
	public:
		Array2D();

		/**
		* A constructor that creates an array of a given size
		* @param anx number of columns of the new array
		* @param any number of rows of the new array
		*/
		Array2D(const unsigned int& anx, const unsigned int& any);

		/**
		* A constructor that creates an array filled with constant values
		* @param anx number of columns of the new array
		* @param any number of rows of the new array
		* @param init initial value to fill the array with
		*/
		Array2D(const unsigned int& anx, const unsigned int& any, const T& init);

		virtual ~Array2D();

		/**
		* A constructor that can be used to create an Array2D object that is contained in the
		* one passed as i_array2D argument. The resulting Array2D object is a by value copy of
		* a subplane of the plane spanned by the i_array2D
		* @param i_array2D array containing to extract the values from
		* @param i_nx lower left corner cell X index
		* @param i_ny lower left corner cell Y index
		* @param i_ncols number of columns of the new array
		* @param i_nrows number of rows of the new array
		*/
		Array2D(const Array2D<T>& i_array2D, const unsigned int& i_nx, const unsigned int& i_ny,
		        const unsigned int& i_ncols, const unsigned int& i_nrows);

		/**
		* @brief A method that can be used to cut out a subplane of an existing Array2D object
		* that is passed as i_array2D argument. The resulting Array2D object is a by value copy of
		* a subplane of the plane spanned by the i_array2D
		* @param i_array2D array containing to extract the values from
		* @param i_nx lower left corner cell X index
		* @param i_ny lower left corner cell Y index
		* @param i_ncols number of columns of the new array
		* @param i_nrows number of rows of the new array
		*/
		void subset(const Array2D<T>& i_array2D, const unsigned int& i_nx, const unsigned int& i_ny,
		            const unsigned int& i_ncols, const unsigned int& i_nrows);

		/**
		* @brief A method that can be used to insert a subplane into an existing Array2D object
		* that is passed as i_array2D argument. This is exactly the opposite of the subset method
		* an can be used to rebuild an array from subsets.
		* @param i_array2D array containing to subset values
		* @param i_nx lower left corner cell X index
		* @param i_ny lower left corner cell Y index
		* @param i_ncols number of columns of the new array
		* @param i_nrows number of rows of the new array
		*/
		void fill(const Array2D<T>& i_array2D, const unsigned int& i_nx, const unsigned int& i_ny,
		          const unsigned int& i_ncols, const unsigned int& i_nrows);

		void fill(const Array2D<T>& i_array2D, const unsigned int& i_nx, const unsigned int& i_ny);

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

		void resize(const unsigned int& nx, const unsigned int& ny);
		void resize(const unsigned int& nx, const unsigned int& ny, const T& init);
		void size(unsigned int& nx, unsigned int& ny) const;
		unsigned int getNx() const;
		unsigned int getNy() const;

		void clear();
		bool isEmpty() const;

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
		const Array2D<T> getAbs() const;
		void abs();

		const std::string toString() const;
		template<class P> friend std::iostream& operator<<(std::iostream& os, const Array2D<P>& array);
		template<class P> friend std::iostream& operator>>(std::iostream& is, Array2D<P>& array);

		bool checkEpsilonEquality(const Array2D<double>& rhs, const double& epsilon) const;
		static bool checkEpsilonEquality(const Array2D<double>& rhs1, const Array2D<double>& rhs2, const double& epsilon);

		T& operator ()(const unsigned int& x, const unsigned int& y);
		const T operator ()(const unsigned int& x, const unsigned int& y) const;
		T& operator ()(const unsigned int& i);
		const T operator ()(const unsigned int& i) const;
		Array2DProxy<T> operator[](const unsigned int& i);

		Array2D<T>& operator =(const Array2D<T>&);
		Array2D<T>& operator =(const T& value);

		Array2D<T>& operator+=(const T& rhs);
		const Array2D<T> operator+(const T& rhs);
		Array2D<T>& operator+=(const Array2D<T>& rhs);
		const Array2D<T> operator+(const Array2D<T>& rhs);

		Array2D<T>& operator-=(const T& rhs);
		const Array2D<T> operator-(const T& rhs);
		Array2D<T>& operator-=(const Array2D<T>& rhs);
		const Array2D<T> operator-(const Array2D<T>& rhs);

		Array2D<T>& operator*=(const T& rhs);
		const Array2D<T> operator*(const T& rhs);
		Array2D<T>& operator*=(const Array2D<T>& rhs);
		const Array2D<T> operator*(const Array2D<T>& rhs);

		Array2D<T>& operator/=(const T& rhs);
		const Array2D<T> operator/(const T& rhs);
		Array2D<T>& operator/=(const Array2D<T>& rhs);
		const Array2D<T> operator/(const Array2D<T>& rhs);

	protected:
		std::vector<T> vecData;
		unsigned int nx;
		unsigned int ny;
		bool keep_nodata;
};

template<class T> inline T& Array2D<T>::operator()(const unsigned int& i) {
#ifndef NOSAFECHECKS
	return vecData.at(i);
#else
	return vecData[i];
#endif
}

template<class T> inline const T Array2D<T>::operator()(const unsigned int& i) const {
#ifndef NOSAFECHECKS
	return vecData.at(i);
#else
	return vecData[i];
#endif
}
template<class T> inline T& Array2D<T>::operator()(const unsigned int& x, const unsigned int& y) {
#ifndef NOSAFECHECKS
	if ((x >= nx) || (y >= ny)) {
		std::stringstream ss;
		ss << "Trying to access array(" << x << "," << y << ")";
		ss << " while array is (" << nx << "," << ny << ")";
		throw IndexOutOfBoundsException(ss.str(), AT);
	}
#endif
	//COLUMN-MAJOR alignment of the vector: fully C-compatible memory layout
	return vecData[x + y*nx];
}

template<class T> inline const T Array2D<T>::operator()(const unsigned int& x, const unsigned int& y) const {
#ifndef NOSAFECHECKS
	if ((x >= nx) || (y >= ny)) {
		std::stringstream ss;
		ss << "Trying to access array(" << x << "," << y << ")";
		ss << " while array is (" << nx << "," << ny << ")";
		throw IndexOutOfBoundsException(ss.str(), AT);
	}
#endif
	return vecData[x + y*nx];
}

template<class T> Array2DProxy<T> Array2D<T>::operator[](const unsigned int& i) {
	return Array2DProxy<T>(*this, i);
}

template<class T> Array2D<T>::Array2D() : vecData(), nx(0), ny(0), keep_nodata(true)
{
}

template<class T> Array2D<T>::~Array2D() { }

template<class T> Array2D<T>::Array2D(const Array2D<T>& i_array2D, const unsigned int& i_nx, const unsigned int& i_ny,
                                      const unsigned int& i_ncols, const unsigned int& i_nrows) :
                                      vecData(i_ncols*i_nrows), nx(i_ncols), ny(i_nrows), keep_nodata(true)
{
	subset(i_array2D, i_nx, i_ny, i_ncols, i_nrows);
}

template<class T> void Array2D<T>::subset(const Array2D<T>& i_array2D, const unsigned int& i_nx, const unsigned int& i_ny,
                                          const unsigned int& i_ncols, const unsigned int& i_nrows)
{
	if (((i_nx+i_ncols) > i_array2D.nx) || ((i_ny+i_nrows) > i_array2D.ny)) {
		std::stringstream ss;
		ss << "Trying to cut an array of size (" << nx << "," << ny << ") ";
		ss << "to size (" << i_ncols << "," << i_nrows << ") starting at (" << i_nx << "," << i_ny << ")";
		throw IndexOutOfBoundsException(ss.str(), AT);
	}

	if ((i_ncols == 0) || (i_nrows == 0)) //the plane to copy has to make sense
		throw IndexOutOfBoundsException("Trying to cut an array into a null sized array!", AT);

	resize(i_ncols, i_nrows); //create new Array2D object
	//Copy by value subspace
	for (unsigned int jj=0; jj<ny; jj++) {
		for (unsigned int ii=0; ii<nx; ii++) {
			operator()(ii,jj) = i_array2D(i_nx+ii, i_ny+jj);
		}
	}
}

template<class T> void Array2D<T>::fill(const Array2D<T>& i_array2D, const unsigned int& i_nx, const unsigned int& i_ny)
{
	unsigned int i_ncols, i_nrows;
	i_array2D.size(i_ncols, i_nrows);
	fill(i_array2D, i_nx, i_ny, i_ncols, i_nrows);
}

template<class T> void Array2D<T>::fill(const Array2D<T>& i_array2D, const unsigned int& i_nx, const unsigned int& i_ny,
                                        const unsigned int& i_ncols, const unsigned int& i_nrows)
{
	if (((i_nx+i_ncols) > nx) || ((i_ny+i_nrows) > ny)) {
		std::stringstream ss;
		ss << "Filling an array of size (" << nx << "," << ny << ") ";
		ss << "with an array of size (" << i_ncols << "," << i_nrows << ") ";
		ss << "starting at (" << i_nx << "," << i_ny << ")";
		throw IndexOutOfBoundsException(ss.str(), AT);
	}

	if ((i_ncols == 0) || (i_nrows == 0)) //the plane to copy has to make sense
		throw IndexOutOfBoundsException("Filling an array with a null sized array!", AT);

	for(unsigned int jj=i_ny; jj<(i_ny+i_nrows); jj++) {
		for(unsigned int ii=i_nx; ii<(i_nx+i_ncols); ii++) {
			const unsigned int ix = ii-i_nx;
			const unsigned int iy = jj-i_ny;
			operator()(ii,jj) = i_array2D(ix, iy);
		}
	}
}

template<class T> Array2D<T>::Array2D(const unsigned int& anx, const unsigned int& any, const T& init) :
                  vecData(anx*any, init), nx(anx), ny(any), keep_nodata(true)
{
	//resize(anx,any,init);
}

template<class T> Array2D<T>::Array2D(const unsigned int& anx, const unsigned int& any) :
                  vecData(anx*any), nx(anx), ny(any), keep_nodata(true)
{
	//resize(anx,any);
}

template<class T> void Array2D<T>::setKeepNodata(const bool i_keep_nodata) {
	keep_nodata = i_keep_nodata;
}

template<class T> bool Array2D<T>::getKeepNodata() {
	return keep_nodata;
}

template<class T> void Array2D<T>::resize(const unsigned int& anx, const unsigned int& any) {
	clear(); //we won't be able to "rescue" old values, so we reset the whole vector
	vecData.resize(anx*any);
	nx = anx;
	ny = any;
}

template<class T> void Array2D<T>::resize(const unsigned int& anx, const unsigned int& any, const T& init) {
	clear(); //we won't be able to "rescue" old values, so we reset the whole vector
	vecData.resize(anx*any, init);
	nx = anx;
	ny = any;
}

template<class T> void Array2D<T>::size(unsigned int& anx, unsigned int& any) const {
	anx=nx;
	any=ny;
}

template<class T> unsigned int Array2D<T>::getNx() const {
	return nx;
}

template<class T> unsigned int Array2D<T>::getNy() const {
	return ny;
}

template<class T> void Array2D<T>::clear() {
	vecData.clear();
	nx=ny=0;
}

template<class T> bool Array2D<T>::isEmpty() const {
	return (nx==0 && ny==0);
}

template<class T> const std::string Array2D<T>::toString() const {
	std::stringstream os;
	os << "<array2d>\n";
	for(unsigned int jj=0; jj<ny; jj++) {
		const unsigned int jnx = jj*nx;
		for (unsigned int ii=0; ii<nx; ii++) {
			os << vecData[ii+jnx] << " "; //COLUMN-MAJOR alignment
		}
		os << "\n";
	}
	os << "</array2d>\n";
	return os.str();
}

template<class P> std::iostream& operator<<(std::iostream& os, const Array2D<P>& array) {
	os.write(reinterpret_cast<const char*>(&array.keep_nodata), sizeof(array.keep_nodata));
	os.write(reinterpret_cast<const char*>(&array.nx), sizeof(array.nx));
	os.write(reinterpret_cast<const char*>(&array.ny), sizeof(array.ny));
	os.write(reinterpret_cast<const char*>(&array.vecData[0]), array.nx*array.ny*sizeof(P));
	return os;
}

template<class P> std::iostream& operator>>(std::iostream& is, Array2D<P>& array) {
	is.read(reinterpret_cast<char*>(&array.keep_nodata), sizeof(array.keep_nodata));
	is.read(reinterpret_cast<char*>(&array.nx), sizeof(array.nx));
	is.read(reinterpret_cast<char*>(&array.ny), sizeof(array.ny));
	array.vecData.resize(array.nx*array.ny);
	is.read(reinterpret_cast<char*>(&array.vecData[0]), array.nx*array.ny*sizeof(P)); //30 times faster than assign() or copy()
	return is;
}


template<class T> T Array2D<T>::getMin() const {

	T min = std::numeric_limits<T>::max();

	const unsigned int nxy = ny*nx;
	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nxy; jj++) {
			const T val = vecData[jj];
			if(val<min) min=val;
		}
		return min;
	} else {
		for (unsigned int jj=0; jj<nxy; jj++) {
			const T val = vecData[jj];
			if(val!=IOUtils::nodata && val<min) min=val;
		}
		if(min!=std::numeric_limits<T>::max()) return min;
		else return (T)IOUtils::nodata;
	}
}

template<class T> T Array2D<T>::getMax() const {

	T max = -std::numeric_limits<T>::max();

	const unsigned int nxy = ny*nx;
	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nxy; jj++) {
			const T val = vecData[jj];
			if(val>max) max=val;
		}
		return max;
	} else {
		for (unsigned int jj=0; jj<nxy; jj++) {
			const T val = vecData[jj];
			if(val!=IOUtils::nodata && val>max) max=val;
		}
		if(max!=-std::numeric_limits<T>::max()) return max;
		else return (T)IOUtils::nodata;
	}
}

template<class T> T Array2D<T>::getMean() const {

	T mean = 0;
	const unsigned int nxy = nx*ny;

	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nxy; jj++) {
			const T val = vecData[jj];
			mean += val;
		}
		if(nxy>0) return mean/(T)(nxy);
		else return (T)0;
	} else {
		unsigned int count = 0;
		for (unsigned int jj=0; jj<nxy; jj++) {
			const T val = vecData[jj];
			if(val!=IOUtils::nodata) {
				mean += val;
				count++;
			}
		}
		if(count>0) return mean/(T)(count);
		else return (T)IOUtils::nodata;
	}
}

template<class T> size_t Array2D<T>::getCount() const
{
	const unsigned int nxy = nx*ny;

	if(keep_nodata==false) {
		return (size_t)nxy;
	} else {
		size_t count = 0;
		for (unsigned int ii=0; ii<nxy; ii++) {
			if(vecData[ii]!=IOUtils::nodata) count++;
		}
		return count;
	}
}

template<class T> void Array2D<T>::abs() {
	if(std::numeric_limits<T>::is_signed) {
		const unsigned int nxy = nx*ny;
		if(keep_nodata==false) {
			for (unsigned int ii=0; ii<nxy; ii++) {
				T& val = vecData[ii];
				if(val<0) val=-val;
			}
		} else {
			for (unsigned int ii=0; ii<nxy; ii++) {
				T& val = vecData[ii];
				if(val<0 && val!=IOUtils::nodata) val=-val;
			}
		}
	}
}

template<class T> const Array2D<T> Array2D<T>::getAbs() const {
	Array2D<T> result = *this; //make a copy
	result.abs(); //already implemented

	return result;
}


//arithmetic operators
template<class T> bool Array2D<T>::checkEpsilonEquality(const Array2D<double>& rhs, const double& epsilon) const {
	if(nx!=rhs.nx || ny!=rhs.ny) return false;

	const unsigned int nxy = nx*ny;
	for (unsigned int jj=0; jj<nxy; jj++)
		if(IOUtils::checkEpsilonEquality(vecData[jj], rhs.vecData[jj], epsilon)==false) return false;

	return true;
}

template<class T> bool Array2D<T>::checkEpsilonEquality(const Array2D<double>& rhs1, const Array2D<double>& rhs2, const double& epsilon) {
	return rhs1.checkEpsilonEquality(rhs2, epsilon);
}

template<class T> Array2D<T>& Array2D<T>::operator=(const Array2D<T>& source) {
	if(this != &source) {
		keep_nodata = source.keep_nodata;
		nx = source.nx;
		ny = source.ny;
		vecData = source.vecData;
	}
	return *this;
}

template<class T> Array2D<T>& Array2D<T>::operator=(const T& value) {
	std::fill(vecData.begin(), vecData.end(), value);
	return *this;
}

template<class T> Array2D<T>& Array2D<T>::operator+=(const Array2D<T>& rhs)
{
	//They have to have equal size
	if ((rhs.nx != nx) || (rhs.ny != ny)) {
		std::stringstream ss;
		ss << "Trying to add two Array2D objects with different dimensions: ";
		ss << "(" << nx << "," << ny << ") + (" << rhs.nx << "," << rhs.ny << ")";
		throw IOException(ss.str(), AT);
	}

	const unsigned int nxy = nx*ny;
	//Add to every single member of the Array2D<T>
	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nxy; jj++)
			vecData[jj] += rhs(jj);
	} else {
		for (unsigned int jj=0; jj<nxy; jj++) {
			if(vecData[jj]==IOUtils::nodata || rhs(jj)==IOUtils::nodata)
				vecData[jj] = IOUtils::nodata;
			else
				vecData[jj] += rhs(jj);
		}
	}

	return *this;
}

template<class T> const Array2D<T> Array2D<T>::operator+(const Array2D<T>& rhs)
{
	Array2D<T> result = *this; //make a copy
	result += rhs; //already implemented

	return result;
}

template<class T> Array2D<T>& Array2D<T>::operator+=(const T& rhs)
{
	//Add to every single member of the Array2D<T>
	const unsigned int nxy = nx*ny;

	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nxy; jj++)
			vecData[jj] += rhs;
	} else {
		for (unsigned int jj=0; jj<nxy; jj++) {
			if(vecData[jj]!=IOUtils::nodata)
				vecData[jj] += rhs;
		}
	}

	return *this;
}

template<class T> const Array2D<T> Array2D<T>::operator+(const T& rhs)
{
	Array2D<T> result = *this;
	result += rhs; //already implemented

	return result;
}

template<class T> Array2D<T>& Array2D<T>::operator-=(const Array2D<T>& rhs)
{
	//They have to have equal size
	if ((rhs.nx != nx) || (rhs.ny != ny)){
		std::stringstream ss;
		ss << "Trying to substract two Array2D objects with different dimensions: ";
		ss << "(" << nx << "," << ny << ") - (" << rhs.nx << "," << rhs.ny << ")";
		throw IOException(ss.str(), AT);
	}
	//Substract to every single member of the Array2D<T>
	const unsigned int nxy = nx*ny;

	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nxy; jj++)
			vecData[jj] -= rhs(jj);
	} else {
		for (unsigned int jj=0; jj<nxy; jj++) {
			if(vecData[jj]==IOUtils::nodata || rhs(jj)==IOUtils::nodata)
				vecData[jj] = IOUtils::nodata;
			else
				vecData[jj] -= rhs(jj);
		}
	}

	return *this;
}

template<class T> const Array2D<T> Array2D<T>::operator-(const Array2D<T>& rhs)
{
	Array2D<T> result = *this; //make a copy
	result -= rhs; //already implemented

	return result;
}

template<class T> Array2D<T>& Array2D<T>::operator-=(const T& rhs)
{
	*this += -rhs;
	return *this;
}

template<class T> const Array2D<T> Array2D<T>::operator-(const T& rhs)
{
	Array2D<T> result = *this;
	result += -rhs; //already implemented

	return result;
}

template<class T> Array2D<T>& Array2D<T>::operator*=(const Array2D<T>& rhs)
{
	//They have to have equal size
	if ((rhs.nx != nx) || (rhs.ny != ny)){
		std::stringstream ss;
		ss << "Trying to multiply two Array2D objects with different dimensions: ";
		ss << "(" << nx << "," << ny << ") * (" << rhs.nx << "," << rhs.ny << ")";
		throw IOException(ss.str(), AT);
	}
	//Add to every single member of the Array2D<T>
	const unsigned int nxy = nx*ny;

	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nxy; jj++)
			vecData[jj] *= rhs(jj);
	} else {
		for (unsigned int jj=0; jj<nxy; jj++) {
			if(vecData[jj]==IOUtils::nodata || rhs(jj)==IOUtils::nodata)
				vecData[jj] = IOUtils::nodata;
			else
				vecData[jj] *= rhs(jj);
		}
	}

	return *this;
}

template<class T> const Array2D<T> Array2D<T>::operator*(const Array2D<T>& rhs)
{
	Array2D<T> result = *this; //make a copy
	result *= rhs; //already implemented

	return result;
}

template<class T> Array2D<T>& Array2D<T>::operator*=(const T& rhs)
{
	//Multiply to every single member of the Array2D<T>
	const unsigned int nxy = nx*ny;

	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nxy; jj++)
			vecData[jj] *= rhs;
	} else {
		for (unsigned int jj=0; jj<nxy; jj++) {
			if(vecData[jj]!=IOUtils::nodata)
				vecData[jj] *= rhs;
		}
	}

	return *this;
}

template<class T> const Array2D<T> Array2D<T>::operator*(const T& rhs)
{
	Array2D<T> result = *this;
	result *= rhs; //already implemented

	return result;
}

template<class T> Array2D<T>& Array2D<T>::operator/=(const Array2D<T>& rhs)
{
	//They have to have equal size
	if ((rhs.nx != nx) || (rhs.ny != ny)){
		std::stringstream ss;
		ss << "Trying to divide two Array2D objects with different dimensions: ";
		ss << "(" << nx << "," << ny << ") / (" << rhs.nx << "," << rhs.ny << ")";
		throw IOException(ss.str(), AT);
	}
	//Divide every single member of the Array2D<T>
	const unsigned int nxy = nx*ny;

	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nxy; jj++)
			vecData[jj] /= rhs(jj);
	} else {
		for (unsigned int jj=0; jj<nxy; jj++) {
			if(vecData[jj]==IOUtils::nodata || rhs(jj)==IOUtils::nodata)
				vecData[jj] = IOUtils::nodata;
			else
				vecData[jj] /= rhs(jj);
		}
	}

	return *this;
}

template<class T> const Array2D<T> Array2D<T>::operator/(const Array2D<T>& rhs)
{
	Array2D<T> result = *this; //make a copy
	result /= rhs; //already implemented

	return result;
}

template<class T> Array2D<T>& Array2D<T>::operator/=(const T& rhs)
{
	*this *= (1./rhs);
	return *this;
}

template<class T> const Array2D<T> Array2D<T>::operator/(const T& rhs)
{
	Array2D<T> result = *this;
	result *= (1./rhs); //already implemented

	return result;
}

} //end namespace mio

#endif
