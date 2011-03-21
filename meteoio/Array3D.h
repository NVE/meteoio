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
#ifndef ARRAY3D_H
#define ARRAY3D_H

#include <meteoio/IOUtils.h>
#include <meteoio/IOExceptions.h>

#include <vector>
#include <limits>
#include <iostream>

#define NOSAFECHECKS

namespace mio {

template <class T> class Array3D;
template <class T> class Array3DProxy2;

/**
 * @class Array3DProxy
 * @brief The template class Array3DProxy is a helper class for the template class Array3D
 *        with the purpose of adding the [][] operator to Array3D
 *
 * @author Thomas Egger
 */
template <class T> class Array3DProxy {
 	public:
		friend class Array3D<T>;
		Array3DProxy2<T> operator[](const unsigned int& _any) {
			return Array3DProxy2<T>(array3D, anx, _any); 
		}

 	private:
 		Array3DProxy(Array3D<T>& _array3D, const unsigned int& _anx) : array3D(_array3D), anx(_anx){}
		Array3D<T>& array3D;
		const unsigned int anx;
};

/**
 * @class Array3DProxy2
 * @brief The template class Array3DProxy2 is a helper class for the template class Array3D
 *        with the purpose of adding the [][][] operator to Array3D
 *
 * @author Thomas Egger
 */
template <class T> class Array3DProxy2 {
 	public:
		friend class Array3DProxy<T>;
		T& operator[](const unsigned int& _anz) {
			return array3D(anx, any, _anz);
		}

	private:
 		Array3DProxy2(Array3D<T>& _array3D, const unsigned int& _anx, 
				    const unsigned int& _any) : array3D(_array3D), anx(_anx), any(_any){}
		Array3D<T>& array3D;
		const unsigned int anx;
		const unsigned int any;
}; 


/**
 * @class Array3D
 * @brief The template class Array3D is a 3D Array (Tensor) able to hold any type of object as datatype. 
 * It relies on the Array3DProxy2 class to provide the [][][] operator (slower than the (i,j,k) call).
 * @ingroup data_str
 * @date  2009-07-19
 * @author Thomas Egger
 */
template<class T> class Array3D {
	public:
		Array3D();

		/**
		* A constructor that can be used to create an Array3D object that is contained in the
		* one passed as _array3D argument. The resulting Array3D object is a by value copy of
		* a subvolume of the volume spanned by the _array3D
		* @param _array3D array containing to extract the values from
		* @param _nx lower left corner cell X index
		* @param _ny lower left corner cell Y index
		* @param _nz lower left corner cell Z index
		* @param _ncols number of columns of the new array
		* @param _nrows number of rows of the new array
		* @param _ndepth number of depths of the new array
		*/
		Array3D(const Array3D<T>& _array3D,
		        const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz,
		        const unsigned int& _ncols, const unsigned int& _nrows, const unsigned int& _ndepth);

		/**
		* A constructor that creates an array of a given size
		* @param anx number of columns of the new array
		* @param any number of rows of the new array
		* @param anz number of rows of the new array
		*/
		Array3D(const unsigned int& anx, const unsigned int& any, const unsigned int& anz);

		/**
		* A constructor that creates an array filled with constant values
		* @param anx number of columns of the new array
		* @param any number of rows of the new array
		* @param anz number of depths of the new array
		* @param init initial value to fill the array with
		*/
		Array3D(const unsigned int& anx, const unsigned int& any, const unsigned int& anz, const T& init);

		/**
		* A method that can be used to create an Array3D object that is contained in the
		* one passed as _array3D argument. The resulting Array3D object is a by value copy of
		* a subvolume of the volume spanned by the _array3D
		* @param _array3D array containing to extract the values from
		* @param _nx lower left corner cell X index
		* @param _ny lower left corner cell Y index
		* @param _nz lower left corner cell Z index
		* @param _ncols number of columns of the new array
		* @param _nrows number of rows of the new array
		* @param _ndepth number of depths of the new array
		*/
		void subset(const Array3D<T>& _array3D,
		            const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz,
		            const unsigned int& _ncols, const unsigned int& _nrows, const unsigned int& _ndepth);

		void resize(const unsigned int& anx, const unsigned int& any, const unsigned int& anz);
		void resize(const unsigned int& anx, const unsigned int& any, const unsigned int& anz, const T& init);
		void size(unsigned int& anx, unsigned int& any, unsigned int& anz) const;
		void clear();

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

		template<class P> friend std::ostream& operator<<(std::ostream& os, const Array3D<P>& array);
		
		T& operator ()(const unsigned int& i);
		const T operator ()(const unsigned int& i) const;
		T& operator ()(const unsigned int& x, const unsigned int& y, const unsigned int& z);
		const T operator ()(const unsigned int& x, const unsigned int& y, const unsigned int& z) const;
		Array3DProxy<T> operator[](const unsigned int& i);

		Array3D<T>& operator =(const Array3D<T>&);
		Array3D<T>& operator =(const T& value);
		
		Array3D<T>& operator+=(const T& rhs);
		const Array3D<T> operator+(const T& rhs);
		Array3D<T>& operator+=(const Array3D<T>& rhs);
		const Array3D<T> operator+(const Array3D<T>& rhs);

		Array3D<T>& operator-=(const T& rhs);
		const Array3D<T> operator-(const T& rhs);
		Array3D<T>& operator-=(const Array3D<T>& rhs);
		const Array3D<T> operator-(const Array3D<T>& rhs);

		Array3D<T>& operator*=(const T& rhs);
		const Array3D<T> operator*(const T& rhs);
		Array3D<T>& operator*=(const Array3D<T>& rhs);
		const Array3D<T> operator*(const Array3D<T>& rhs);

		Array3D<T>& operator/=(const T& rhs);
		const Array3D<T> operator/(const T& rhs);
		Array3D<T>& operator/=(const Array3D<T>& rhs);
		const Array3D<T> operator/(const Array3D<T>& rhs);

		void abs();
		const Array3D<T> getAbs() const;

	protected:
		std::vector<T> vecData; ///< The actual objects are stored in a one-dimensional vector
		unsigned int nx;
		unsigned int ny;
		unsigned int nz;
		unsigned int nxny; //nx times ny
};

template<class T> T& Array3D<T>::operator()(const unsigned int& i) {
	return vecData[i];
}

template<class T> const T Array3D<T>::operator()(const unsigned int& i) const {
	return vecData[i];
}

template<class T> T& Array3D<T>::operator()(const unsigned int& x, const unsigned int& y, const unsigned int& z) {
#ifndef NOSAFECHECKS
	if ((x >= nx) || (y >= ny) || (z >= nz)) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif

	//ROW-MAJOR alignment of the vector: fully C-compatible memory layout
	return vecData[x + y*nx + z*nxny];
}

template<class T> const T Array3D<T>::operator()(const unsigned int& x, const unsigned int& y, const unsigned int& z) const {
#ifndef NOSAFECHECKS
	if ((x >= nx) || (y >= ny) || (z >= nz)) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif
	return vecData[x + y*nx + z*nxny];
}

template<class T> Array3DProxy<T> Array3D<T>::operator[](const unsigned int& i) {
	return Array3DProxy<T>(*this, i); 
}


template<class T> Array3D<T>::Array3D() {
	nx = ny = nz = nxny = 0;
}

template<class T> Array3D<T>::Array3D(const Array3D<T>& _array3D,
                                      const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz,
                                      const unsigned int& _ncols, const unsigned int& _nrows, const unsigned int& _ndepth)
{
	subset(_array3D, _nx, _ny, _nz, _ncols, _nrows, _ndepth);
}

template<class T> void Array3D<T>::subset(const Array3D<T>& _array3D,
                                     const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz,
                                     const unsigned int& _ncols, const unsigned int& _nrows, const unsigned int& _ndepth)
{
	
	if (((_nx+_ncols) > _array3D.nx) || ((_ny+_nrows) > _array3D.ny) || ((_nz+_ndepth) > _array3D.nz))
		throw IndexOutOfBoundsException("", AT);

	if ((_ncols == 0) || (_nrows == 0) || (_ndepth == 0)) //the space has to make sense
		throw IndexOutOfBoundsException("", AT);

	resize(_ncols, _nrows, _ndepth); //create new Array3D object

	//Copy by value subspace
	for (unsigned int ii=0; ii<nz; ii++) { 
		for (unsigned int jj=0; jj<ny; jj++) {
			for (unsigned int kk=0; kk<nx; kk++) {
				//Running through the vector in order of memory alignment
				operator()(kk,jj,ii) = _array3D(_nx+kk, _ny+jj, _nz+ii); 
			}
		}
	}
}

template<class T> Array3D<T>::Array3D(const unsigned int& anx, const unsigned int& any, const unsigned int& anz) {
	resize(anx, any, anz);
}

template<class T> Array3D<T>::Array3D(const unsigned int& anx, const unsigned int& any, const unsigned int& anz, const T& init) {
	resize(anx, any, anz, init);
}

template<class T> void Array3D<T>::resize(const unsigned int& anx, const unsigned int& any, const unsigned int& anz) {
	clear();  //we won't be able to "rescue" old values, so we reset the whole vector
	vecData.resize(anx*any*anz);
	nx = anx;
	ny = any;
	nz = anz;
	nxny = nx*ny;
}

template<class T> void Array3D<T>::resize(const unsigned int& anx, const unsigned int& any, const unsigned int& anz, const T& init) {
	resize(anx, any, anz);
	std::fill(vecData.begin(), vecData.end(), init);
}

template<class T> void Array3D<T>::size(unsigned int& anx, unsigned int& any, unsigned int& anz) const {
	anx=nx;
	any=ny;
	anz=nz;
}

template<class T> void Array3D<T>::clear() {
	vecData.clear();
	nx = ny = nz = nxny = 0;
}

template<class T> std::ostream& operator<<(std::ostream& os, const Array3D<T>& array) {
	os << "<array3d>\n";
	for (unsigned int kk=0; kk<array.nz; kk++) {
		os << "depth[" << kk << "]\n";
		for(unsigned int ii=0; ii<array.nx; ii++) {
			for (unsigned int jj=0; jj<array.ny; jj++) {
				os << array(ii,jj,kk) << " ";
			}
			os << "\n";
		}
	}
	os << "</array3d>\n";
	return os;
}

template<class T> T Array3D<T>::getMin(const IOUtils::nodata_handling flag_nodata) const {

	T min = std::numeric_limits<T>::max();
	const unsigned int nxyz = ny*nx*nz;
	
	if(flag_nodata==IOUtils::RAW_NODATA) {
		for (unsigned int jj=0; jj<nxyz; jj++) {
			const T val = operator()(jj);
			if(val<min) min=val;
		}
		return min;
	} else if(flag_nodata==IOUtils::PARSE_NODATA) {
		for (unsigned int jj=0; jj<nxyz; jj++) {
			const T val = operator()(jj);
			if(val!=IOUtils::nodata && val<min) min=val;
		}
		if(min!=std::numeric_limits<T>::max()) return min;
		else return (T)IOUtils::nodata;
	} else {
		throw InvalidArgumentException("Unknown nodata_handling flag",AT);
	}
}

template<class T> T Array3D<T>::getMax(const IOUtils::nodata_handling flag_nodata) const {

	T max = -std::numeric_limits<T>::max();
	const unsigned int nxyz = ny*nx*nz;
	
	if(flag_nodata==IOUtils::RAW_NODATA) {
		for (unsigned int jj=0; jj<nxyz; jj++) {
			const T val = operator()(jj);
			if(val>max) max=val;
		}
		return max;
	} else if(flag_nodata==IOUtils::PARSE_NODATA) {
		for (unsigned int jj=0; jj<nxyz; jj++) {
			const T val = operator()(jj);
			if(val!=IOUtils::nodata && val>max) max=val;
		}
		if(max!=-std::numeric_limits<T>::max()) return max;
		else return (T)IOUtils::nodata;
	} else {
		throw InvalidArgumentException("Unknown nodata_handling flag",AT);
	}
}

template<class T> T Array3D<T>::getMean(const IOUtils::nodata_handling flag_nodata) const {

	T mean = 0;
	const unsigned int nxyz = nx*ny*nz;
	
	if(flag_nodata==IOUtils::RAW_NODATA) {
		for (unsigned int jj=0; jj<nxyz; jj++) {
			const T val = operator()(jj);
			mean += val;
		}
		if(nxyz>0) return mean/(T)(nxyz);
		else return (T)0;
	} else if(flag_nodata==IOUtils::PARSE_NODATA) {
		unsigned int count = 0;
		for (unsigned int jj=0; jj<nxyz; jj++) {
			const T val = operator()(jj);
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

//arithmetic operators
template<class T> Array3D<T>& Array3D<T>::operator=(const Array3D<T>& source) {
	if(this != &source) {
		nx = source.nx;
		ny = source.ny;
		nz = source.nz;
		nxny = source.nxny;
		vecData = source.vecData;
	}
	return *this;
}

template<class T> Array3D<T>& Array3D<T>::operator=(const T& value) {
	std::fill(vecData.begin(), vecData.end(), value);
	return *this;
}

template<class T> Array3D<T>& Array3D<T>::operator+=(const Array3D<T>& rhs)
{
	//They have to have equal size
	if ((rhs.nx != nx) || (rhs.ny != ny) || (rhs.nz != nz))
		throw IOException("Trying to add two Array3D objects with different dimensions", AT);

	//Add to every single member of the Array3D<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			for(unsigned int kk=0; kk<nz; kk++) {
				operator()(ii,jj,kk) += rhs(ii,jj,kk);
			}
		}
	}

	return *this;
}

template<class T> const Array3D<T> Array3D<T>::operator+(const Array3D<T>& rhs)
{
	Array3D<T> result = *this; //make a copy
	result += rhs; //already implemented

	return result;
}

template<class T> Array3D<T>& Array3D<T>::operator+=(const T& rhs)
{
	//Add to every single member of the Array3D<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			for(unsigned int kk=0; kk<nz; kk++) {
				operator()(ii,jj,kk) += rhs;
			}
		}
	}

	return *this;
}

template<class T> const Array3D<T> Array3D<T>::operator+(const T& rhs)
{
	Array3D<T> result = *this;
	result += rhs; //already implemented

	return result;
}

template<class T> Array3D<T>& Array3D<T>::operator-=(const Array3D<T>& rhs)
{
	//They have to have equal size
	if ((rhs.nx != nx) || (rhs.ny != ny) || (rhs.nz != nz))
		throw IOException("Trying to substract two Array3D objects with different dimensions", AT);

	//Substract to every single member of the Array3D<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			for(unsigned int kk=0; kk<nz; kk++) {
				operator()(ii,jj,kk) -= rhs(ii,jj,kk);
			}
		}
	}

	return *this;
}

template<class T> const Array3D<T> Array3D<T>::operator-(const Array3D<T>& rhs)
{
	Array3D<T> result = *this; //make a copy
	result -= rhs; //already implemented

	return result;
}

template<class T> Array3D<T>& Array3D<T>::operator-=(const T& rhs)
{
	//Substract to every single member of the Array3D<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			for(unsigned int kk=0; kk<nz; kk++) {
				operator()(ii,jj,kk) -= rhs;
			}
		}
	}

	return *this;
}

template<class T> const Array3D<T> Array3D<T>::operator-(const T& rhs)
{
	Array3D<T> result = *this;
	result -= rhs; //already implemented

	return result;
}

template<class T> Array3D<T>& Array3D<T>::operator*=(const Array3D<T>& rhs)
{
	//They have to have equal size
	if ((rhs.nx != nx) || (rhs.ny != ny) || (rhs.nz != nz))
		throw IOException("Trying to multiply two Array3D objects with different dimensions", AT);

	//Multiply every single member of the Array3D<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			for(unsigned int kk=0; kk<nz; kk++) {
				operator()(ii,jj,kk) *= rhs(ii,jj,kk);
			}
		}
	}

	return *this;
}

template<class T> const Array3D<T> Array3D<T>::operator*(const Array3D<T>& rhs)
{
	Array3D<T> result = *this; //make a copy
	result *= rhs; //already implemented

	return result;
}

template<class T> Array3D<T>& Array3D<T>::operator*=(const T& rhs)
{
	//Multiply every single member of the Array3D<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			for(unsigned int kk=0; kk<nz; kk++) {
				operator()(ii,jj,kk) *= rhs;
			}
		}
	}

	return *this;
}

template<class T> const Array3D<T> Array3D<T>::operator*(const T& rhs)
{
	Array3D<T> result = *this;
	result *= rhs; //already implemented

	return result;
}

template<class T> Array3D<T>& Array3D<T>::operator/=(const Array3D<T>& rhs)
{
	//They have to have equal size
	if ((rhs.nx != nx) || (rhs.ny != ny) || (rhs.nz != nz))
		throw IOException("Trying to divide two Array3D objects with different dimensions", AT);

	//Divide every single member of the Array3D<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			for(unsigned int kk=0; kk<nz; kk++) {
				operator()(ii,jj,kk) /= rhs(ii,jj,kk);
			}
		}
	}

	return *this;
}

template<class T> const Array3D<T> Array3D<T>::operator/(const Array3D<T>& rhs)
{
	Array3D<T> result = *this; //make a copy
	result /= rhs; //already implemented

	return result;
}

template<class T> Array3D<T>& Array3D<T>::operator/=(const T& rhs)
{
	//Divide every single member of the Array3D<T>
	for (unsigned int ii=0; ii<nx; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			for(unsigned int kk=0; kk<nz; kk++) {
				operator()(ii,jj,kk) /= rhs;
			}
		}
	}

	return *this;
}

template<class T> const Array3D<T> Array3D<T>::operator/(const T& rhs)
{
	Array3D<T> result = *this;
	result /= rhs; //already implemented

	return result;
}

template<class T> void Array3D<T>::abs() {
	if(std::numeric_limits<T>::is_signed) {
		const unsigned int nxyz = nx*ny*nz;
		for (unsigned int ii=0; ii<nxyz; ii++) {
			T& val = operator()(ii);
			if(val<0) val=-val;
		}
	}
}


template<class T> const Array3D<T> Array3D<T>::getAbs() const {
	Array3D<T> result = *this; //make a copy
	result.abs(); //already implemented

	return result;
}

} //end namespace mio

#endif
