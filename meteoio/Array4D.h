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
#ifndef ARRAY4D_H
#define ARRAY4D_H

#include <meteoio/IOUtils.h>
#include <meteoio/IOExceptions.h>

#include <vector>
#include <limits>
#include <iostream>

namespace mio {

/**
 * @class Array4D
 * @brief The template class Array4D is a 4D Array (Tensor) able to hold any type of object as datatype.
 * @date  2012-12-20
 * @author Till Schumann
 */
template<class T> class Array4D {
	public:
		Array4D();

		/**
		* A constructor that can be used to create an Array4D object that is contained in the
		* one passed as i_array4D argument. The resulting Array4D object is a by value copy of
		* a subvolume of the volume spanned by the i_array4D
		* @param i_array4D array containing to extract the values from
		* @param i_nw lower left corner cell W index
		* @param i_nx lower left corner cell X index
		* @param i_ny lower left corner cell Y index
		* @param i_nz lower left corner cell Z index
		* @param i_sizeW length of the w dimension of the new array
		* @param i_sizeX length of the x dimension of the new array
		* @param i_sizeY length of the y dimension of the new array
		* @param i_sizeZ length of the z dimension of the new array
		*/
		Array4D(const Array4D<T>& i_array4D,
		        const unsigned int& i_nw, const unsigned int& i_nx, const unsigned int& i_ny, const unsigned int& i_nz,
		        const unsigned int& i_sizeW, const unsigned int& i_sizeX, const unsigned int& i_sizeY, const unsigned int& i_sizeZ);

		Array4D(const unsigned int& anw, const unsigned int& anx, const unsigned int& any, const unsigned int& anz);

		/**
		* A constructor that creates an array filled with constant values
		* @param i_sizeW length of the w dimension of the new array
		* @param i_sizeX length of the x dimension of the new array
		* @param i_sizeY length of the y dimension of the new array
		* @param i_sizeZ length of the z dimension of the new array
		* @param init initial value to fill the array with
		*/
		Array4D(const unsigned int& i_sizeW, const unsigned int& i_sizeX, const unsigned int& i_sizeY, const unsigned int& i_sizeZ, const T& init);

		/**
		* A method that can be used to create an Array4D object that is contained in the
		* one passed as i_array4D argument. The resulting Array4D object is a by value copy of
		* a subvolume of the volume spanned by the i_array4D
		* @param i_array4D array containing to extract the values from
		* @param i_nw lower left corner cell W index
		* @param i_nx lower left corner cell X index
		* @param i_ny lower left corner cell Y index
		* @param i_nz lower left corner cell Z index
		* @param i_sizeW length of the w dimension of the new array
		* @param i_sizeX length of the x dimension of the new array
		* @param i_sizeY length of the y dimension of the new array
		* @param i_sizeZ length of the z dimension of the new array
		*/
		void subset(const Array4D<T>& i_array4D,
		            const unsigned int& i_nw, const unsigned int& i_nx, const unsigned int& i_ny, const unsigned int& i_nz,
		            const unsigned int& i_sizeW, const unsigned int& i_sizeX, const unsigned int& i_sizeY, const unsigned int& i_sizeZ);

		/**
		* @brief A method that can be used to insert a subplane into an existing Array4D object
		* that is passed as i_array4D argument. This is exactly the opposite of the subset method
		* an can be used to rebuild an array from subsets.
		* @param i_array4D array containing to extract the values from
		* @param i_nw lower left corner cell W index
		* @param i_nx lower left corner cell X index
		* @param i_ny lower left corner cell Y index
		* @param i_nz lower left corner cell Z index
		* @param i_sizeW length of the w dimension of the new array
		* @param i_sizeX length of the x dimension of the new array
		* @param i_sizeY length of the y dimension of the new array
		* @param i_sizeZ length of the z dimension of the new array
		*/
		void fill(const Array4D<T>& i_array4D,
		          const unsigned int& i_nw, const unsigned int& i_nx, const unsigned int& i_ny, const unsigned int& i_nz,
		          const unsigned int& i_sizeW, const unsigned int& i_sizeX, const unsigned int& i_sizeY, const unsigned int& i_sizeZ);

		void fill(const Array4D<T>& i_array4D, const unsigned int& i_nw, const unsigned int& i_nx, const unsigned int& i_ny, const unsigned int& i_nz);

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

		void resize(const unsigned int& anw, const unsigned int& anx, const unsigned int& any, const unsigned int& anz);
		void resize(const unsigned int& anw, const unsigned int& anx, const unsigned int& any, const unsigned int& anz, const T& init);
		void size(unsigned int& anw, unsigned int& anx, unsigned int& any, unsigned int& anz) const;
		unsigned int getNw() const;
		unsigned int getNx() const;
		unsigned int getNy() const;
		unsigned int getNz() const;
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
		const Array4D<T> getAbs() const;
		void abs();

		const std::string toString() const;
		template<class P> friend std::iostream& operator<<(std::iostream& os, const Array4D<P>& array);
		template<class P> friend std::iostream& operator>>(std::iostream& is, Array4D<P>& array);

		bool checkEpsilonEquality(const Array4D<double>& rhs, const double& epsilon) const;
		static bool checkEpsilonEquality(const Array4D<double>& rhs1, const Array4D<double>& rhs2, const double& epsilon);

		T& operator ()(const unsigned int& i);
		const T operator ()(const unsigned int& i) const;
		T& operator ()(const unsigned int& w, const unsigned int& x, const unsigned int& y, const unsigned int& z);
		const T operator ()(const unsigned int& w, const unsigned int& x, const unsigned int& y, const unsigned int& z) const;

		Array4D<T>& operator =(const Array4D<T>&);
		Array4D<T>& operator =(const T& value);

		Array4D<T>& operator+=(const T& rhs);
		const Array4D<T> operator+(const T& rhs);
		Array4D<T>& operator+=(const Array4D<T>& rhs);
		const Array4D<T> operator+(const Array4D<T>& rhs);

		Array4D<T>& operator-=(const T& rhs);
		const Array4D<T> operator-(const T& rhs);
		Array4D<T>& operator-=(const Array4D<T>& rhs);
		const Array4D<T> operator-(const Array4D<T>& rhs);

		Array4D<T>& operator*=(const T& rhs);
		const Array4D<T> operator*(const T& rhs);
		Array4D<T>& operator*=(const Array4D<T>& rhs);
		const Array4D<T> operator*(const Array4D<T>& rhs);

		Array4D<T>& operator/=(const T& rhs);
		const Array4D<T> operator/(const T& rhs);
		Array4D<T>& operator/=(const Array4D<T>& rhs);
		const Array4D<T> operator/(const Array4D<T>& rhs);

	protected:
		std::vector<T> vecData; ///< The actual objects are stored in a one-dimensional vector
		unsigned int nw;
		unsigned int nx;
		unsigned int ny;
		unsigned int nz;
		unsigned int nwnx; //xw time nx
		unsigned int nwnxny; //xw time nx times ny
		bool keep_nodata;
};

template<class T> inline T& Array4D<T>::operator()(const unsigned int& i) {
#ifndef NOSAFECHECKS
	return vecData.at(i);
#else
	return vecData[i];
#endif
}

template<class T> inline const T Array4D<T>::operator()(const unsigned int& i) const {
#ifndef NOSAFECHECKS
	return vecData.at(i);
#else
	return vecData[i];
#endif
}

template<class T> inline T& Array4D<T>::operator()(const unsigned int& w, const unsigned int& x, const unsigned int& y, const unsigned int& z) {
#ifndef NOSAFECHECKS
	if ((w >= nw) || (x >= nx) || (y >= ny) || (z >= nz))  {
		std::stringstream ss;
		ss << "Trying to access array(" << w << "," << x << "," << y << "," << z << ")";
		ss << " while array is ("<< nw << "," << nx << "," << ny << "," << nz << ")";
		throw IndexOutOfBoundsException(ss.str(), AT);
	}
#endif

	//ROW-MAJOR alignment of the vector: fully C-compatible memory layout
	return vecData[w + x*nw + y*nwnx + z*nwnxny];
}

template<class T> inline const T Array4D<T>::operator()(const unsigned int& w, const unsigned int& x, const unsigned int& y, const unsigned int& z) const {
#ifndef NOSAFECHECKS
	if ((w >= nw) || (x >= nx) || (y >= ny) || (z >= nz))  {
		std::stringstream ss;
		ss << "Trying to access array("<< w << "," << x << "," << y << "," << z << ")";
		ss << " while array is ("<< nw << "," << nx << "," << ny << "," << nz << ")";
		throw IndexOutOfBoundsException(ss.str(), AT);
	}
#endif
	return vecData[w + nw*x + y*nwnx + z*nwnxny];
}


template<class T> Array4D<T>::Array4D() : vecData(), nw(0), nx(0), ny(0), nz(0), nwnx(0), nwnxny(0), keep_nodata(true)
{
}

template<class T> Array4D<T>::Array4D(const Array4D<T>& i_array4D,
                                      const unsigned int& i_nw,const unsigned int& i_nx, const unsigned int& i_ny, const unsigned int& i_nz,
                                      const unsigned int& i_sizeW, const unsigned int& i_sizeX, const unsigned int& i_sizeY, const unsigned int& i_sizeZ)
                               : vecData(i_sizeW*i_sizeX*i_sizeY*i_sizeZ), nw(i_sizeW), nx(i_sizeX), ny(i_sizeY), nz(i_sizeZ), nwnx(i_sizeW*i_sizeX), nwnxny(i_sizeW*i_sizeX*i_sizeY), keep_nodata(true)
{
	subset(i_array4D, i_nw, i_nx, i_ny, i_nz, i_sizeW, i_sizeX, i_sizeY, i_sizeZ);
}

template<class T> void Array4D<T>::subset(const Array4D<T>& i_array4D,
                                     const unsigned int& i_nw, const unsigned int& i_nx, const unsigned int& i_ny, const unsigned int& i_nz,
                                     const unsigned int& i_sizeW, const unsigned int& i_sizeX, const unsigned int& i_sizeY, const unsigned int& i_sizeZ)
{

	if (((i_nw+i_sizeW) > i_array4D.nw) || ((i_nx+i_sizeX) > i_array4D.nx) || ((i_ny+i_sizeY) > i_array4D.ny) || ((i_nz+i_sizeZ) > i_array4D.nz)) {
		std::stringstream ss;
		ss << "Trying to cut an array of size ("<< nw << "," << nx << "," << ny << "," << nz << ") ";
		ss << "to size (" << i_sizeW << "," << i_sizeX << "," << i_sizeY << "," << i_sizeZ << ") ";
		ss << "starting at ("<< i_nw << "," << i_nx << "," << i_ny << "," << i_nz << ")";
		throw IndexOutOfBoundsException(ss.str(), AT);
	}

	if ((i_sizeW == 0) || (i_sizeX == 0) || (i_sizeY == 0) || (i_sizeZ == 0)) //the space has to make sense
		throw IndexOutOfBoundsException("Trying to cut an array into a null sized array!", AT);

	resize(i_sizeW, i_sizeX, i_sizeY, i_sizeZ); //create new Array4D object
	//Copy by value subspace
	for (unsigned int ii=0; ii<nz; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			for (unsigned int kk=0; kk<nx; kk++) {
				for (unsigned int ll=0; ll<nw; ll++)
				//Running through the vector in order of memory alignment
					operator()(ll,kk,jj,ii) = i_array4D(i_nw+ll, i_nx+kk, i_ny+jj, i_nz+ii);
			}
		}
	}
}

template<class T> void Array4D<T>::fill(const Array4D<T>& i_array4D, const unsigned int& i_nw, const unsigned int& i_nx, const unsigned int& i_ny, const unsigned int& i_nz)
{
	unsigned int i_sizeW, i_sizeX, i_sizeY, i_sizeZ;
	i_array4D.size(i_sizeW, i_sizeX, i_sizeY, i_sizeZ);
	fill(i_array4D, i_nw, i_nx, i_ny, i_nz, i_sizeW, i_sizeX, i_sizeY, i_sizeZ);
}

template<class T> void Array4D<T>::fill(const Array4D<T>& i_array4D,
                                     const unsigned int& i_nw, const unsigned int& i_nx, const unsigned int& i_ny, const unsigned int& i_nz,
                                     const unsigned int& i_sizeW, const unsigned int& i_sizeX, const unsigned int& i_sizeY, const unsigned int& i_sizeZ)
{

	if (((i_nw+i_sizeW) > i_array4D.nw) || ((i_nx+i_sizeX) > i_array4D.nx) || ((i_ny+i_sizeY) > i_array4D.ny) || ((i_nz+i_sizeZ) > i_array4D.nz)) {
		std::stringstream ss;
		ss << "Filling an array of size (" << nw << "," << nx << "," << ny << "," << nz << ") ";
		ss << "with an array of size (" << i_sizeW << "," << i_sizeX << "," << i_sizeY << "," << i_sizeZ << ") ";
		ss << "starting at (" << i_nw << "," << i_nx << "," << i_ny << "," << i_nz << ")";
		throw IndexOutOfBoundsException(ss.str(), AT);
	}

	if ((i_sizeW == 0) || (i_sizeX == 0) || (i_sizeY == 0) || (i_sizeZ == 0)) //the space has to make sense
		throw IndexOutOfBoundsException("Filling an array with a null sized array!", AT);

	//Copy by value subspace
	for (unsigned int ii=i_nz; ii<(i_nz+i_sizeZ); ii++) {
		for (unsigned int jj=i_ny; jj<(i_ny+i_sizeY); jj++) {
			for (unsigned int kk=i_nx; kk<(i_nx+i_sizeX); kk++) {
				for (unsigned int ll=i_nw; ll<(i_nw+i_sizeW); ll++) {
					const unsigned int iw = ll-i_nw;
					const unsigned int ix = kk-i_nx;
					const unsigned int iy = jj-i_ny;
					const unsigned int iz = ii-i_nz;
					operator()(ll,kk,jj,ii) = i_array4D(iw,ix, iy, iz);
				}
			}
		}
	}
}


template<class T> Array4D<T>::Array4D(const unsigned int& anw, const unsigned int& anx, const unsigned int& any, const unsigned int& anz)
                             : vecData(anw*anx*any*anz), nw(anw), nx(anx), ny(any), nz(anz), nwnx(anw*anx), nwnxny(anw*anx*any),  keep_nodata(true)
{
	//resize(anw, anx, any, anz);
}

template<class T> Array4D<T>::Array4D(const unsigned int& anw, const unsigned int& anx, const unsigned int& any, const unsigned int& anz, const T& init)
                             : vecData(anw*anx*any*anz, init), nw(anw), nx(anx), ny(any), nz(anz), nwnx(anw*anx), nwnxny(anw*anx*any),  keep_nodata(true)
{
	//resize(anw, anx, any, anz, init);
}

template<class T> void Array4D<T>::setKeepNodata(const bool i_keep_nodata) {
	keep_nodata = i_keep_nodata;
}

template<class T> bool Array4D<T>::getKeepNodata() {
	return keep_nodata;
}

template<class T> void Array4D<T>::resize(const unsigned int& anw, const unsigned int& anx, const unsigned int& any, const unsigned int& anz) {
	clear();  //we won't be able to "rescue" old values, so we reset the whole vector
	vecData.resize(anw*anx*any*anz);
	nw = anw;
	nx = anx;
	ny = any;
	nz = anz;
	nwnx = nw*nx;
	nwnxny = nw*nx*ny;
}

template<class T> void Array4D<T>::resize(const unsigned int& anw, const unsigned int& anx, const unsigned int& any, const unsigned int& anz, const T& init) {
	clear();  //we won't be able to "rescue" old values, so we reset the whole vector
	vecData.resize(anw*anx*any*anz, init);
	nw = anw;
	nx = anx;
	ny = any;
	nz = anz;
	nwnx = nw*nx;
	nwnxny = nw*nx*ny;
}

template<class T> void Array4D<T>::size(unsigned int& anw, unsigned int& anx, unsigned int& any, unsigned int& anz) const {
	anw=nw;
	anx=nx;
	any=ny;
	anz=nz;
}

template<class T> unsigned int Array4D<T>::getNw() const {
	return nw;
}

template<class T> unsigned int Array4D<T>::getNx() const {
	return nx;
}

template<class T> unsigned int Array4D<T>::getNy() const {
	return ny;
}

template<class T> unsigned int Array4D<T>::getNz() const {
	return nz;
}

template<class T> void Array4D<T>::clear() {
	vecData.clear();
	nw = nx = ny = nz = nwnx = nwnxny = 0;
}

template<class T> bool Array4D<T>::isEmpty() const {
	return (nw==0 && nx==0 && ny==0 && nz==0);
}

template<class T> const std::string Array4D<T>::toString() const {
	std::stringstream os;
	os << "<array4d>\n";
	for (unsigned int ll=0; ll<nw; ll++) {
		os << "dim4[" << ll << "]\n";
		for (unsigned int kk=0; kk<nz; kk++) {
			os << "depth[" << kk << "]\n";
			for(unsigned int ii=0; ii<nx; ii++) {
				for (unsigned int jj=0; jj<ny; jj++) {
					os << operator()(ii,jj,kk,ll) << " ";
				}
				os << "\n";
			}
		}
	}
	os << "</array4d>\n";
	return os.str();
}

template<class P> std::iostream& operator<<(std::iostream& os, const Array4D<P>& array) {
	os.write(reinterpret_cast<const char*>(&array.keep_nodata), sizeof(array.keep_nodata));
	os.write(reinterpret_cast<const char*>(&array.nx), sizeof(array.nx));
	os.write(reinterpret_cast<const char*>(&array.ny), sizeof(array.ny));
	os.write(reinterpret_cast<const char*>(&array.nz), sizeof(array.nz));
	os.write(reinterpret_cast<const char*>(&array.nw), sizeof(array.nw));
	os.write(reinterpret_cast<const char*>(&array.vecData[0]), array.nx*array.ny*array.nz*array.nw*sizeof(P));
	return os;
}

template<class P> std::iostream& operator>>(std::iostream& is, Array4D<P>& array) {
	is.read(reinterpret_cast<char*>(&array.keep_nodata), sizeof(array.keep_nodata));
	is.read(reinterpret_cast<char*>(&array.nx), sizeof(array.nx));
	is.read(reinterpret_cast<char*>(&array.ny), sizeof(array.ny));
	is.read(reinterpret_cast<char*>(&array.nz), sizeof(array.nz));
	is.read(reinterpret_cast<char*>(&array.nw), sizeof(array.nw));
	array.vecData.resize(array.nx*array.ny*array.nz*array.nw);
	is.read(reinterpret_cast<char*>(&array.vecData[0]), array.nx*array.ny*array.nz*array.nw*sizeof(P)); //30 times faster than assign() or copy()
	return is;
}

template<class T> T Array4D<T>::getMin() const {

	T min = std::numeric_limits<T>::max();
	const unsigned int nwyz = nwnxny*nz;

	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nwyz; jj++) {
			const T val = operator()(jj);
			if(val<min) min=val;
		}
		return min;
	} else {
		for (unsigned int jj=0; jj<nwyz; jj++) {
			const T val = operator()(jj);
			if(val!=IOUtils::nodata && val<min) min=val;
		}
		if(min!=std::numeric_limits<T>::max()) return min;
		else return (T)IOUtils::nodata;
	}
}

template<class T> T Array4D<T>::getMax() const {

	T max = -std::numeric_limits<T>::max();
	const unsigned int nwyz = nwnxny*nz;

	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nwyz; jj++) {
			const T val = operator()(jj);
			if(val>max) max=val;
		}
		return max;
	} else {
		for (unsigned int jj=0; jj<nwyz; jj++) {
			const T val = operator()(jj);
			if(val!=IOUtils::nodata && val>max) max=val;
		}
		if(max!=-std::numeric_limits<T>::max()) return max;
		else return (T)IOUtils::nodata;
	}
}

template<class T> T Array4D<T>::getMean() const {

	T mean = 0;
	const unsigned int nwyz = nwnxny*nz;

	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nwyz; jj++) {
			const T val = operator()(jj);
			mean += val;
		}
		if(nwyz>0) return mean/(T)(nwyz);
		else return (T)0;
	} else {
		unsigned int count = 0;
		for (unsigned int jj=0; jj<nwyz; jj++) {
			const T val = operator()(jj);
			if(val!=IOUtils::nodata) {
				mean += val;
				count++;
			}
		}
		if(count>0) return mean/(T)(count);
		else return (T)IOUtils::nodata;
	}
}

template<class T> size_t Array4D<T>::getCount() const
{
	const unsigned int nwyz = nwnxny*nz;

	if(keep_nodata==false) {
		return (size_t)nwyz;
	} else {
		size_t count = 0;
		for (unsigned int ii=0; ii<nwyz; ii++) {
			if(vecData[ii]!=IOUtils::nodata) count++;
		}
		return count;
	}
}

template<class T> void Array4D<T>::abs() {
	if(std::numeric_limits<T>::is_signed) {
		const unsigned int nwyz = nwnxny*nz;
		if(keep_nodata==false) {
			for (unsigned int ii=0; ii<nwyz; ii++) {
				T& val = operator()(ii);
				if(val<0) val=-val;
			}
		} else {
			for (unsigned int ii=0; ii<nwyz; ii++) {
				T& val = operator()(ii);
				if(val<0 && val!=IOUtils::nodata) val=-val;
			}
		}
	}
}


template<class T> const Array4D<T> Array4D<T>::getAbs() const {
	Array4D<T> result = *this; //make a copy
	result.abs(); //already implemented

	return result;
}

//arithmetic operators
template<class T> bool Array4D<T>::checkEpsilonEquality(const Array4D<double>& rhs, const double& epsilon) const {
	if(nw!=rhs.nw || nx!=rhs.nx || ny!=rhs.ny || nz!=rhs.nz) return false;

	const unsigned int nwyz = nwnxny*nz;
	for (unsigned int jj=0; jj<nwyz; jj++)
		if(IOUtils::checkEpsilonEquality(vecData[jj], rhs.vecData[jj], epsilon)==false) return false;

	return true;
}

template<class T> bool Array4D<T>::checkEpsilonEquality(const Array4D<double>& rhs1, const Array4D<double>& rhs2, const double& epsilon) {
	return rhs1.checkEpsilonEquality(rhs2, epsilon);
}

template<class T> Array4D<T>& Array4D<T>::operator=(const Array4D<T>& source) {
	if(this != &source) {
		keep_nodata = source.keep_nodata;
		nw = source.nw;
		nx = source.nx;
		ny = source.ny;
		nz = source.nz;
		nwnx = source.nwnx;
		nwnxny = source.nwnxny;
		vecData = source.vecData;
	}
	return *this;
}

template<class T> Array4D<T>& Array4D<T>::operator=(const T& value) {
	std::fill(vecData.begin(), vecData.end(), value);
	return *this;
}

template<class T> Array4D<T>& Array4D<T>::operator+=(const Array4D<T>& rhs)
{
	//They have to have equal size
	if ((rhs.nw != nw) || (rhs.nx != nx) || (rhs.ny != ny) || (rhs.nz != nz)) {
		std::stringstream ss;
		ss << "Trying to add two Array4D objects with different dimensions: ";
		ss << "(" << nw << "," << nx << "," << ny << "," << nz << ") + (" << rhs.nw << "," << rhs.nx << "," << rhs.ny << "," << rhs.nz << ")";
		throw IOException(ss.str(), AT);
	}
	//Add to every single member of the Array4D<T>
	const unsigned int nwyz = nwnxny*nz;
	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nwyz; jj++) {
			operator()(jj) += rhs(jj);
		}
	} else {
		for (unsigned int jj=0; jj<nwyz; jj++) {
			if(operator()(jj)==IOUtils::nodata || rhs(jj)==IOUtils::nodata)
				operator()(jj) = IOUtils::nodata;
			else
				operator()(jj) += rhs(jj);
		}
	}

	return *this;
}

template<class T> const Array4D<T> Array4D<T>::operator+(const Array4D<T>& rhs)
{
	Array4D<T> result = *this; //make a copy
	result += rhs; //already implemented

	return result;
}

template<class T> Array4D<T>& Array4D<T>::operator+=(const T& rhs)
{
	//Add to every single member of the Array4D<T>
	const unsigned int nwyz = nwnxny*nz;
	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nwyz; jj++) {
			operator()(jj) += rhs;
		}
	} else {
		for (unsigned int jj=0; jj<nwyz; jj++) {
			if(operator()(jj)!=IOUtils::nodata)
				operator()(jj) += rhs;
		}
	}

	return *this;
}

template<class T> const Array4D<T> Array4D<T>::operator+(const T& rhs)
{
	Array4D<T> result = *this;
	result += rhs; //already implemented

	return result;
}

template<class T> Array4D<T>& Array4D<T>::operator-=(const Array4D<T>& rhs)
{
	//They have to have equal size
	if ((rhs.nw != nw) || (rhs.nx != nx) || (rhs.ny != ny) || (rhs.nz != nz)) {
		std::stringstream ss;
		ss << "Trying to substract two Array4D objects with different dimensions: ";
		ss << "(" << nw << "," << nx << "," << ny << "," << nz << ") - (" << rhs.nw << "," << rhs.nx << "," << rhs.ny << "," << rhs.nz << ")";
		throw IOException(ss.str(), AT);
	}
	//Substract to every single member of the Array4D<T>
	const unsigned int nwyz = nwnxny*nz;
	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nwyz; jj++) {
			operator()(jj) -= rhs(jj);
		}
	} else {
		for (unsigned int jj=0; jj<nwyz; jj++) {
			if(operator()(jj)==IOUtils::nodata || rhs(jj)==IOUtils::nodata)
				operator()(jj) = IOUtils::nodata;
			else
				operator()(jj) -= rhs(jj);
		}
	}

	return *this;
}

template<class T> const Array4D<T> Array4D<T>::operator-(const Array4D<T>& rhs)
{
	Array4D<T> result = *this; //make a copy
	result -= rhs; //already implemented

	return result;
}

template<class T> Array4D<T>& Array4D<T>::operator-=(const T& rhs)
{
	//Substract to every single member of the Array4D<T>
	const unsigned int nwxyz = nwnxny*nz;
	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nwxyz; jj++) {
			operator()(jj) -= rhs;
		}
	} else {
		for (unsigned int jj=0; jj<nwxyz; jj++) {
			if(operator()(jj)!=IOUtils::nodata)
				operator()(jj) -= rhs;
		}
	}

	return *this;
}

template<class T> const Array4D<T> Array4D<T>::operator-(const T& rhs)
{
	Array4D<T> result = *this;
	result -= rhs; //already implemented

	return result;
}

template<class T> Array4D<T>& Array4D<T>::operator*=(const Array4D<T>& rhs)
{
	//They have to have equal size
	if ((rhs.nw != nw) || (rhs.nx != nx) || (rhs.ny != ny) || (rhs.nz != nz)) {
		std::stringstream ss;
		ss << "Trying to multiply two Array4D objects with different dimensions: ";
		ss << "("<< nw << "," << nx << "," << ny << "," << nz << ") * (" << rhs.nw << "," << rhs.nx << "," << rhs.ny << "," << rhs.nz << ")";
		throw IOException(ss.str(), AT);
	}
	//Multiply every single member of the Array4D<T>
	const unsigned int nwxyz = nwnxny*nz;
	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nwxyz; jj++) {
			operator()(jj) *= rhs(jj);
		}
	} else {
		for (unsigned int jj=0; jj<nwxyz; jj++) {
			if(operator()(jj)==IOUtils::nodata || rhs(jj)==IOUtils::nodata)
				operator()(jj) = IOUtils::nodata;
			else
				operator()(jj) *= rhs(jj);
		}
	}

	return *this;
}

template<class T> const Array4D<T> Array4D<T>::operator*(const Array4D<T>& rhs)
{
	Array4D<T> result = *this; //make a copy
	result *= rhs; //already implemented

	return result;
}

template<class T> Array4D<T>& Array4D<T>::operator*=(const T& rhs)
{
	//Multiply every single member of the Array4D<T>
	const unsigned int nwxyz = nwnxny*nz;
	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nwxyz; jj++) {
			operator()(jj) *= rhs;
		}
	} else {
		for (unsigned int jj=0; jj<nwxyz; jj++) {
			if(operator()(jj)!=IOUtils::nodata)
				operator()(jj) *= rhs;
		}
	}

	return *this;
}

template<class T> const Array4D<T> Array4D<T>::operator*(const T& rhs)
{
	Array4D<T> result = *this;
	result *= rhs; //already implemented

	return result;
}

template<class T> Array4D<T>& Array4D<T>::operator/=(const Array4D<T>& rhs)
{
	//They have to have equal size
	if ((rhs.nw != nw) || (rhs.nx != nx) || (rhs.ny != ny) || (rhs.nz != nz)) {
		std::stringstream ss;
		ss << "Trying to divide two Array4D objects with different dimensions: ";
		ss << "(" << nw << "," << nx << "," << ny << "," << nz << ") / (" << rhs.nw << "," << rhs.nx << "," << rhs.ny << "," << rhs.nz << ")";
		throw IOException(ss.str(), AT);
	}
	//Divide every single member of the Array4D<T>
	const unsigned int nwxyz = nwnxny*nz;
	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nwxyz; jj++) {
			operator()(jj) /= rhs(jj);
		}
	} else {
		for (unsigned int jj=0; jj<nwxyz; jj++) {
			if(operator()(jj)==IOUtils::nodata || rhs(jj)==IOUtils::nodata)
				operator()(jj) = IOUtils::nodata;
			else
				operator()(jj) /= rhs(jj);
		}
	}

	return *this;
}

template<class T> const Array4D<T> Array4D<T>::operator/(const Array4D<T>& rhs)
{
	Array4D<T> result = *this; //make a copy
	result /= rhs; //already implemented

	return result;
}

template<class T> Array4D<T>& Array4D<T>::operator/=(const T& rhs)
{
	//Divide every single member of the Array4D<T>
	const unsigned int nwxyz = nwnxny*nz;
	if(keep_nodata==false) {
		for (unsigned int jj=0; jj<nwxyz; jj++) {
			operator()(jj) /= rhs;
		}
	} else {
		for (unsigned int jj=0; jj<nwxyz; jj++) {
			if(operator()(jj)!=IOUtils::nodata)
				operator()(jj) /= rhs;
		}
	}

	return *this;
}

template<class T> const Array4D<T> Array4D<T>::operator/(const T& rhs)
{
	Array4D<T> result = *this;
	result /= rhs; //already implemented

	return result;
}

} //end namespace mio

#endif
