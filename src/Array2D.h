#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <vector>
#include "IOExceptions.h"
#include "IOUtils.h"
#include <limits>

#define NOSAFECHECKS

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
		Array2DProxy(Array2D<T>& _array2D, const unsigned int& _anx) : array2D(_array2D), anx(_anx){}
		Array2D<T>& array2D;
		const unsigned int anx;
}; 

/**
 * @class Array2D
 * @brief The template class Array2D is a 2D Array (Matrix) able to hold any type of object as datatype
 *
 * @author Thomas Egger
 */
template<class T> class Array2D {
	public:
		Array2D();
		Array2D(const unsigned int& anx, const unsigned int& any);
		Array2D(const unsigned int& anx, const unsigned int& any, const T& init);

		/**
		* A constructor that can be used to create an Array2D object that is contained in the
		* one passed as _array2D argument. The resulting Array2D object is a by value copy of
		* a subplane of the plane spanned by the _array2D
		*/
		Array2D(const Array2D<T>& _array2D, const unsigned int& _nx, const unsigned int& _ny, 
			    const unsigned int& _ncols, const unsigned int& _nrows);

		void resize(const unsigned int& nx, const unsigned int& ny);
		void resize(const unsigned int& nx, const unsigned int& ny, const T& init);
		void size(unsigned int& nx, unsigned int& ny) const;
		T getMin();
		T getMax();

		void clear();
		T& operator ()(const unsigned int& x, const unsigned int& y);
		const T operator ()(const unsigned int& x, const unsigned int& y) const;
		Array2DProxy<T> operator[](const unsigned int& i);

	protected:
		std::vector<T> vecData;
		unsigned int nx;
		unsigned int ny;
};


template<class T> T& Array2D<T>::operator()(const unsigned int& x, const unsigned int& y) {
#ifndef NOSAFECHECKS
	if ((x >= nx) || (y >= ny)) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif

	//the 2D array is stored by columns in a 1D vector. Each column follows the previous one.
	//This matches the data usage in Alpine3D, data access being done by columns.
	return vecData[x*ny + y];
}

template<class T> const T Array2D<T>::operator()(const unsigned int& x, const unsigned int& y) const {
#ifndef NOSAFECHECKS
	if ((x >= nx) || (y >= ny)) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif
	return vecData[x*ny + y];
}

template<class T> Array2DProxy<T> Array2D<T>::operator[](const unsigned int& i) {
	return Array2DProxy<T>(*this, i); 
}


template<class T> Array2D<T>::Array2D() {
	nx = ny = 0;
}

template<class T> Array2D<T>::Array2D(const Array2D<T>& _array2D, const unsigned int& _nx, const unsigned int& _ny, 
			    const unsigned int& _ncols, const unsigned int& _nrows)
{
	if (((_nx+_ncols) > _array2D.nx) || ((_ny+_nrows) > _array2D.ny))
		throw IndexOutOfBoundsException("", AT);

	if ((_ncols == 0) || (_nrows == 0)) //the plane to copy has to make sense
		throw IndexOutOfBoundsException("", AT);

	resize(_ncols, _nrows); //create new Array2D object

	//Copy by value subspace
	for (unsigned int ii=0; ii<ny; ii++) { 
		for (unsigned int jj=0; jj<nx; jj++) {
			//Running through the vector in order of memory alignment HACK
			operator()(jj,ii) = _array2D(_nx+jj, _ny+ii);
		}
	}

}

template<class T> Array2D<T>::Array2D(const unsigned int& anx, const unsigned int& any, const T& init) {
	nx = ny = 0;
	resize(anx,any,init);
}

template<class T> Array2D<T>::Array2D(const unsigned int& anx, const unsigned int& any) {
	nx = ny = 0;
	resize(anx,any);
}

template<class T> void Array2D<T>::resize(const unsigned int& anx, const unsigned int& any) {
	clear();

	if ((anx > 0) && (any > 0)) {
		vecData.resize(anx*any);
		nx = anx;
		ny = any;
	} else {
		throw IndexOutOfBoundsException("", AT);    
	}
}

template<class T> void Array2D<T>::resize(const unsigned int& anx, const unsigned int& any, const T& init) {
	resize(anx, any);

	for (unsigned int ii=0; ii<nx; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			operator()(ii,jj) = init;
		}
	}
}

template<class T> void Array2D<T>::size(unsigned int& anx, unsigned int& any) const{
	anx=nx;
	any=ny;
}

template<class T> void Array2D<T>::clear() {
	vecData.clear();
	nx=ny=0;
}

template<class T> T Array2D<T>::getMin() {

	T min = std::numeric_limits<T>::max();

	for (unsigned int ii=0; ii<nx; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			const T val = vecData[ii*ny + jj];
			if(val<min) min=val;
		}
	}
	
	return min;
}

template<class T> T Array2D<T>::getMax() {

	T max = -std::numeric_limits<T>::max();

	for (unsigned int ii=0; ii<nx; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			const T val = vecData[ii*ny + jj];
			if(val>max) max=val;
		}
	}

	return max;
}

#endif
