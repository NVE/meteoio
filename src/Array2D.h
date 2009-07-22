#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <vector>
#include "IOExceptions.h"

#define NOSAFECHECKS

template <class T> class CArray2D;

/**
 * @class CArray2DProxy
 * @brief The template class CArray2DProxy is a helper class for the template class CArray2D
 *        with the purpose of adding the [][] operator to CArray2D
 *
 * @author Thomas Egger
 */
template <class T> class CArray2DProxy {
	public:
		friend class CArray2D<T>;
		T& operator[](const unsigned int& j) {
			return array2D(anx, j);
		}

	private:
		CArray2DProxy(CArray2D<T>& _array2D, const unsigned int& _anx) : array2D(_array2D), anx(_anx){}
		CArray2D<T>& array2D;
		const unsigned int anx;
}; 

/**
 * @class CArray2D
 * @brief The template class CArray2D is a 2D Array (Matrix) able to hold any type of object as datatype
 *
 * @author Thomas Egger
 */
template<class T> class CArray2D {
	public:
		CArray2D();
		CArray2D(const unsigned int& anx, const unsigned int& any);
		CArray2D(const unsigned int& anx, const unsigned int& any, const T& init);

		/**
		* A constructor that can be used to create an CArray2D object that is contained in the
		* one passed as _array2D argument. The resulting CArray2D object is a by value copy of
		* a subplane of the plane spanned by the _array2D
		*/
		CArray2D(const CArray2D<T>& _array2D, const unsigned int& _nx, const unsigned int& _ny, 
			    const unsigned int& _ncols, const unsigned int& _nrows);

		void Create(const unsigned int& nx, const unsigned int& ny);
		void Create(const unsigned int& nx, const unsigned int& ny, const T& init);
		void GetSize(unsigned int& nx, unsigned int& ny) const;

		void Destroy();
		T& operator ()(const unsigned int& x, const unsigned int& y);
		const T operator ()(const unsigned int& x, const unsigned int& y) const;
		CArray2DProxy<T> operator[](const unsigned int& i);

	protected:
		std::vector<T> vecData;
		unsigned int nx;
		unsigned int ny;
};


template<class T> T& CArray2D<T>::operator()(const unsigned int& x, const unsigned int& y) {
#ifndef NOSAFECHECKS
	if ((x >= nx) || (y >= ny)) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif

	//the 2D array is stored by columns in a 1D vector. Each column follows the previous one.
	//This matches the data usage in Alpine3D, data access being done by columns.
	return vecData[x*ny + y];
}

template<class T> const T CArray2D<T>::operator()(const unsigned int& x, const unsigned int& y) const {
#ifndef NOSAFECHECKS
	if ((x >= nx) || (y >= ny)) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif
	return vecData[x*ny + y];
}

template<class T> CArray2DProxy<T> CArray2D<T>::operator[](const unsigned int& i) {
	return CArray2DProxy<T>(*this, i); 
}


template<class T> CArray2D<T>::CArray2D() {
	nx = ny = 0;
}

template<class T> CArray2D<T>::CArray2D(const CArray2D<T>& _array2D, const unsigned int& _nx, const unsigned int& _ny, 
			    const unsigned int& _ncols, const unsigned int& _nrows)
{
	if (((_nx+_ncols) > _array2D.nx) || ((_ny+_nrows) > _array2D.ny))
		throw IndexOutOfBoundsException("", AT);

	if ((_ncols == 0) || (_nrows == 0)) //the plane to copy has to make sense
		throw IndexOutOfBoundsException("", AT);

	Create(_ncols, _nrows); //create new Array3D object

	//Copy by value subspace
	for (unsigned int ii=0; ii<ny; ii++) { 
		for (unsigned int jj=0; jj<nx; jj++) {
			//Running through the vector in order of memory alignment
			operator()(jj,ii) = _array2D(_nx+jj, _ny+ii);
		}
	}

}

template<class T> CArray2D<T>::CArray2D(const unsigned int& anx, const unsigned int& any, const T& init) {
	nx = ny = 0;
	Create(anx,any,init);
}

template<class T> CArray2D<T>::CArray2D(const unsigned int& anx, const unsigned int& any) {
	nx = ny = 0;
	Create(anx,any);
}

template<class T> void CArray2D<T>::Create(const unsigned int& anx, const unsigned int& any) {
	Destroy();

	if ((anx > 0) && (any > 0)) {
		vecData.resize(anx*any);
		nx = anx;
		ny = any;
	} else {
		throw IndexOutOfBoundsException("", AT);    
	}
}

template<class T> void CArray2D<T>::Create(const unsigned int& anx, const unsigned int& any, const T& init) {
	Create(anx, any);

	for (unsigned int ii=0; ii<nx; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			operator()(ii,jj) = init;
		}
	}
}

template<class T> void CArray2D<T>::GetSize(unsigned int& anx, unsigned int& any) const{
	anx=nx;
	any=ny;
}

template<class T> void CArray2D<T>::Destroy() {
	vecData.clear();
	nx=ny=0;
}

#endif
