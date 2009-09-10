#ifndef ARRAY3D_H
#define ARRAY3D_H

#include <vector>
#include <limits>
#include "IOExceptions.h"

#define NOSAFECHECKS

template <class T> class Array3D;
template <class T> class Array3DProxy2;

/**
 * @class Array3DProxy
 * @brief The template class Array3DProxy is a helper class for the template class Array3D
 *        with the purpose of adding the [][][] operator to Array3D
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
 * @brief The template class Array3D is a 3D Array (Tensor) able to hold any type of object as datatype
 *
 * @date  2009-07-19
 * @author Thomas Egger
 */
template<class T> class Array3D {
	public:
		Array3D();

		/**
		* A constructor that can be used to create an Array3D object that is contained in the
		* one passed as _array3D argument. The resulting Array3D object is a by value copy of
		* a subspace of the space spanned by the _array3D
		*/
		Array3D(const Array3D<T>& _array3D,
			   const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz,
			   const unsigned int& _cols, const unsigned int& _nrows, const unsigned int& _ndepth);
		Array3D(const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz);
		Array3D(const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz, const T& _init);

		void resize(const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz);
		void resize(const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz, const T& _init);
		void size(unsigned int& _nx, unsigned int& _ny, unsigned int& _nz) const;
		void clear();
		T getMin();
		T getMax();

		T& operator ()(const unsigned int& x, const unsigned int& y, const unsigned int& z);
		const T operator ()(const unsigned int& x, const unsigned int& y, const unsigned int& z) const;
		Array3DProxy<T> operator[](const unsigned int& i);

	protected:
		std::vector<T> vecData; ///< The actual objects are stored in a one-dimensional vector
		unsigned int nx;
		unsigned int ny;
		unsigned int nz;
		unsigned int nxny; //nx times ny
};


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

template<class T> Array3D<T>::Array3D(const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz) {
	resize(_nx, _ny, _nz);
}

template<class T> Array3D<T>::Array3D(const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz, const T& _init) {
	resize(_nx, _ny, _nz, _init);
}

template<class T> void Array3D<T>::resize(const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz) {
	clear();

	if ((_nx > 0) && (_ny > 0) && (_nz > 0)) {
		vecData.resize(_nx*_ny*_nz);
		nx = _nx;
		ny = _ny;
		nz = _nz;
		nxny = nx*ny;
	} else {
		throw IndexOutOfBoundsException("", AT);    
	}
}

template<class T> void Array3D<T>::resize(const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz, const T& _init) {
	resize(_nx, _ny, _nz);

	for (unsigned int ii=0; ii<nz; ii++) { 
		for (unsigned int jj=0; jj<ny; jj++) {
			for (unsigned int kk=0; kk<nx; kk++) {
				operator()(kk,jj,ii) = _init; //Running through the vector in order of memory alignment
			}
		}
	}
}

template<class T> void Array3D<T>::size(unsigned int& _nx, unsigned int& _ny, unsigned int& _nz) const {
	_nx=nx;
	_ny=ny;
	_nz=nz;
}

template<class T> void Array3D<T>::clear() {
	vecData.clear();
	nx = ny = nz = nxny = 0;
}

template<class T> T Array3D<T>::getMin() {

	T min = std::numeric_limits<T>::max();

	for (unsigned int ii=0; ii<nx; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			for (unsigned int kk=0; kk<nx; kk++) {
				const T val = vecData[ii + jj*nx + kk*nxny];
				if(val<min) min=val;
			}
		}
	}
	
	return min;
}

template<class T> T Array3D<T>::getMax() {

	T max = -std::numeric_limits<T>::max();

	for (unsigned int ii=0; ii<nx; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			for (unsigned int kk=0; kk<nx; kk++) {
				const T val = vecData[ii + jj*nx + kk*nxny];
				if(val>max) max=val;
			}
		}
	}

	return max;
}

#endif
