#ifndef ARRAY3D_H
#define ARRAY3D_H

#include "IOExceptions.h"
#include "Array2D.h"
#include <vector>

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
		Array3DProxy2<T> operator[](unsigned int _any) {
			return Array3DProxy2<T>(array3D, anx, _any); 
		}

 	private:
 		Array3DProxy(Array3D<T>& _array3D, unsigned int& _anx) : array3D(_array3D), anx(_anx){}
		Array3D<T>& array3D;
		unsigned int anx;
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
		T& operator[](unsigned int _anz) {
			return array3D(anx, any, _anz);
		}

	private:
 		Array3DProxy2(Array3D<T>& _array3D, unsigned int& _anx, unsigned int& _any) : array3D(_array3D), anx(_anx), any(_any){}
		Array3D<T>& array3D;
		unsigned int anx;
		unsigned int any;
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
		Array3D(const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz);
		Array3D(const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz, const T& _init);

		void resize(const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz);
		void resize(const unsigned int& _nx, const unsigned int& _ny, const unsigned int& _nz, const T& _init);
		void size(unsigned int& _nx, unsigned int& _ny, unsigned int& _nz) const;
		void clear();

		T& operator ()(const unsigned int& x, const unsigned int& y, const unsigned int& z);
		const T operator ()(const unsigned int& x, const unsigned int& y, const unsigned int& z) const;
		Array3DProxy<T> operator[](unsigned int i);

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
	return vecData[x*ny + y + z*nxny];
}

template<class T> const T Array3D<T>::operator()(const unsigned int& x, const unsigned int& y, const unsigned int& z) const {
#ifndef NOSAFECHECKS
	if ((x >= nx) || (y >= ny) || (z >= nz)) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif
	return vecData[x*ny + y + z*nxny];
}

template<class T> Array3DProxy<T> Array3D<T>::operator[](unsigned int i) {
	return Array3DProxy<T>(*this, i); 
}


template<class T> Array3D<T>::Array3D() {
	nx = ny = nz = nxny = 0;
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

	for (unsigned int ii=0; ii<nx; ii++) {
		for (unsigned int jj=0; jj<ny; jj++) {
			for (unsigned int kk=0; kk<nz; kk++) {
				operator()(ii,jj,kk) = _init;
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

#endif
