#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <vector>
#include "IOExceptions.h"

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
	T& operator[](unsigned int j) { return array2D(anx, j); }
	~CArray2DProxy(){}

 private:
 CArray2DProxy(CArray2D<T>& _array2D, int _anx) : array2D(_array2D), anx(_anx){}
	CArray2D<T>& array2D;
	unsigned int anx;
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
	CArray2D(int anx, int any);
	~CArray2D();

	void Create(int nx, int ny);
	void Create(int nx, int ny, T init);
	void GetSize(int &nx, int &ny);

	void Destroy();
	CArray2D<T>& operator=(CArray2D &val);
	T& operator ()(unsigned int x, unsigned int y);
	const T operator ()(unsigned int x, unsigned int y) const;
	CArray2DProxy<T> operator[](unsigned int i);
	
 protected:
	std::vector<T> vecData;
	unsigned int nx;
	unsigned int ny;
};


template<class T> T& CArray2D<T>::operator()(unsigned int x, unsigned int y){
#ifndef NOSAFECHECKS
	if ((x >= nx) || (y >= ny))
		THROW IndexOutOfBoundsException("", AT);
#endif

	//the 2D array is stored by columns in a 1D vector. Each column follows the previous one.
	//This matches the data usage in Alpine3D, data access being done by columns.
	return vecData[x*ny + y];
}

template<class T> const T CArray2D<T>::operator()(unsigned int x, unsigned int y) const {
#ifndef NOSAFECHECKS
	if ((x >= nx) || (y >= ny))
		THROW IndexOutOfBoundsException("", AT);
#endif
	return vecData[x*ny + y];
}

template<class T> CArray2DProxy<T> CArray2D<T>::operator[](unsigned int i){
	return CArray2DProxy<T>(*this, i); 
}


template<class T> CArray2D<T>::CArray2D(){
	nx = ny = 0;
}

template<class T> CArray2D<T>::CArray2D(int anx, int any){
	nx = ny = 0;
	Create(anx,any);
}

template<class T> CArray2D<T>::~CArray2D(){
	Destroy();
}

template<class T> void CArray2D<T>::Create(int anx, int any){
	Destroy();

	if ((anx > 0) && (any > 0)){
		vecData.resize(anx*any);
		nx = anx;
		ny = any;
	} else {
		THROW IndexOutOfBoundsException("", AT);    
	}
}

template<class T> void CArray2D<T>::Create(int anx, int any, T init){
	Create(anx, any);

	for (unsigned int ii=0; ii<nx; ii++){
		for (unsigned int jj=0; jj<ny; jj++){
			operator()(ii,jj) = init;
		}
	}
}

template<class T> void CArray2D<T>::GetSize(int &anx, int &any){
	anx=nx;
	any=ny;
}

template<class T> void CArray2D<T>::Destroy(){
	vecData.clear();
	nx=ny=0;
}

template<class T> CArray2D<T>& CArray2D<T>::operator=(CArray2D<T>& val){
	int anx,any;
	val.GetSize(anx,any);
	
	vecData = val.vecData;
	nx = anx;
	ny = any;

	return *this;
}

#endif
