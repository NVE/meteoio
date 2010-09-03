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

#include <meteoio/Matrix.h>

namespace mio {

Matrix::Matrix() {
	nx = ny = 0;
}

Matrix::Matrix(const int& anx, const int& any) {
	if(anx<0 || any<0) {
		std::stringstream tmp;
		tmp << "Trying construct a matrix with negative dimensions: ";
		tmp << "(" << anx << "," << any << ")";
		throw IOException(tmp.str(), AT);
	}
	nx = ny = 0;
	resize((unsigned)anx,(unsigned)any);
}

Matrix::Matrix(const unsigned int& anx, const unsigned int& any, const double& init) {
	nx = ny = 0;
	resize(anx,any,init);
}

Matrix::Matrix(const unsigned int& n, const double& init) {
	nx = ny = 0;
	resize(n,n,0.);
	for(unsigned int jj=0; jj<ny; jj++) operator()(jj,jj) = init;
}

std::ostream& operator<<(std::ostream& os, const Matrix& data) {
	//os << "<matrix>\n";
	const unsigned int wd=5;
	os << "\n┌ ";
	for(unsigned int jj=0; jj<(data.ny*(wd+1)); jj++)
		os << " ";
	os << " ┐\n";
	for(unsigned int jj=0; jj<data.ny; jj++) {
		os << "│ ";
		for (unsigned int ii=0; ii<data.nx; ii++) {
			os << std::setw(wd) << std::fixed << std::setprecision(2) << data(ii,jj) << " ";
		}
		os << " │\n";
	}
	os << "└ ";
	for(unsigned int jj=0; jj<(data.ny*(wd+1)); jj++)
		os << " ";
	os << " ┘\n";
	//os << "\n";
	return os;
}

/*const double Matrix::operator()(const unsigned int& x, const unsigned int& y) const {
#ifndef NOSAFECHECKS
	if ((x >= nx) || (y >= ny)) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif
	return vecData[ (x-1) + (y-1)*nx];
}*/

Matrix& Matrix::operator*=(const Matrix& rhs) {

	//check dimensions compatibility
	if(nx!=rhs.ny) {
		std::stringstream tmp;
		tmp << "Trying to multiply two matrix with incompatible dimensions: ";
		tmp << "(" << nx << "," << ny << ") * ";
		tmp << "(" << rhs.nx << "," << rhs.ny << ")";
		throw IOException(tmp.str(), AT);
	}

	//create new matrix
	Matrix result(rhs.nx,ny);

	//fill product matrix
	for(unsigned int jj=0; jj<result.ny; jj++) {
		for(unsigned int ii=0; ii<result.nx; ii++) {
			double sum=0.;
			for(unsigned int idx=0; idx<nx; idx++)
				sum+= operator()(idx,jj) * rhs(ii,idx);
			result(ii,jj) = sum;
		}
	}
	*this = result;
	return *this;
}

const Matrix Matrix::operator*(const Matrix& rhs) {
	Matrix result = *this;
	result *= rhs; //already implemented

	return result;
}

/*void Matrix::T() {

}*/

const Matrix Matrix::T() {
//other possibility: create a "transpose" flag that simply swaps the data reading...
	Matrix result(ny, nx);

	for(unsigned int jj=0; jj<result.ny; jj++) {
		for(unsigned int ii=0; ii<result.nx; ii++) {
			result(ii,jj) = operator()(jj,ii);
		}
	}
	return result;
}

double Matrix::det() const {
	if(nx!=ny) {
		std::stringstream tmp;
		tmp << "Trying to calculate the determinant of a non-square matrix ";
		tmp << "(" << nx << "," << ny << ") !";
		throw IOException(tmp.str(), AT);
	}
	Matrix U;
	LU(U);

	double product=1.;
	for(unsigned int i=0; i<nx; i++) product *= U(i,i);
	
	return product;
}

const Matrix Matrix::LU(Matrix& U) const {
//Dolittle algorithm, cf http://math.fullerton.edu/mathews/numerical/linear/dol/dol.html
//HACK: there is no permutation matrix, so it might not be able to give a decomposition...
//This is implemented for matrix storage from 0 to n-1
	if(nx!=ny) {
		std::stringstream tmp;
		tmp << "Trying to calculate the LU decomposition of a non-square matrix ";
		tmp << "(" << nx << "," << ny << ") !";
		throw IOException(tmp.str(), AT);
	}

	U = *this;
	const Matrix& A = *this;
	Matrix L(nx, 1.); //initialized as diagonal matrix, then populated

	const unsigned int n = nx;
	for(unsigned int k=1; k<=n; k++) {
		//compute U elements
		for(unsigned int j=1; j<k; j++) {
			U(j-1,k-1) = 0.;
		}
		for(unsigned int j=k; j<=n; j++) {
			double sum=0.;
			for(unsigned int m=1; m<=(k-1); m++) sum +=L(m-1,k-1)*U(j-1,m-1);
			U(j-1,k-1) = A(j-1,k-1) - sum;
		}

		//compute L elements
		for(unsigned int i=k+1; i<=n; i++) {
			double sum=0.;
			for(unsigned int m=1; m<=(k-1); m++) sum += L(m-1,i-1)*U(k-1,m-1);
			L(k-1,i-1) = (A(k-1,i-1) - sum) / U(k-1,k-1);
		}
	}

	return L;
}

/*void Matrix::inv() {
	
}*/

const Matrix Matrix::inv() {
//This uses an LU decomposition followed by backward and forward solving for the inverse
	if(nx!=ny) {
		std::stringstream tmp;
		tmp << "Trying to invert a non-square matrix ";
		tmp << "(" << nx << "," << ny << ") !";
		throw IOException(tmp.str(), AT);
	}

	Matrix U;
	Matrix L=LU(U); //LU decomposition of current object

	//checking that the matrix can be inversed
	double product=1.;
	for(unsigned int i=0; i<nx; i++) product *= U(i,i);
	if(product==0.) {
		throw IOException("The given matrix is singular and can not be inversed", AT);
	}

	const unsigned int n = nx;
	//we solve AX=I with X=A-1. Since A=LU, then LUX = I
	//we start by forward solving LY=I with Y=UX
	Matrix Y(n, n);
	for(unsigned int i=1; i<=n; i++) {
		Y(i-1,i-1) = 1./L(i-1,i-1); //j==i
		for(unsigned int j=1; j<i; j++) { //j<i
			double sum=0.;
			for(unsigned int k=i-1; k>=1; k--) {
				sum += L(k-1,i-1) * Y(j-1,k-1);
			}
			Y(j-1,i-1) = -1./L(i-1,i-1)*sum;
		}
		for(unsigned int j=i+1; j<=n; j++) { //j>i
			Y(j-1,i-1) = 0.;
		}
	}

	//now, we backward solve UX=Y
	Matrix X(n,n);
	for(unsigned int i=n; i>=1; i--) { //lines
		for(unsigned int j=1; j<=n; j++) { //lines
			double sum=0.;
			for(unsigned int k=i+1; k<=n; k++) {
				sum += U(k-1,i-1) * X(j-1,k-1);
			}
			X(j-1,i-1) = (Y(j-1,i-1) - sum)/U(i-1,i-1);
		}
	}

	return X;
}

} //end namespace
