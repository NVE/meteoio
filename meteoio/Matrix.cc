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
	nrows = ncols = 0;
}

Matrix::Matrix(const int& rows, const int& cols) {
	if(rows<0 || cols<0) {
		std::stringstream tmp;
		tmp << "Trying construct a matrix with negative dimensions: ";
		tmp << "(" << rows << "," << cols << ")";
		throw IOException(tmp.str(), AT);
	}
	nrows = ncols = 0;
	resize((unsigned)rows,(unsigned)cols);
}

Matrix::Matrix(const unsigned int& rows, const unsigned int& cols) {
	nrows = ncols = 0;
	resize(rows,cols);
}

Matrix::Matrix(const unsigned int& rows, const unsigned int& cols, const double& init) {
	nrows = ncols = 0;
	resize(rows,cols,init);
}

Matrix::Matrix(const unsigned int& n, const double& init) {
	nrows = ncols = 0;
	resize(n,n,0.);
	for(unsigned int ii=1; ii<=n; ii++) operator()(ii,ii) = init;
}

void Matrix::resize(const unsigned int& rows, const unsigned int& cols) {
	clear();

	if ((rows > 0) && (cols > 0)) {
		vecData.resize(rows*cols);
		ncols = cols;
		nrows = rows;
	} else {
		throw IndexOutOfBoundsException("Can not resize a matrix to negative sizes!", AT);
	}
}

void Matrix::resize(const unsigned int& rows, const unsigned int& cols, const double& init) {
	resize(rows, cols);

	for (unsigned int ii=1; ii<=nrows; ii++) {
		for (unsigned int jj=1; jj<=ncols; jj++) {
			operator()(ii,jj) = init;
		}
	}
}

void Matrix::size(unsigned int& rows, unsigned int& cols) const{
	rows=nrows;
	cols=ncols;
}

void Matrix::clear() {
	vecData.clear();
	nrows=ncols=0;
}

double& Matrix::operator ()(const unsigned int& i, const unsigned int& j) {
#ifndef NOSAFECHECKS
	if ((i<1) || (i > nrows) || (j<1) || (j > ncols)) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif
	return vecData[(j-1) + (i-1)*ncols];
}

const double Matrix::operator ()(const unsigned int& i, const unsigned int& j) const {
#ifndef NOSAFECHECKS
	if ((i<1) || (i > nrows) || (j<1) || (j > ncols)) {
		throw IndexOutOfBoundsException("", AT);
	}
#endif
	return vecData[(j-1) + (i-1)*ncols];
}

std::ostream& operator<<(std::ostream& os, const Matrix& data) {
	const unsigned int wd=6;
	os << "\n┌ ";
	for(unsigned int jj=1; jj<=(data.ncols*(wd+1)); jj++)
		os << " ";
	os << " ┐\n";
	for(unsigned int ii=1; ii<=data.nrows; ii++) {
		os << "│ ";
		for (unsigned int jj=1; jj<=data.ncols; jj++) {
			os << std::setw(wd) << std::fixed << std::setprecision(2) << data(ii,jj) << " ";
		}
		os << " │\n";
	}
	os << "└ ";
	for(unsigned int jj=1; jj<=(data.ncols*(wd+1)); jj++)
		os << " ";
	os << " ┘\n";
	return os;
}

bool Matrix::operator==(const Matrix& in) const {
	unsigned int in_nrows, in_ncols;
	in.size(in_nrows, in_ncols);

	if(nrows!=in_nrows || ncols!=in_ncols)
		return false;

	for(unsigned int i=1; i<=nrows; i++) {
		for(unsigned int j=1; j<=ncols; j++) {
			if( operator()(i,j) != in(i,j) ) return false;
		}
	}

	return true;
}

bool Matrix::operator!=(const Matrix& in) const {
	return !(*this==in);
}

Matrix& Matrix::operator+=(const Matrix& rhs) {
	//check dimensions compatibility
	if(nrows!=rhs.nrows || ncols!=rhs.ncols) {
		std::stringstream tmp;
		tmp << "Trying to add two matrix with incompatible dimensions: ";
		tmp << "(" << nrows << "," << ncols << ") * ";
		tmp << "(" << rhs.nrows << "," << rhs.ncols << ")";
		throw IOException(tmp.str(), AT);
	}

	//fill sum matrix
	for(unsigned int i=1; i<=nrows; i++) {
		for(unsigned int j=1; j<=ncols; j++) {
			operator()(i,j) += rhs(i,j);
		}
	}

	return *this;
}

const Matrix Matrix::operator+(const Matrix& rhs) {
	Matrix result = *this;
	result += rhs; //already implemented

	return result;
}

Matrix& Matrix::operator+=(const double& rhs) {
	//fill sum matrix
	for(unsigned int i=1; i<=nrows; i++) {
		for(unsigned int j=1; j<=ncols; j++) {
			operator()(i,j) += rhs;
		}
	}

	return *this;
}

const Matrix Matrix::operator+(const double& rhs) {
	Matrix result = *this;
	result += rhs; //already implemented

	return result;
}

Matrix& Matrix::operator-=(const Matrix& rhs) {
	//check dimensions compatibility
	if(nrows!=rhs.nrows || ncols!=rhs.ncols) {
		std::stringstream tmp;
		tmp << "Trying to substract two matrix with incompatible dimensions: ";
		tmp << "(" << nrows << "," << ncols << ") * ";
		tmp << "(" << rhs.nrows << "," << rhs.ncols << ")";
		throw IOException(tmp.str(), AT);
	}

	//fill sum matrix
	for(unsigned int i=1; i<=nrows; i++) {
		for(unsigned int j=1; j<=ncols; j++) {
			operator()(i,j) -= rhs(i,j);
		}
	}

	return *this;
}

const Matrix Matrix::operator-(const Matrix& rhs) {
	Matrix result = *this;
	result -= rhs; //already implemented

	return result;
}

Matrix& Matrix::operator-=(const double& rhs) {
	*this += -rhs;

	return *this;
}

const Matrix Matrix::operator-(const double& rhs) {
	Matrix result = *this;
	result += -rhs; //already implemented

	return result;
}

Matrix& Matrix::operator*=(const Matrix& rhs) {
	//check dimensions compatibility
	if(ncols!=rhs.nrows) {
		std::stringstream tmp;
		tmp << "Trying to multiply two matrix with incompatible dimensions: ";
		tmp << "(" << nrows << "," << ncols << ") * ";
		tmp << "(" << rhs.nrows << "," << rhs.ncols << ")";
		throw IOException(tmp.str(), AT);
	}

	//create new matrix
	Matrix result(nrows,rhs.ncols);

	//fill product matrix
	for(unsigned int i=1; i<=result.nrows; i++) {
		for(unsigned int j=1; j<=result.ncols; j++) {
			double sum=0.;
			for(unsigned int idx=1; idx<=ncols; idx++) {
				sum+= operator()(i,idx) * rhs(idx,j);
			}
			result(i,j) = sum;
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

Matrix& Matrix::operator*=(const double& rhs) {
	for(unsigned int i=1; i<=nrows; i++) {
		for(unsigned int j=1; j<=ncols; j++) {
			operator()(i,j) *= rhs;
		}
	}

	return *this;
}

const Matrix Matrix::operator*(const double& rhs) {
	Matrix result = *this;
	result *= rhs; //already implemented

	return result;
}

Matrix& Matrix::operator/=(const double& rhs) {
	*this *= (1./rhs);
	return *this;
}

const Matrix Matrix::operator/(const double& rhs) {
	Matrix result = *this;
	result *= 1./rhs; //already implemented

	return result;
}

/*void Matrix::T() {

}*/

const Matrix Matrix::T() {
//other possibility: create a "transpose" flag that simply swaps the data reading...
	Matrix result(ncols, nrows);
	for(unsigned int i=1; i<=result.nrows; i++) {
		for(unsigned int j=1; j<=result.ncols; j++) {
			result(i,j) = operator()(j,i);
		}
	}
	return result;
}

double Matrix::det() const {
	if(nrows!=ncols) {
		std::stringstream tmp;
		tmp << "Trying to calculate the determinant of a non-square matrix ";
		tmp << "(" << nrows << "," << ncols << ") !";
		throw IOException(tmp.str(), AT);
	}
	Matrix U;
	LU(U);

	double product=1.;
	for(unsigned int i=1; i<=nrows; i++) product *= U(i,i);
	
	return product;
}

const Matrix Matrix::LU(Matrix& U) const {
//Dolittle algorithm, cf http://math.fullerton.edu/mathews/numerical/linear/dol/dol.html
//HACK: there is no permutation matrix, so it might not be able to give a decomposition...
//This is implemented for matrix storage from 0 to n-1
	if(nrows!=ncols) {
		std::stringstream tmp;
		tmp << "Trying to calculate the LU decomposition of a non-square matrix ";
		tmp << "(" << nrows << "," << ncols << ") !";
		throw IOException(tmp.str(), AT);
	}

	const unsigned int n = nrows;
	U = *this;
	const Matrix& A = *this;
	Matrix L(n, 1.); //initialized as diagonal matrix, then populated

	for(unsigned int k=1; k<=n; k++) {
		//compute U elements
		for(unsigned int j=1; j<k; j++) {
			U(k,j) = 0.;
		}
		for(unsigned int j=k; j<=n; j++) {
			double sum=0.;
			for(unsigned int m=1; m<=(k-1); m++) sum += L(k,m)*U(m,j);
			U(k,j) = A(k,j) - sum;
		}

		//compute L elements
		for(unsigned int i=k+1; i<=n; i++) {
			double sum=0.;
			for(unsigned int m=1; m<=(k-1); m++) sum += L(i,m)*U(m,k);
			L(i,k) = (A(i,k) - sum) / (U(k,k)+1e-12);
		}
	}

	return L;
}

/*void Matrix::inv() {
	
}*/

const Matrix Matrix::inv() {
//This uses an LU decomposition followed by backward and forward solving for the inverse
	if(nrows!=ncols) {
		std::stringstream tmp;
		tmp << "Trying to invert a non-square matrix ";
		tmp << "(" << nrows << "," << ncols << ") !";
		throw IOException(tmp.str(), AT);
	}
	const unsigned int n = nrows;

	Matrix U;
	Matrix L=LU(U); //LU decomposition of current object

	//checking that the matrix can be inversed
	double product=1.;
	for(unsigned int i=1; i<=n; i++) product *= U(i,i);
	if(product==0.) {
		throw IOException("The given matrix is singular and can not be inversed", AT);
	}

	//we solve AX=I with X=A-1. Since A=LU, then LUX = I
	//we start by forward solving LY=I with Y=UX
	Matrix Y(n, n);
	for(unsigned int i=1; i<=n; i++) {
		Y(i,i) = 1./L(i,i); //j==i
		for(unsigned int j=1; j<i; j++) { //j<i
			double sum=0.;
			for(unsigned int k=i-1; k>=1; k--) {
				sum += L(i,k) * Y(k,j);
			}
			Y(i,j) = -1./(L(i,i)+1e-12) * sum;
		}
		for(unsigned int j=i+1; j<=n; j++) { //j>i
			Y(i,j) = 0.;
		}
	}

	//now, we backward solve UX=Y
	Matrix X(n,n);
	for(unsigned int i=n; i>=1; i--) { //lines
		for(unsigned int j=1; j<=n; j++) { //lines
			double sum=0.;
			for(unsigned int k=i+1; k<=n; k++) {
				sum += U(i,k) * X(k,j);
			}
			X(i,j) = (Y(i,j) - sum) / (U(i,i)+1e-12);
		}
	}

	return X;
}

bool Matrix::isIdentity() const {
	const double eps=1e-6;

	if(nrows!=ncols) {
		std::stringstream tmp;
		tmp << "A non-square matrix ";
		tmp << "(" << nrows << "," << ncols << ") can not be the identity matrix!";
		throw IOException(tmp.str(), AT);
	}
	
	bool is_identity=true;
	for(unsigned int i=1; i<=nrows; i++) {
		for(unsigned int j=1; j<=ncols; j++) {
			const double val = operator()(i,j);
			if(i!=j) {
				if(IOUtils::checkEpsilonEquality(val,0.,eps)==false) {
					is_identity=false;
					break;
				}
			} else {
				if(IOUtils::checkEpsilonEquality(val,1.,eps)==false) {
					is_identity=false;
					break;
				}
			}
		}
	}

	return is_identity;
}

const bool Matrix::isIdentity(const Matrix& A) {
	return A.isIdentity();
}

} //end namespace
