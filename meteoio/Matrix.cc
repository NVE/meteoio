/***********************************************************************************/
/*  Copyright 2011 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <time.h> //needed for random()
#include <cmath> //needed for fabs()

namespace mio {

const double Matrix::epsilon = 1e-9; //for considering a determinant to be zero, etc
const double Matrix::epsilon_mtr = 1e-6; //for comparing two matrix

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
	identity(n, init);
}

Matrix::Matrix(const Matrix& init) {
	unsigned int tmprows, tmpcols;
	init.size(tmprows, tmpcols);
	nrows = ncols = 0;
	resize(tmprows, tmpcols);

	for (unsigned int ii=1; ii<=nrows; ii++) {
		for (unsigned int jj=1; jj<=ncols; jj++) {
			operator()(ii,jj) = init(ii,jj);
		}
	}
}

void Matrix::identity(const unsigned int& n, const double& init) {
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

void Matrix::random(const double& range) {
	srand((unsigned)time(0));
	for(unsigned int i=1; i<=nrows; i++) {
		for(unsigned int j=1; j<=ncols; j++) {
			operator()(i,j) = (double)rand()/(double)RAND_MAX*range;
		}
	}
}

double& Matrix::operator ()(const unsigned int& i, const unsigned int& j) {
#ifndef NOSAFECHECKS
	if ((i<1) || (i > nrows) || (j<1) || (j > ncols)) {
		std::stringstream ss;
		ss << "Trying to access matrix[" << i << "," << j << "]";
		throw IndexOutOfBoundsException(ss.str(), AT);
	}
#endif
	return vecData[(j-1) + (i-1)*ncols];
}

double Matrix::operator ()(const unsigned int& i, const unsigned int& j) const {
#ifndef NOSAFECHECKS
	if ((i<1) || (i > nrows) || (j<1) || (j > ncols)) {
		std::stringstream ss;
		ss << "Trying to access matrix[" << i << "," << j << "]";
		throw IndexOutOfBoundsException(ss.str(), AT);
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
			//if( operator()(i,j) != in(i,j) ) return false;
			if( !IOUtils::checkEpsilonEquality( operator()(i,j) , in(i,j), epsilon_mtr) ) return false;
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

const Matrix Matrix::operator+(const Matrix& rhs) const {
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

const Matrix Matrix::operator+(const double& rhs) const {
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

const Matrix Matrix::operator-(const Matrix& rhs) const {
	Matrix result = *this;
	result -= rhs; //already implemented

	return result;
}

Matrix& Matrix::operator-=(const double& rhs) {
	*this += -rhs;

	return *this;
}

const Matrix Matrix::operator-(const double& rhs) const {
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
				sum += operator()(i,idx) * rhs(idx,j);
			}
			result(i,j) = sum;
		}
	}

	*this = result;
	return *this;
}

const Matrix Matrix::operator*(const Matrix& rhs) const {
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

const Matrix Matrix::operator*(const double& rhs) const {
	Matrix result = *this;
	result *= rhs; //already implemented

	return result;
}

Matrix& Matrix::operator/=(const double& rhs) {
	*this *= (1./rhs);
	return *this;
}

const Matrix Matrix::operator/(const double& rhs) const {
	Matrix result = *this;
	result *= 1./rhs; //already implemented

	return result;
}

double Matrix::scalar(const Matrix& m) {
	return m.scalar();
}

double Matrix::scalar() const {
	if(ncols!=1 || nrows!=1) {
		std::stringstream tmp;
		tmp << "Trying to get scalar value of a non (1x1) matrix ";
		tmp << "(" << nrows << "," << ncols << ") !";
		throw IOException(tmp.str(), AT);
	}
	return operator()(1,1);
}

double Matrix::dot(const Matrix& A, const Matrix& B) {
	unsigned int Acols, Arows, Bcols, Brows;
	A.size(Arows, Acols);
	B.size(Brows, Bcols);

	if(Acols!=1 || Bcols!=1) {
		std::stringstream tmp;
		tmp << "Trying to get dot product of non vector matrix ";
		tmp << "(" << Arows << "," << Acols << ") · ";
		tmp << "(" << Brows << "," << Bcols << ") · ";
		throw IOException(tmp.str(), AT);
	}
	if(Arows!=Brows) {
		std::stringstream tmp;
		tmp << "Trying to get dot product of incompatible matrix ";
		tmp << "(" << Arows << "," << Acols << ") · ";
		tmp << "(" << Brows << "," << Bcols << ") · ";
		throw IOException(tmp.str(), AT);
	}

	double sum=0.;
	for(unsigned int i=1; i<=Arows; i++) {
		sum += A(i,1)*B(i,1);
	}

	return sum;
}

Matrix Matrix::T(const Matrix& m) {
	return m.getT();
}

Matrix Matrix::getT() const {
//other possibility: create a "transpose" flag that simply swaps the data reading...
	Matrix result(ncols, nrows);
	for(unsigned int i=1; i<=result.nrows; i++) {
		for(unsigned int j=1; j<=result.ncols; j++) {
			result(i,j) = operator()(j,i);
		}
	}
	return result;
}

void Matrix::T() {
	Matrix tmp(*this);
	*this = tmp.getT();
}

double Matrix::det() const {
	if(nrows!=ncols) {
		std::stringstream tmp;
		tmp << "Trying to calculate the determinant of a non-square matrix ";
		tmp << "(" << nrows << "," << ncols << ") !";
		throw IOException(tmp.str(), AT);
	}
	Matrix L,U;
	if(LU(L,U)==false) return 0.;

	double product=1.;
	for(unsigned int i=1; i<=nrows; i++) product *= U(i,i);

	return product;
}

bool Matrix::LU(Matrix& L, Matrix& U) const {
//Dolittle algorithm, cf http://math.fullerton.edu/mathews/numerical/linear/dol/dol.html
//HACK: there is no permutation matrix, so it might not be able to give a decomposition...
	if(nrows!=ncols) {
		std::stringstream tmp;
		tmp << "Trying to calculate the LU decomposition of a non-square matrix ";
		tmp << "(" << nrows << "," << ncols << ") !";
		throw IOException(tmp.str(), AT);
	}

	const unsigned int n = nrows;
	U.clear();
	U = *this;
	L.identity(n, 1.); //initialized as identity matrix, then populated
	const Matrix& A = *this;

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

		if( k<n && IOUtils::checkEpsilonEquality(U(k,k), 0., epsilon) ) return false; //we can not compute L
		//compute L elements
		for(unsigned int i=k+1; i<=n; i++) {
			double sum=0.;
			for(unsigned int m=1; m<=(k-1); m++) sum += L(i,m)*U(m,k);
			L(i,k) = (A(i,k) - sum) / U(k,k);
		}
	}
	return true;
}

Matrix Matrix::getInv() const {
//This uses an LU decomposition followed by backward and forward solving for the inverse
//See for example Press, William H.; Flannery, Brian P.; Teukolsky, Saul A.; Vetterling, William T. (1992), "LU Decomposition and Its Applications", Numerical Recipes in FORTRAN: The Art of Scientific Computing (2nd ed.), Cambridge University Press, pp. 34–42
	if(nrows!=ncols) {
		std::stringstream tmp;
		tmp << "Trying to invert a non-square matrix ";
		tmp << "(" << nrows << "," << ncols << ") !";
		throw IOException(tmp.str(), AT);
	}
	const unsigned int n = nrows;

	Matrix U;
	Matrix L;
	if(LU(L, U)==false) {
		throw IOException("LU decomposition of given matrix not possible", AT);
	}

	//we solve AX=I with X=A-1. Since A=LU, then LUX = I
	//we start by forward solving LY=I with Y=UX
	Matrix Y(n, n);
	for(unsigned int i=1; i<=n; i++) {
		if(IOUtils::checkEpsilonEquality(L(i,i), 0., epsilon)) {
			throw IOException("The given matrix can not be inversed", AT);
		}
		Y(i,i) = 1./L(i,i); //j==i
		for(unsigned int j=1; j<i; j++) { //j<i
			double sum=0.;
			for(unsigned int k=i-1; k>=1; k--) { //equivalent to 1 -> i-1
				sum += L(i,k) * Y(k,j);
			}
			Y(i,j) = -1./L(i,i) * sum;
		}
		for(unsigned int j=i+1; j<=n; j++) { //j>i
			Y(i,j) = 0.;
		}
	}

	//now, we backward solve UX=Y
	Matrix X(n,n);
	for(unsigned int i=n; i>=1; i--) { //lines
		if(IOUtils::checkEpsilonEquality(U(i,i), 0., epsilon)) { //HACK: actually, only U(n,n) needs checking
			throw IOException("The given matrix is singular and can not be inversed", AT);
		}
		for(unsigned int j=1; j<=n; j++) { //lines
			double sum=0.;
			for(unsigned int k=i+1; k<=n; k++) {
				sum += U(i,k) * X(k,j);
			}
			X(i,j) = (Y(i,j) - sum) / U(i,i);
		}
	}

	return X;
}

void Matrix::inv() {
//same as getInv() const but we write the final result on top of the input matrix
	if(nrows!=ncols) {
		std::stringstream tmp;
		tmp << "Trying to invert a non-square matrix ";
		tmp << "(" << nrows << "," << ncols << ") !";
		throw IOException(tmp.str(), AT);
	}
	const unsigned int n = nrows;

	Matrix U;
	Matrix L;
	if(LU(L, U)==false) {
		throw IOException("LU decomposition of given matrix not possible", AT);
	}

	//we solve AX=I with X=A-1. Since A=LU, then LUX = I
	//we start by forward solving LY=I with Y=UX
	Matrix Y(n, n);
	for(unsigned int i=1; i<=n; i++) {
		if(IOUtils::checkEpsilonEquality(L(i,i), 0., epsilon)) {
			throw IOException("The given matrix can not be inversed", AT);
		}
		Y(i,i) = 1./L(i,i); //j==i
		for(unsigned int j=1; j<i; j++) { //j<i
			double sum=0.;
			for(unsigned int k=i-1; k>=1; k--) { //equivalent to 1 -> i-1
				sum += L(i,k) * Y(k,j);
			}
			Y(i,j) = -1./L(i,i) * sum;
		}
		for(unsigned int j=i+1; j<=n; j++) { //j>i
			Y(i,j) = 0.;
		}
	}

	//now, we backward solve UX=Y
	Matrix& X = *this; //we write the solution over the input matrix
	for(unsigned int i=n; i>=1; i--) { //lines
		if(IOUtils::checkEpsilonEquality(U(i,i), 0., epsilon)) { //HACK: actually, only U(n,n) needs checking
			throw IOException("The given matrix is singular and can not be inversed", AT);
		}
		for(unsigned int j=1; j<=n; j++) { //lines
			double sum=0.;
			for(unsigned int k=i+1; k<=n; k++) {
				sum += U(i,k) * X(k,j);
			}
			X(i,j) = (Y(i,j) - sum) / U(i,i);
		}
	}
}

void Matrix::solve(const Matrix& A, const Matrix& B, Matrix& X) {
//This uses an LU decomposition followed by backward and forward solving for A·X=B
	unsigned int Anrows,Ancols, Bnrows, Bncols;
	A.size(Anrows, Ancols);
	if(Anrows!=Ancols) {
		std::stringstream tmp;
		tmp << "Trying to solve A·X=B with A non square matrix ";
		tmp << "(" << Anrows << "," << Ancols << ") !";
		throw IOException(tmp.str(), AT);
	}
	B.size(Bnrows, Bncols);
	if(Anrows!=Bnrows)  {
		std::stringstream tmp;
		tmp << "Trying to solve A·X=B with A and B of incompatible dimensions ";
		tmp << "(" << Anrows << "," << Ancols << ") and (";
		tmp << "(" << Bnrows << "," << Bncols << ") !";
		throw IOException(tmp.str(), AT);
	}
	const unsigned int n = Anrows;
	const unsigned int m = Bncols;

	Matrix U;
	Matrix L;
	if(A.LU(L, U)==false) {
		throw IOException("LU decomposition of A matrix not possible", AT);
	}

	//we solve AX=B. Since A=LU, then LUX = B
	//we start by forward solving LY=B with Y=UX
	Matrix Y(n, m);
	for(unsigned int i=1; i<=n; i++) {
		if(IOUtils::checkEpsilonEquality(L(i,i), 0., epsilon)) {
			throw IOException("The given matrix can not be inversed", AT);
		}
		for(unsigned int j=1; j<=m; j++) {
			double sum=0.;
			for(unsigned int k=1; k<i; k++) {
				sum += L(i,k) * Y(k,j);
			}
			Y(i,j) = (B(i,j) - sum) / L(i,i);
		}
	}

	//now, we backward solve UX=Y
	X.resize(n,m); //we need to ensure that X has the correct dimensions
	for(unsigned int i=n; i>=1; i--) { //lines
		if(IOUtils::checkEpsilonEquality(U(i,i), 0., epsilon)) { //HACK: actually, only U(n,n) needs checking
			std::stringstream ss;
			ss << "The given " << Ancols << "*" << Anrows << " matrix is singular and can not be inversed";
			throw IOException(ss.str(), AT);
		}
		for(unsigned int j=1; j<=m; j++) {
			double sum = 0.;
			for(unsigned int k=i+1; k<=n; k++) {
				sum += U(i,k) * X(k,j);
			}
			X(i,j) = (Y(i,j) - sum) / U(i,i);
		}
	}
}

Matrix Matrix::solve(const Matrix& A, const Matrix& B) {
//This uses an LU decomposition followed by backward and forward solving for A·X=B
	Matrix X;
	solve(A, B, X);
	return X;
}

bool Matrix::isIdentity() const {
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
				if(IOUtils::checkEpsilonEquality(val,0.,epsilon_mtr)==false) {
					is_identity=false;
					break;
				}
			} else {
				if(IOUtils::checkEpsilonEquality(val,1.,epsilon_mtr)==false) {
					is_identity=false;
					break;
				}
			}
		}
	}

	return is_identity;
}

bool Matrix::isIdentity(const Matrix& A) {
	return A.isIdentity();
}

void Matrix::partialPivoting(std::vector<unsigned int>& pivot_idx) {
	pivot_idx.clear();

	//bad luck: if a row has several elements that are max of their columns,
	//we don't optimize its position. Ie: we can end up with a small element
	//on the diagonal
	for(unsigned int j=1; j<=ncols; j++) {
		const unsigned int old_i = j;
		const unsigned int new_i = findMaxInCol(j);
		if(new_i!=j) { //ie: pivoting needed
			swapRows(old_i, new_i);
			pivot_idx.push_back(new_i);
		} else
			pivot_idx.push_back(old_i);
	}
}

void Matrix::partialPivoting() {
	std::vector<unsigned int> pivot_idx;
	partialPivoting(pivot_idx);
}

void Matrix::maximalPivoting() {
	std::vector<unsigned int> pivot_idx;
	Matrix tmp( *this );

	for(unsigned int i=1; i<=nrows; i++) {
		const double scale = operator()(i,findMaxInRow(i));
		for(unsigned int j=1; j<=ncols; j++) {
			operator()(i,j) /= scale;
		}
	}
	tmp.partialPivoting(pivot_idx);

	//pivot on original matrix //HACK: not finished yet!
	throw IOException("Method not implemented yet!!", AT);
}

/*void Matrix::bidiagonalize() {
	//Matrix e(1,ncols);
	std::vector<double> e(ncols+1); //so we remain compatible with matrix index
	double g=0., x=0.;

	for(unsigned int i=1; i<=ncols; i++) {
		e[i]=g; s=0.; l=i+1;
		for(unsigned int j=i; j<=m; j++) s += ( operator()(i,j)*operator()(i,j) );
	}
}*/

//return the index of the line containing the highest absolute value at column col
unsigned int Matrix::findMaxInCol(const unsigned int &col) {
	unsigned int row_idx = 0;
	double max_val=0.;

	for(unsigned int i=1; i<=nrows; i++) {
		const double val = fabs( operator()(i,col) );
		if( val>max_val) {
			max_val=val;
			row_idx=i;
		}
	}
	return row_idx;
}

//return the index of the line containing the highest absolute value at column col
unsigned int Matrix::findMaxInRow(const unsigned int &row) {
	unsigned int col_idx = 0;
	double max_val=0.;

	for(unsigned int j=1; j<=ncols; j++) {
		const double val = fabs( operator()(row,j) );
		if( val>max_val) {
			max_val=val;
			col_idx=j;
		}
	}
	return col_idx;
}


void Matrix::swapRows(const unsigned int &i1, const unsigned int &i2) {
	for(unsigned int j=1; j<=ncols; j++) {
		const double tmp = operator()(i2,j);
		operator()(i2,j) = operator()(i1,j);
		operator()(i1,j) = tmp;
	}
}

} //end namespace
