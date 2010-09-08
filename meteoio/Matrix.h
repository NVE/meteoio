/***********************************************************************************/
/*  Copyright 2010 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef MATRIX_H
#define MATRIX_H

#include <meteoio/IOUtils.h>
#include <meteoio/IOExceptions.h>

#include <vector>
#include <iostream>

namespace mio {

#define NOSAFECHECKS

/**
 * @class Matrix
 * @brief This class implements the basic operations on matrices.
 * Elements are access in matrix notation: that is A(1,2) represents the second element of the
 * first line. Index go from 1 to nrows/ncols.
 * 
 * It might not be the best ever such implementation, but the goal is to provide a standalone matrix class.
 * It might be later possible to chose between using the embedded implementation or to act as a
 * front end to BLAS for those who have it installed on their system.
 *
 * @author Mathias Bavay
 */
class Matrix {
	public:
		Matrix();

		/**
		* @brief A constructor that creates a matrix of a given size
		* @param rows number of rows of the new matrix
		* @param cols number of columns of the new matrix
		*/
		Matrix(const int& rows, const int& cols);
		Matrix(const unsigned int& rows, const unsigned int& cols);

		/**
		* @brief A constructor that creates a matrix filled with constant values
		* @param rows number of rows of the new matrix
		* @param cols number of columns of the new matrix
		* @param init initial value to fill the matrix with
		*/
		Matrix(const unsigned int& rows, const unsigned int& cols, const double& init);

		/**
		* @brief A constructor that creates a diagonal matrix of size n
		* @param n dimension of the new square matrix
		* @param init initial value to fill the matrix with
		*/
		Matrix(const unsigned int& n, const double& init);

		/**
		* @brief Convert the current matrix to a identity matrix of size n
		* @param n dimension of the new square matrix
		* @param init initial value to fill the matrix with
		*/
		void identity(const unsigned int& n, const double& init);

		void resize(const unsigned int& rows, const unsigned int& cols);
		void resize(const unsigned int& rows, const unsigned int& cols, const double& init);

		/**
		* @brief get the dimensions of the current object
		* @param rows number of rows of the matrix
		* @param cols number of columns of the matrix
		*/
		void size(unsigned int& rows, unsigned int& cols) const;

		/**
		* @brief free the memory and set the matrix dimensions to (0,0)
		*/
		void clear();

		double& operator ()(const unsigned int& x, const unsigned int& y);
		double operator ()(const unsigned int& x, const unsigned int& y) const;

		/**
		* @brief matrix transpose
		* @return transposed matrix
		*/
		const Matrix T();
		//void T();

		/**
		* @brief matrix invert. 
		* It first performs LU decomposition and then computes the inverse by
		* backward and forward solving of LU * A-1 = I
		* @return inversed matrix
		*/
		const Matrix inv();
		//void inv();

		/**
		* @brief matrix determinant
		* @return determinant
		*/
		double det() const;

		/**
		* @brief matrix LU decomposition. 
		* Perform LU decomposition by the Dolittle algorithm, 
		* (cf http://math.fullerton.edu/mathews/numerical/linear/dol/dol.html)
		* HACK: there is no permutation matrix, so it might not be able to give a decomposition...
		* @param L lower diagonal matrix
		* @param U upper diagonal matrix
		*/
		void LU(Matrix& L, Matrix& U) const;

		friend std::ostream& operator<<(std::ostream& os, const Matrix& data);

		Matrix& operator+=(const Matrix& rhs);
		const Matrix operator+(const Matrix& rhs);
		Matrix& operator+=(const double& rhs);
		const Matrix operator+(const double& rhs);

		Matrix& operator-=(const Matrix& rhs);
		const Matrix operator-(const Matrix& rhs);
		Matrix& operator-=(const double& rhs);
		const Matrix operator-(const double& rhs);

		Matrix& operator*=(const Matrix& rhs);
		const Matrix operator*(const Matrix& rhs);
		Matrix& operator*=(const double& rhs);
		const Matrix operator*(const double& rhs);

		Matrix& operator/=(const double& rhs);
		const Matrix operator/(const double& rhs);

		bool operator==(const Matrix&) const; ///<Operator that tests for equality
		bool operator!=(const Matrix&) const; ///<Operator that tests for inequality

		/**
		* @brief check if a matrix is the identity matrix
		* @return true if it is I
		*/
		bool isIdentity() const;
		static bool isIdentity(const Matrix& A);

	protected:
		std::vector<double> vecData;
		unsigned int ncols;
		unsigned int nrows;
};

} //end namespace

#endif
