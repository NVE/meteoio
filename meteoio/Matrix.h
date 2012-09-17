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
 * If the compilation flag NOSAFECHECKS is used, bounds check is turned off (leading to increased performances).
 *
 * @ingroup data_str
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
		* @brief Copy constructor
		* @param init matrix to copy
		*/
		Matrix(const Matrix& init);

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

		/**
		* @brief fill the matrix with random numbers.
		* @param range range of the randoms numbers (they will be between -range and +range)
		*/
		void random(const double& range);

		double& operator ()(const unsigned int& x, const unsigned int& y);
		double operator ()(const unsigned int& x, const unsigned int& y) const;

		/**
		* @brief Converts a 1x1 matrix to a scalar.
		* @return scalar value
		*/
		double scalar() const;
		static double scalar(const Matrix& m);

		/**
		* @brief Dot product.
		* @return scalar value
		*/
		static double dot(const Matrix& A, const Matrix& B);

		/**
		* @brief matrix transpose.
		* @return transposed matrix
		*/
		Matrix getT() const;
		static Matrix T(const Matrix& m);
		void T();

		/**
		* @brief matrix invert.
		* It first performs LU decomposition and then computes the inverse by
		* backward and forward solving of LU * A-1 = I
		* see Press, William H.; Flannery, Brian P.; Teukolsky, Saul A.; Vetterling, William T. (1992), "LU Decomposition and Its Applications", Numerical Recipes in FORTRAN: The Art of Scientific Computing (2nd ed.), Cambridge University Press, pp. 34–42
		* @return inversed matrix
		*/
		Matrix getInv() const;
		void inv();

		/**
		* @brief matrix solving for A·X=B.
		* It first performs LU decomposition and then solves A·X=B by
		* backward and forward solving of LU * X = B
		* @param A A matrix
		* @param B B matrix
		* @return solution matrix
		*/
		static Matrix solve(const Matrix& A, const Matrix& B);

		/**
		* @brief matrix solving for A·X=B.
		* It first performs LU decomposition and then solves A·X=B by
		* backward and forward solving of LU * X = B
		* @param A A matrix
		* @param B B matrix
		* @param X solution matrix
		*/
		static void solve(const Matrix& A, const Matrix& B, Matrix& X);

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
		* @return false if the decomposition can not be performed (division by zero)
		*/
		bool LU(Matrix& L, Matrix& U) const;

		/**
		* @brief matrix partial pivoting.
		* This reorders the rows so that each diagonal element is the maximum in its column
		* (see https://secure.wikimedia.org/wikipedia/en/wiki/Pivot_element)
		* @param pivot_idx new indices (to apply when solving A * X = B, for example
		*/
		void partialPivoting(std::vector<unsigned int>& pivot_idx);
		void partialPivoting();

		void maximalPivoting();

		/**
		* @brief matrix bidiagonalization
		* This uses Householder's reduction, see Golub, 1970.
		* (see https://secure.wikimedia.org/wikipedia/en/wiki/Bidiagonalization)
		*/
		//void bidiagonalize();

		friend std::ostream& operator<<(std::ostream& os, const Matrix& data);

		Matrix& operator+=(const Matrix& rhs);
		const Matrix operator+(const Matrix& rhs) const;
		Matrix& operator+=(const double& rhs);
		const Matrix operator+(const double& rhs) const;

		Matrix& operator-=(const Matrix& rhs);
		const Matrix operator-(const Matrix& rhs) const;
		Matrix& operator-=(const double& rhs);
		const Matrix operator-(const double& rhs) const;

		Matrix& operator*=(const Matrix& rhs);
		const Matrix operator*(const Matrix& rhs) const;
		Matrix& operator*=(const double& rhs);
		const Matrix operator*(const double& rhs) const;

		Matrix& operator/=(const double& rhs);
		const Matrix operator/(const double& rhs) const;

		bool operator==(const Matrix&) const; ///<Operator that tests for equality
		bool operator!=(const Matrix&) const; ///<Operator that tests for inequality

		/**
		* @brief check if a matrix is the identity matrix
		* @return true if it is I
		*/
		bool isIdentity() const;
		static bool isIdentity(const Matrix& A);

		static const double epsilon, epsilon_mtr;

	protected:
		std::vector<double> vecData;
		unsigned int ncols;
		unsigned int nrows;

		unsigned int findMaxInCol(const unsigned int &col);
		unsigned int findMaxInRow(const unsigned int &row);
		void swapRows(const unsigned int &i1, const unsigned int &i2);
};

} //end namespace

#endif
