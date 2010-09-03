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
#ifndef MATRIX_H
#define MATRIX_H

#include <meteoio/IOUtils.h>
#include <meteoio/IOExceptions.h>

#include <vector>
#include <iostream>

namespace mio {

#define NOSAFECHECKS

//HACK: indices are WRONG!! we should start at 1, not 0!!
class Matrix {
	public:
		Matrix();

		/**
		* A constructor that creates an array of a given size
		* @param anx number of columns of the new array
		* @param any number of rows of the new array
		*/
		Matrix(const int& anx, const int& any);

		/**
		* A constructor that creates an array filled with constant values
		* @param anx number of columns of the new array
		* @param any number of rows of the new array
		* @param init initial value to fill the array with
		*/
		Matrix(const unsigned int& anx, const unsigned int& any, const double& init);

		/**
		* A constructor that creates a diagonal matrix of size n
		* @param n number of columns of the new array
		* @param any number of rows of the new array
		*/
		Matrix(const unsigned int& n, const double& init);

		void resize(const unsigned int& nx, const unsigned int& ny);
		void resize(const unsigned int& nx, const unsigned int& ny, const double& init);
		void size(unsigned int& nx, unsigned int& ny) const;
		void clear();

		double& operator ()(const unsigned int& x, const unsigned int& y);
		const double operator ()(const unsigned int& x, const unsigned int& y) const;

		//void T();
		const Matrix T();
		//void inv();
		const Matrix inv();
		double det() const;
		const Matrix LU(Matrix& U) const;

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

		bool isIdentity() const;
		static const bool isIdentity(const Matrix& A);

	protected:
		std::vector<double> vecData;
		unsigned int nx;
		unsigned int ny;
};

} //end namespace

#endif
