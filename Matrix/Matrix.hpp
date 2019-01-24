#pragma once

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <deque>
#include <algorithm>
#include <string>

#define CHECK_BOUNDS(_row, _col) if(!InBounds(_row, _col)) return

#define _VOID 
#define _BOOL false
#define _INT 0

typedef std::deque<std::deque<double>> doubleMatrix;
typedef unsigned int uint;

class Matrix {
public:
	Matrix(uint rows, uint cols);
	Matrix(const Matrix& other);

	~Matrix();

	void PrintMatrix();
	friend bool Delta(uint right, uint left);


	// Constructional functions
	friend void Zero(Matrix& mat);
	friend void Identity(Matrix& mat);
	

	// Setters / Getters

	void SetNumber(uint row, uint col, double val);
	double GetNumber(uint row, uint col) const;

	uint GetRows() const;
	uint GetColumns() const;

	bool IsInvertable();
	bool IsSquare();


	// Transformation / arithmetic functions

	void Transpose();
	void MultiplyRow(uint row, double factor);
	void SwapRows(uint left, uint right);
	void AddMultiplesToRow(uint base, uint target, double factor);

	void TransformToEchelonForm();
	void TransformToReducedEchelonForm();

	friend double det(const Matrix& mat);
	friend Matrix mul(const Matrix& left, const Matrix& right);
	
	// Operator overloads

	Matrix operator+(const Matrix& other);
	Matrix operator-(const Matrix& other);
	Matrix operator*(const double& other);
	Matrix operator/(const double& other);

	void operator+=(const Matrix& other);
	void operator-=(const Matrix& other);
	void operator*=(const double& other);
	void operator/=(const double& other);

	Matrix operator*(const Matrix& other);
	void operator*=(const Matrix& other);

	friend std::ostream& operator<<(std::ostream& stream, const Matrix& other);
private:
	doubleMatrix m_matrix;
	uint m_rows, m_cols;

	bool InBounds(uint rows, uint cols) const;
	bool DimensionsFitting(const Matrix& left, const Matrix& right);
	void Resize(uint rows, uint cols, doubleMatrix& matrix);


};


////////////////////////////////////////////////////
/// \brief Constructs the matrix by setting the dimensions and filling it up with zeros
///
/// \param rows Rows of the matrix
/// \param cols Cols of the matrix
///
////////////////////////////////////////////////////
inline Matrix::Matrix(uint rows, uint cols) :
	m_rows(rows),
	m_cols(cols)
{
	Resize(rows, cols, m_matrix);
}


////////////////////////////////////////////////////
/// \brief Copy constructor
///
/// \param other The matrix to set this one to
///
////////////////////////////////////////////////////
inline Matrix::Matrix(const Matrix& other)
{
	*this = other;
}


Matrix::~Matrix() { /* EMPTY */ }


////////////////////////////////////////////////////
/// \brief Fills the matrix up with zeros
///
////////////////////////////////////////////////////
inline void Zero(Matrix& mat)
{
	mat = Matrix(mat.m_rows, mat.m_cols);
}

////////////////////////////////////////////////////
/// \brief Tries to create the identity matrix
/// The matrix must be square, or it will not work
///
////////////////////////////////////////////////////
inline void Identity(Matrix& mat)
{
	if (!mat.IsSquare())
		return;

	Zero(mat);
	for (uint k = 0; k < mat.m_cols; k++)
		mat.m_matrix[k][k] = 1;
}


////////////////////////////////////////////////////
/// \brief Sets the number at a given position to a given value
/// 
/// \param row The row of the number to change
/// \param col The column of the number to change
/// \param val The new value of the number
///
////////////////////////////////////////////////////
inline void Matrix::SetNumber(uint row, uint col, double val)
{
	CHECK_BOUNDS(row, col) _VOID;

	m_matrix[row][col] = val;
}


////////////////////////////////////////////////////
/// \brief Gets the number at a given position to a given value
///
/// \return The number at the position
///
////////////////////////////////////////////////////
inline double Matrix::GetNumber(uint row, uint col) const
{
	CHECK_BOUNDS(row, col) _INT;

	return m_matrix[row][col];
}


////////////////////////////////////////////////////
/// \brief Returns the number of rows the matrix has
///
////////////////////////////////////////////////////
uint Matrix::GetRows() const
{
	return m_rows;
}

////////////////////////////////////////////////////
/// \brief Returns the number of columns the matrix has
///
////////////////////////////////////////////////////
uint Matrix::GetColumns() const
{
	return m_cols;
}

////////////////////////////////////////////////////
/// \brief Returns wether you could invert the matrix or not
///
////////////////////////////////////////////////////
bool Matrix::IsInvertable()
{
	return det(*this) == 0 ? false : true;
}

////////////////////////////////////////////////////
/// \brief Returns wether the columns and rows match
///
////////////////////////////////////////////////////
bool Matrix::IsSquare()
{
	return m_rows == m_cols;
}


////////////////////////////////////////////////////
/// \brief Prints the whole matrix to the console
///
////////////////////////////////////////////////////
inline void Matrix::PrintMatrix()
{
	for (uint row = 0; row < m_rows; row++)
	{
		for (uint col = 0; col < m_cols; col++)
		{
			std::cout << m_matrix[row][col];
			std::cout << " ";
		}

		std::cout << std::endl;
	}
}


inline bool Delta(uint right, uint left)
{
	return (right == left);
}


////////////////////////////////////////////////////
/// \brief Transposes the matrix
///
////////////////////////////////////////////////////
inline void Matrix::Transpose()
{
	doubleMatrix transposed;
	Resize(m_cols, m_rows, transposed);

	for (int row = 0; row < m_rows; row++)
		for (int col = 0; col < m_cols; col++)
			transposed[col][row] = m_matrix[row][col];
	
	m_matrix = transposed;

	std::swap(m_rows, m_cols);
}


////////////////////////////////////////////////////
/// \brief Multiplies each value of a given row by a given factor
///
/// \param row The row that the multiplication should be applied to
/// \param factor The factor to multiply the row by
///
////////////////////////////////////////////////////
inline void Matrix::MultiplyRow(uint row, double factor)
{
	CHECK_BOUNDS(row, 0) _VOID;

	for (int col = 0; col < m_cols; col++)
		m_matrix[row][col] *= factor;
}


////////////////////////////////////////////////////
/// \brief Swaps two rows in the matrix
///
/// \param left One row
/// \param right The other row
///
////////////////////////////////////////////////////
inline void Matrix::SwapRows(uint left, uint right)
{
	CHECK_BOUNDS(left, right) _VOID;

	std::swap(m_matrix[left], m_matrix[right]);
}


////////////////////////////////////////////////////
/// \brief Swaps two rows in the matrix
///
/// \param left One row
/// \param right The other row
///
////////////////////////////////////////////////////
inline void Matrix::AddMultiplesToRow(uint base, uint target, double factor)
{
	int djfids = 1;

	CHECK_BOUNDS(base, 0) _VOID;
	CHECK_BOUNDS(target, 0) _VOID;

	for(int col = 0; col < m_cols; col++)
		m_matrix[target][col] += m_matrix[base][col] * factor;
}


inline double det(const Matrix& mat)
{
	Matrix holder = mat;
	holder.TransformToEchelonForm();

	double det = 1;

	for (int col = 0; col < mat.m_cols; col++) {
		for (int row = col; row < mat.m_rows; row++) {
			if (holder.m_matrix[row][col] != 0) {
				det *= holder.m_matrix[row][col];
				break;
			}
		}
	}

	return det;
}


////////////////////////////////////////////////////
/// \brief Performs matrix multiplication
///
/// \param Two matrices to multiply
///
/// \return A multiplied matrix when success. An empty matrix if failed
///
////////////////////////////////////////////////////
inline Matrix mul(const Matrix& right, const Matrix& left)
{
	if (right.m_cols != left.m_rows)
		return Matrix(0, 0);

	Matrix returnMatrix = Matrix(right.m_rows, left.m_cols);

	for (uint row = 0; row < returnMatrix.m_rows; row++)
		for (uint col = 0; col < returnMatrix.m_cols; col++)
			for (uint k = 0; k < right.m_cols; k++)
				returnMatrix.m_matrix[row][col] += right.m_matrix[row][k] * left.m_matrix[k][col];

	return returnMatrix;
}


////////////////////////////////////////////////////
/// \brief Transforms a matrix to reduced row echelon form py performing Gaussian elemination
///
////////////////////////////////////////////////////
inline void Matrix::TransformToEchelonForm()
{
	// Quick description of the algorithm used:
	//
	// When using the Gaussian algorithm you start in the first column.
	// The first number of the upper row must be a non-zero value, so if it happens to be a zero,
	// swap this row with any other row starting with a non-zero integer.
	// After that you must make any values BELOW the upper row equal zero, by adding fitting multiples
	// of the upper row to each row. As soon as the upper row is the only row in the column to have a non-zero
	// value, continue on to the next colum, but "cross out" the upper row (the second row now becomes the upper row)
	// Repeat until every column has been transformed.
	//
	// After that, the matrix should be in "echelon form", which means that the matrix looks something like this:
	//
	// [ 3 4 9 2 8 ]
	// [ 0 4 1 8 2 ]
	// [ 0 0 8 3 5 ]
	// [ 0 0 0 0 9 ]
	// [ 0 0 0 0 0 ]
	//
	// PROS to using this algorithm:
	// - It's the only one I know (yet)
	//
	// CONS to using this algorithm:
	// - Probably slow and hideous

	int currentRow = 0;

	// The algorithm will move through each column
	for (int col = 0; col < m_cols; col++)
	{
		if (currentRow > m_rows)
			break;

		// Check if some value in the first column is zero and count how many
		int nonZeroRow = -1;
		bool allZeros = true;
		for (int row = currentRow; row < m_rows; row++)
		{
			// If one is found, set rowWithZero to it and zeros++
			if (m_matrix[row][col] != 0) 
			{
				nonZeroRow = row;
				allZeros = false;
			}
		}
		
		// If all rows are zeros, then skip
		if(allZeros)
		{
			continue;
		}

		// If no row is zero, make all except the current one to zero!
		else if (nonZeroRow < 0)
		{
			for (int row = currentRow + 1; row < m_rows; row++)
			{
				AddMultiplesToRow(currentRow, row, -(m_matrix[row][col] / m_matrix[currentRow][col]));
			}
		}

		// If current row is zero, swap a nonZero row with it, then make them all zero
		else if(nonZeroRow != currentRow)
		{
			// Swap rows
			std::swap(m_matrix[nonZeroRow], m_matrix[currentRow]);
			
			// Make all rows non-zero
			for (int row = currentRow + 1; row < m_rows; row++)
			{
				if (m_matrix[row][col] == 0) {

					continue;
				}

				AddMultiplesToRow(currentRow, row, -(m_matrix[row][col] / m_matrix[currentRow][col]));
			}
		}

		for (int row = currentRow + 1; row < m_rows; row++)
		{
			if (m_matrix[row][col] == 0) {
				
				continue;
			}

			AddMultiplesToRow(currentRow, row, -(m_matrix[row][col] / m_matrix[currentRow][col]));
		}

		currentRow++;
	}
}

inline void Matrix::TransformToReducedEchelonForm()
{
	TransformToEchelonForm();



	// The matrix is now in echelon form. It would be better if it was in reduced echelon form.
	//
	// This means that all the first numbers in each row must be equal to one, and all first
	// values in each row must be the only non-zero value in their column.
	// It looks like this:
	//
	// [ 1 0 0 2 0 ]
	// [ 0 1 0 8 0 ]
	// [ 0 0 1 3 0 ]
	// [ 0 0 0 0 1 ]
	// [ 0 0 0 0 0 ]

	// Starting from the right, add fitting multiples of the first non-zero row to all the ones above
	for (int col = m_cols - 1; col >= 0; col--)
	{
		// Find first non-zero number (starting from the bottom)
		int nonZeroRow = m_rows - 1;
		for (int row = m_rows - 1; row >= 0; row--)
		{
			if (m_matrix[row][col] != 0)
			{
				nonZeroRow = row;
				break;
			}
		}

		// Divide the bottom-most row by itself to make it = 1
		MultiplyRow(nonZeroRow, 1 / m_matrix[nonZeroRow][col]);

		// Add fitting multiples of that row to all the others
		for (int row = nonZeroRow - 1; row >= 0; row--)
		{
			AddMultiplesToRow(nonZeroRow, row, -(m_matrix[row][col] / m_matrix[nonZeroRow][col]));
		}
	}

	// Convert all the "-0" to "0" (idk why c++ does this)
	for (int row = 0; row < m_rows; row++)
		for (int col = 0; col < m_cols; col++)
			if (m_matrix[row][col] == -0) m_matrix[row][col] = 0;
}


inline bool Matrix::InBounds(uint row, uint col) const
{
	return ((row < m_rows) && (col < m_cols)) ? true : false;
}

inline bool Matrix::DimensionsFitting(const Matrix& right, const Matrix& left)
{
	return (right.GetRows() == left.GetRows()) && (right.GetColumns() == left.GetColumns());
}

inline void Matrix::Resize(uint rows, uint cols, doubleMatrix& matrix)
{
	rows = rows == 0 ? 1 : rows;
	cols = cols == 0 ? 1 : cols;

	matrix.resize(rows);
	for (int row = 0; row < rows; row++)
		matrix[row].resize(cols);
}




////////////////////// OPERATOR OVERLOADING //////////////////////////

Matrix Matrix::operator+(const Matrix& other)
{
	if (!DimensionsFitting(*this, other))
		return Matrix(0, 0);

	Matrix returnMatrix(m_rows, m_cols);

	for (int row = 0; row < m_rows; row++)
		for (int col = 0; col < m_cols; col++)
			returnMatrix.SetNumber(row, col, m_matrix[row][col] + other.GetNumber(row, col));

	return returnMatrix;
}


Matrix Matrix::operator-(const Matrix& other)
{
	if (!DimensionsFitting(*this, other))
		return Matrix(0, 0);

	Matrix returnMatrix(m_rows, m_cols);

	for (int row = 0; row < m_rows; row++)
		for (int col = 0; col < m_cols; col++)
			returnMatrix.SetNumber(row, col, m_matrix[row][col] - other.GetNumber(row, col));

	return returnMatrix;
}


Matrix Matrix::operator*(const double& other)
{
	Matrix returnMatrix(m_rows, m_cols);

	for (int row = 0; row < m_rows; row++)
		for (int col = 0; col < m_cols; col++)
			returnMatrix.SetNumber(row, col, m_matrix[row][col] * other);

	return returnMatrix;
}


Matrix Matrix::operator/(const double& other)
{
	return (*this * (1 / other));
}


void Matrix::operator+=(const Matrix& other)
{
	*this = *this + other;
}

void Matrix::operator-=(const Matrix& other)
{
	*this = *this - other;
}

void Matrix::operator*=(const double& other)
{
	*this = *this * other;
}

void Matrix::operator/=(const double& other)
{
	*this = *this / other;
}


Matrix Matrix::operator*(const Matrix& other)
{
	return mul(*this, other);
}

void Matrix::operator*=(const Matrix& other)
{
	*this = mul(*this, other);
}


std::ostream& operator<<(std::ostream& stream, const Matrix& other)
{
	for (uint row = 0; row < other.m_rows; row++)
	{
		for (uint col = 0; col < other.m_cols; col++)
		{
			stream << other.m_matrix[row][col] << " ";
		}

		stream << std::endl;
	}

	return stream;
}
#endif // MATRIX_h