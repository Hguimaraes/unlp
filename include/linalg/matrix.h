#ifndef MATRIX_H_
#define MATRIX_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

namespace unlp{

/**
 * @Description: 
 *  Class to represent the mathematical definition of matrix.
 * It is intended to be used with 2x2 matrix in this assignment,
 * but can be used with  m x n matrix (The determinant will not be available).
 */
class Matrix {
  public:
    /**
     *  Construct a matrix object with the specified dimensions (rows x cols)
     * initialized to value.
     */
    Matrix(unsigned m, unsigned n, double value = 0.0);
    
    // Construct a column vector from a vector.
    Matrix(const vector<double>&);

    // Construct a matrix from a vector of vectors.
    Matrix(const vector<vector<double> >&);
        
    // Set Aij to value.
    void set(unsigned i, unsigned j, double value);

    // Usual matrix operations.
    Matrix operator+(const Matrix&) const;
    Matrix operator-(const Matrix&) const;
    Matrix operator*(const Matrix&) const;
    Matrix operator*(double) const;
    Matrix operator/(double) const;

    // Return the transpose of the matrix
    Matrix t() const;
    
    // Determinant of a 2 x 2 Matrix
    double det2() const;
    
    /**
     * If is a column or row matrix, return the ith element
     */
    double at(unsigned i) const;
    
    /**
     * Return the Aij element of the matrix
     */
    double at(unsigned i, unsigned j) const;

    /**
     * Debug/Auxiliar methods
     */
    void printMatrix() const;
    unsigned getCols() const;
    unsigned getRows() const;

  private:
    /** 
     * m = number of lines
     * n = number of columns
     * v = internal representation of the matrix
     */
    unsigned m, n;
    vector< vector<double> > v;

};

}

#endif