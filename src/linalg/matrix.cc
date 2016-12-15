#include "lingalg/matrix.h"

namespace unlp {

Matrix::Matrix(unsigned m, unsigned n, double value): m(m), n(n) {
  v = vector< vector<double> >(m, vector<double>(n, value));
}

Matrix::Matrix(const vector<double>& o): m(o.size()), n(1) {
  for (int i = 0; i < o.size(); ++i)
    v.push_back(vector<double>(1, o[i]));
}

Matrix::Matrix(const vector<vector<double> >& o): 
  m(o.size()), n(o[0].size()), v(o)
{
}

void Matrix::set(unsigned i, unsigned j, double value) {
  v[i-1][j-1] = value;
}

Matrix Matrix::t() const {
  Matrix w(getCols(), getRows());
  for (unsigned i = 1; i <= getRows(); ++i)
    for (unsigned j = 1; j <= getCols(); ++j)
      w.set(j,i,at(i,j));
  return w;
}

double Matrix::det2() const {
  if (m != 2 && n != 2)
    throw std::invalid_argument("ERROR: Can't apply det2 to a non 2x2 matrix");
  return at(1,1) * at(2,2) - at(1,2) * at(2,1);
}

double Matrix::at(unsigned i) const {
 if(n == 1) return v[i-1][0];
 else
  throw std::invalid_argument("ERROR: This is not a row or column vector");
}

double Matrix::at(unsigned i, unsigned j) const {
  return v[i-1][j-1];
}

void Matrix::printMatrix() const {
  std::cout << "INFO: Matrix debug" << std::endl;
  std::cout << "\t" << "#rows=" << m << ", #cols=" <<  n << std::endl;
  for (unsigned i = 1; i <= m; ++i) {
    std::cout << "\t";
    for (unsigned j = 1; j <= n; ++j)
      std::cout << at(i,j) << " ";
    std::cout << std::endl;
  }
}

unsigned Matrix::getCols() const {
  return n;
}

unsigned Matrix::getRows() const {
  return m;
}

Matrix Matrix::operator+(const Matrix& o) const {
  Matrix a(m, n);
  for (unsigned i = 1; i <= m; ++i)
    for (unsigned j = 1; j <= n; ++j)
      a.set(i, j, at(i,j) + o.at(i,j));
  return a;
}

Matrix Matrix::operator-(const Matrix& o) const {
  Matrix a(m, n);
  for (unsigned i = 1; i <= m; ++i)
    for (unsigned j = 1; j <= n; ++j)
      a.set(i, j, at(i,j) - o.at(i,j));
  return a;
}

Matrix Matrix::operator*(const Matrix& o) const {
  Matrix w(getRows(), o.getCols());
  
  // Check the dimensions
  if (getCols() != o.getRows())
    throw std::invalid_argument("ERROR: Invalid matrix multiplication");
  
  for (unsigned i = 1; i <= getRows(); ++i){
    for (unsigned j = 1; j <= o.getCols(); ++j) {
      double sum = 0.0;
      for (unsigned k = 1; k <= getCols(); ++k) {
        sum += at(i,k) * o.at(k,j);
      }
      w.set(i,j,sum);
    }
  }

  return w;
}

Matrix Matrix::operator*(double k) const {
  Matrix a(m, n);
  for (unsigned i = 1; i <= m; ++i)
    for (unsigned j = 1; j <= n; ++j)
      a.set(i, j, at(i,j) * k);
  return a;
}

Matrix Matrix::operator/(double k) const {
  Matrix a(m, n);
  for (unsigned i = 1; i <= m; ++i)
    for (unsigned j = 1; j <= n; ++j)
      a.set(i, j, at(i,j) / k);
  return a;
}

}