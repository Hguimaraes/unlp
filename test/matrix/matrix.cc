#include "matrix.h"

namespace unlp {

// OK
Matrix::Matrix(unsigned m, unsigned n, double value): m(m), n(n) {
  v = vector< vector<double> >(m, vector<double>(n, value));
}

// OK
Matrix::Matrix(const vector<double>& o): m(o.size()), n(1) {
  for (int i = 0; i < o.size(); ++i)
    v.push_back(vector<double>(1, o[i]));
}

// OK
Matrix::Matrix(const vector<vector<double> >& o): 
  m(o.size()), n(o[0].size()), v(o)
{
}

// OK
void Matrix::set(unsigned i, unsigned j, double value) {
  v[i-1][j-1] = value;
}


// OK
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
 if(n == 1) return v[i][0];
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
  return m;
}

unsigned Matrix::getRows() const {
  return n;
}
}