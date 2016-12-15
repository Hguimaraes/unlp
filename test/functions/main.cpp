#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

namespace unlp {

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

/**
 * @Description: 
 *  Abstract Class to represent the mathematical functions of this assignment
 */
class Functions {
  public:
    /**
     * Return the value of the function in the given points
     */
    virtual double fun(Matrix points) = 0;
    
    /**
     * Return a Vector of the gradient in the given points
     */
    virtual Matrix gradient(Matrix points) = 0;
    
    /**
     * Return the Matrix representing the Hessian of the given points
     */
    virtual Matrix hessian(Matrix points) = 0;
};

/**
 * @Description: 
 *  Class to represent the mathematical function (A) of this assignment
 */
class NLAFunction : public Functions {
  public:
    double fun(Matrix points);
    Matrix gradient(Matrix points);
    Matrix hessian(Matrix points);
};

double NLAFunction::fun(Matrix points){
  double a = points.at(1);
  double b = points.at(2);
  return log(a*b*(1-a)*(1-b));
}

Matrix NLAFunction::gradient(Matrix points){
  vector<double> tmp;
  double a = points.at(1);
  double b = points.at(2);
  
  double grad_a = (1 - 2*a)/(a - pow(a,2));
  double grad_b = (1 - 2*b)/(b - pow(b,2));

  tmp.push_back(grad_a);
  tmp.push_back(grad_b);

  Matrix w(tmp);
  return w;
}

Matrix NLAFunction::hessian(Matrix points){
  double a = points.at(1);
  double b = points.at(2);

  vector<vector<double> > tmp;
  vector<double> grad_a;
  vector<double> grad_b;

  grad_a.push_back((-2*pow(a,2)+2*a-1)/(pow(1-a,2)*pow(a,2)));
  grad_a.push_back(0);

  grad_b.push_back(0);
  grad_b.push_back((-2*pow(b,2)+2*b-1)/(pow(1-b,2)*pow(b,2)));

  tmp.push_back(grad_a);
  tmp.push_back(grad_b);

  Matrix w(tmp);
  return w;
}

/**
 * @Description: 
 *  Class to represent the mathematical function (B) of this assignment
 */
class NLBFunction : public Functions {
  public:
    double fun(Matrix points);
    Matrix gradient(Matrix points);
    Matrix hessian(Matrix points);
};

double NLBFunction::fun(Matrix points){
  double a = points.at(1);
  double b = points.at(2);
  return log(1 + pow(log(a),2) + pow(log(b),2));
}

Matrix NLBFunction::gradient(Matrix points){
  vector<double> tmp;
  double a = points.at(1);
  double b = points.at(2);
  
  double grad_a = (2*log(a)) / (a*fun(points));
  double grad_b = (2*log(b)) / (b*fun(points));

  tmp.push_back(grad_a);
  tmp.push_back(grad_b);

  Matrix w(tmp);
  return w;
}

Matrix NLBFunction::hessian(Matrix points){
  double a = points.at(1);
  double b = points.at(2);

  vector<vector<double> > tmp;
  vector<double> grad_a;
  vector<double> grad_b;

  double mixedDerivatives = -(4*log(a)*log(b))/(a*b*pow(fun(points),2));

  grad_a.push_back(
    -(2*(log(a)*(pow(log(b),2) + 1) + pow(log(a),3) + pow(log(a),2) - pow(log(b),2) - 1))
    / (pow(a,2)*pow(fun(points),2))
  );
  
  grad_a.push_back(mixedDerivatives);

  grad_b.push_back(mixedDerivatives);
  
  grad_b.push_back(
     -(2*(pow(log(a),2)*(log(b) - 1) + pow(log(b),3) + pow(log(b),2) + log(b) - 1))
    /(pow(b,2)*pow(fun(points),2))
  );

  tmp.push_back(grad_a);
  tmp.push_back(grad_b);

  Matrix w(tmp);
  return w;
}

}

int main(int argc, char **argv){
  vector<double> tmp;
  tmp.push_back(0.8);
  tmp.push_back(0.1);

  unlp::Matrix point(tmp);
  unlp::NLAFunction fa;

  cout << fa.fun(point) << endl << endl;

  unlp::Matrix grad = fa.gradient(point);
  unlp::Matrix hess = fa.hessian(point);
  grad.printMatrix();
  hess.printMatrix();

  cout << endl << endl << endl;

  unlp::NLBFunction fb;

  cout << fb.fun(point) << endl << endl;

  unlp::Matrix gradB = fb.gradient(point);
  unlp::Matrix hessB = fb.hessian(point);
  gradB.printMatrix();
  hessB.printMatrix();

  return 0;
}