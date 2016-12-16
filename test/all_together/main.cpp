#include <iostream>
#include <vector>
#include <cmath>

#define VERBOSE true
#define NI 0.25
#define GAMMA 0.8

using namespace std;

namespace unlp {

/**
 * +---------------------------+
 * |          Matrix           |
 * +---------------------------+
 */

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
    
    /// Copy a existing matrix.
    Matrix(const Matrix&);

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
    Matrix operator+(double) const;
    Matrix operator-(double) const;
    Matrix operator*(double) const;
    Matrix operator/(double) const;
    friend Matrix operator+(double, const Matrix&);
    friend Matrix operator-(double, const Matrix&);
    friend Matrix operator*(double, const Matrix&);
    friend Matrix operator/(double, const Matrix&);

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

Matrix::Matrix(const Matrix& o) :
  m(o.m),
  n(o.n),
  v(o.v) {
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

Matrix Matrix::operator+(double k) const {
  Matrix a(m, n);
  for (unsigned i = 1; i <= m; ++i)
    for (unsigned j = 1; j <= n; ++j)
      a.set(i, j, at(i,j) + k);
  return a;
}

Matrix Matrix::operator-(double k) const {
  Matrix a(m, n);
  for (unsigned i = 1; i <= m; ++i)
    for (unsigned j = 1; j <= n; ++j)
      a.set(i, j, at(i,j) - k);
  return a;
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

Matrix operator+(double k, const Matrix& o) {
  return o + k;
}

Matrix operator-(double k, const Matrix& o) {
  return o - k;
}

Matrix operator*(double s, const Matrix& o) {
  return o * s;
}

Matrix operator/(double s, const Matrix& o) {
  return o / s;
}

Matrix inv2(Matrix x){
  
  if(x.getCols() != 2 || x.getRows() != 2)
    throw std::invalid_argument("ERROR: The Matrix should have the size 2x2");

  Matrix tmp = x;
  tmp.set(1,1,x.at(2,2));
  tmp.set(2,1,-x.at(2,1));
  tmp.set(1,2,-x.at(1,2));
  tmp.set(2,2,x.at(1,1));
  tmp = tmp/x.det2();
  return tmp;
}

Matrix identity(unsigned dim){
  Matrix id(dim, dim, 0.0);
    
  // Fill the identity matrix
  for (int i = 1; i <= dim; ++i)
    id.set(i,i, 1.0);

  return id;
}

/**
 * +---------------------------+
 * |         Functions         |
 * +---------------------------+
 */


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
  return -log(a*b*(1-a)*(1-b));
}

Matrix NLAFunction::gradient(Matrix points){
  vector<double> tmp;
  double a = points.at(1);
  double b = points.at(2);
  
  double grad_a = -(1 - 2*a)/(a - pow(a,2));
  double grad_b = -(1 - 2*b)/(b - pow(b,2));

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

  grad_a.push_back(-(-2*pow(a,2)+2*a-1)/(pow(1-a,2)*pow(a,2)));
  grad_a.push_back(0);

  grad_b.push_back(0);
  grad_b.push_back(-(-2*pow(b,2)+2*b-1)/(pow(1-b,2)*pow(b,2)));

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
  double grad_a = (2*log(a)) / (a*(1 + pow(log(a),2) + pow(log(b),2)));
  double grad_b = (2*log(b)) / (b*(1 + pow(log(a),2) + pow(log(b),2)));

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

  double mixedDerivatives = -(4*log(a)*log(b))/(a*b*pow((1 + pow(log(a),2) + pow(log(b),2)),2));

  grad_a.push_back(
    -(2*(log(a)*(pow(log(b),2) + 1) + pow(log(a),3) + pow(log(a),2) - pow(log(b),2) - 1))/(pow(a,2)*pow((1 + pow(log(a),2) + pow(log(b),2)),2))
  );
  
  grad_a.push_back(mixedDerivatives);

  grad_b.push_back(mixedDerivatives);
  
  grad_b.push_back(
     -(2*(pow(log(a),2)*(log(b) - 1) + pow(log(b),3) + pow(log(b),2) + log(b) - 1))/(pow(b,2)*pow((1 + pow(log(a),2) + pow(log(b),2)),2))
  );

  tmp.push_back(grad_a);
  tmp.push_back(grad_b);

  Matrix w(tmp);
  return w;
}

/**
 * @Description: 
 *  Class to represent a mathematical function (Example)
 */
class NLExFunction : public Functions {
  public:
    double fun(Matrix points);
    Matrix gradient(Matrix points);
    Matrix hessian(Matrix points);
};

double NLExFunction::fun(Matrix points){
  double a = points.at(1);
  double b = points.at(2);
  return pow(a,4) + pow(b,4);
}

Matrix NLExFunction::gradient(Matrix points){
  vector<double> tmp;
  double a = points.at(1);
  double b = points.at(2);
  
  double grad_a = 4*pow(a,3);
  double grad_b = 4*pow(b,3);

  tmp.push_back(grad_a);
  tmp.push_back(grad_b);

  Matrix w(tmp);
  return w;
}

Matrix NLExFunction::hessian(Matrix points){
  double a = points.at(1);
  double b = points.at(2);

  vector<vector<double> > tmp;
  vector<double> grad_a;
  vector<double> grad_b;

  double mixedDerivatives = 0;

  grad_a.push_back(12*pow(a,2));
  
  grad_a.push_back(mixedDerivatives);

  grad_b.push_back(mixedDerivatives);
  
  grad_b.push_back(12*pow(b,2));

  tmp.push_back(grad_a);
  tmp.push_back(grad_b);

  Matrix w(tmp);
  return w;
}

/**
 * +---------------------------+
 * |        Algorithms         |
 * +---------------------------+
 */

bool isClose(Matrix x, Matrix y, double error_max){
  if(x.getRows() != y.getRows() || x.getCols() != y.getCols())
    throw std::invalid_argument("ERROR: The arrays should have the same size");
  
  double sum_distance = 0;
  for (int i = 1; i <= x.getRows(); ++i){
    sum_distance += pow(abs(x.at(i) - y.at(i)),2);
  }
  
  if(sqrt(sum_distance) <= error_max) return true;
  return false;

}

double armijo(Functions& func, Matrix x,
  Matrix d, double n, double gamma, unsigned maxIter){
  double step = 1;
  bool wolfe_condition = false;
  int it = 0;

  while(!wolfe_condition){
    double leftCond = func.fun(x + step*d);
    double rightCond = func.fun(x) + n*step*(func.gradient(x).t()*d).at(1);
    wolfe_condition = leftCond <= rightCond;

    if(!wolfe_condition)
      step *= gamma;

    it++;
    if(it == maxIter) break;
  }

  return step;
}

Matrix BFGS(Functions& func, Matrix x, Matrix x_old, Matrix hessian){
  Matrix p = x - x_old;
  Matrix q = func.gradient(x) - func.gradient(x_old);

  double denominator = (p.t()*q).at(1);
  if(denominator != 0){
    Matrix tmp_a = ((q.t()*hessian)*q)/denominator;
    Matrix tmp_b = (p.t()*p)/denominator;
    Matrix tmp_c = p*q.t()*hessian;
    Matrix tmp_d = hessian*q*p.t();
    Matrix tmp_e = tmp_c + tmp_d;
    
    hessian = hessian + ((1.0 + tmp_a.at(1))*tmp_b).at(1) - (tmp_e/denominator);
  } else {
    hessian = Matrix(2,2,0);
  }

  return hessian;

}

void gradient_descent(Functions& func, Matrix init_point,
  double error_max = 1e-7, unsigned maxIter = 100000)
{
  // Initial parameters
  Matrix x = init_point;
  Matrix prev = x + 1000;
  int it = 0;
  double step = 0.00001;

  while(!isClose(x, prev, error_max) && it != maxIter){
    Matrix direction = -1*func.gradient(x);
    step = armijo(func, x, direction, NI, GAMMA, maxIter);
    prev = x;
    x = x + step*direction;
    it++;
  }

  cout << "number of it = " << it << endl;
  cout << "optimal value = " << func.fun(x) << endl;

  // Print Report
  x.printMatrix();
}

void newton(unlp::Functions& func, unlp::Matrix init_point,
  double error_max = 1e-5, unsigned maxIter = 100000)
{
  // Starting the algorithm in the initial points, setting the previous
  // value to Null and starting at iteration 0
  Matrix x = init_point;
  Matrix prev = x + 1000;
  int it = 0;
  double step = 0.00001;

  while(!isClose(x, prev, error_max) && it != maxIter){
    Matrix direction = -1*inv2(func.hessian(x))*func.gradient(x);
    step = armijo(func, x, direction, NI, GAMMA, maxIter);
    prev = x;
    x = x + step*direction;
    it++;
  }

  cout << "number of it = " << it << endl;
  cout << "optimal value = " << func.fun(x) << endl;

  // Print Report
  x.printMatrix();
}

void quasi_newton(unlp::Functions& func, unlp::Matrix init_point,
  double error_max = 1e-5, unsigned maxIter = 100000)
{
  // Starting the algorithm in the initial points, setting the previous
  // value to Null and starting at iteration 0
  Matrix x = init_point;
  Matrix prev = x + 1000;
  int it = 0;
  double step = 0.00001;
  Matrix hessian = identity(2);

  while(!isClose(x, prev, error_max) && it != maxIter){
    Matrix direction = -1*hessian*func.gradient(x);
    step = armijo(func, x, direction, NI, GAMMA, maxIter);
    prev = x;
    x = x + step*direction;
    it++;

    hessian = BFGS(func, x, prev, hessian);

  }

  cout << "number of it = " << it << endl;
  cout << "optimal value = " << func.fun(x) << endl;

  // Print Report
  x.printMatrix();
}

}

int main(int argc, char **argv){
  vector<double> tmp;
  tmp.push_back(1);
  tmp.push_back(0.8);

  unlp::Matrix point(tmp);
  
  /*// fA
  unlp::NLAFunction fa;

  cout << "GRADIENT METHOD: " << endl;
  unlp::gradient_descent(fa, point);

  cout << endl;
  cout << "NEWTON METHOD: " << endl;
  unlp::newton(fa, point);

  cout << endl;
  cout << "QUASI-NEWTON METHOD: " << endl;
  unlp::quasi_newton(fa, point);

  cout << endl;
  cout << endl;
  cout << endl;*/
  
  // fB
  unlp::NLBFunction fb;

  cout << "GRADIENT METHOD: " << endl;
  unlp::gradient_descent(fb, point);

  cout << endl;
  cout << "NEWTON METHOD: " << endl;
  unlp::newton(fb, point);

  cout << endl;
  cout << "QUASI-NEWTON METHOD: " << endl;
  unlp::quasi_newton(fb, point);
  return 0;
}