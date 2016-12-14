#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "linalg/matrix.h"

using namespace std;

namespace unlp{

/**
 * @Description: 
 *  Abstract Class to represent the mathematical functions of this assignment
 */
class Functions {
  /**
   * Return the second derivative
   */
  virtual long double secondOrderDerivative(Matrix points) = 0;
  
  /**
   * Return the second order mixed derivative of our function at the "point"
   */
  virtual double secondOrderMixedDerivative(Matrix points) = 0;

public:
  /**
   * Return the value of the function in the given points
   */
  virtual double function(Matrix points) = 0;
  
  /**
   * Return a Vector of the gradient in the given points
   */
  virtual Matrix gradient(Matrix points) = 0;
  
  /**
   * Return the Matrix representing the Hessian of the given points
   */
  virtual Matrix hessian(Matrix points) = 0;
};

}

#endif