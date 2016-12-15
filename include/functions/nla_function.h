#ifndef NLA_FUNCTIONS_H_
#define NLA_FUNCTIONS_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "functions.h"
#include "linalg/matrix.h"

using namespace std;

namespace unlp{

/**
 * @Description: 
 *  Class to represent the mathematical function (A) of this assignment
 */
class NLAFunction : public Functions {
  long double secondOrderDerivative(Matrix points);
  double secondOrderMixedDerivative(Matrix points);
  
  public:
    double function(Matrix points);
    Matrix gradient(Matrix points);
    Matrix hessian(Matrix points);
};

}

#endif