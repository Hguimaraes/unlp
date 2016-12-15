#ifndef NLB_FUNCTIONS_H_
#define NLB_FUNCTIONS_H_

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
 *  Class to represent the mathematical function (B) of this assignment
 */
class NLBFunction : public Functions {
  long double secondOrderDerivative(Matrix points);
  double secondOrderMixedDerivative(Matrix points);
  
  public:
    double function(Matrix points);
    Matrix gradient(Matrix points);
    Matrix hessian(Matrix points);
};

}

#endif