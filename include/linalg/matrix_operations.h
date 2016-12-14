#ifndef MATRIX_OPERATIONS_H_
#define MATRIX_OPERATIONS_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "matrix.h"

using namespace std;

namespace unlp{

/**
 * @Description: 
 *  Collection of methods useful to handle matrix operations for this
 * homework.
 */

unlp::Matrix identity(unsigned dim);
unlp::Matrix inv(unlp::Matrix x);
double dot(unlp::Matrix x, unlp::Matrix y);
}

#endif