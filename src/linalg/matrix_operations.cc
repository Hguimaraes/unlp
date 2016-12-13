#include "matrix_operations.h"

namespace unlp {

Matrix identity(unsigned dim){
  Matrix id(0.0);
    
  // Fill the identity matrix
  for (int i = 0; i < dim; ++i)
    id.set(i,i, 1.0);

  return id;
}

}