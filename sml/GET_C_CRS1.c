#include <stdlib.h>
#include <complex.h>
#include "SML.h"

C_CRS1 *GET_C_CRS1(int dim, long max) {
   
   C_CRS1 *Matrix = malloc(sizeof(C_CRS1));
   
   Matrix->max_val = max;
   Matrix->max_row = dim+1;
   
   Matrix->row_dim = 0;
   Matrix->col_dim = 0;
   Matrix->Val = GET_ARRAY_C_DOUBLE1(max);
   Matrix->Col = GET_ARRAY_INT1(max);
   Matrix->Row = GET_ARRAY_LINT1(dim+1);
   
   return Matrix;
   
}
