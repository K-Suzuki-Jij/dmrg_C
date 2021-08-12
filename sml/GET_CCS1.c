#include <stdlib.h>
#include "SML.h"

CCS1 *GET_CCS1(int dim, long max) {
   
   CCS1 *Matrix = malloc(sizeof(CCS1));
   
   Matrix->max_val = max;
   Matrix->max_col = dim + 1;
   
   Matrix->row_dim = 0;
   Matrix->col_dim = 0;
   Matrix->Val = GET_ARRAY_DOUBLE1(max);
   Matrix->Col = GET_ARRAY_LINT1(dim + 1);
   Matrix->Row = GET_ARRAY_INT1(max);
   
   return Matrix;
   
}
