#include <stdlib.h>
#include "SML.h"

LCRS1 *GET_LCRS1(long dim, long max) {
   
   LCRS1 *Matrix = malloc(sizeof(LCRS1));
   
   Matrix->max_val = max;
   Matrix->max_row = dim+1;
   
   Matrix->row_dim = 0;
   Matrix->col_dim = 0;
   Matrix->Val = GET_ARRAY_DOUBLE1(max);
   Matrix->Col = GET_ARRAY_LINT1(max);
   Matrix->Row = GET_ARRAY_LINT1(dim+1);
   
   return Matrix;
   
}
