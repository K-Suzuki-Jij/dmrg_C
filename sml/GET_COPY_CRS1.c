#include <stdlib.h>
#include "SML.h"

CRS1 *GET_COPY_CRS1(CRS1 *M) {
   
   int  dim      = M->row_dim;
   long elem_num = M->Row[dim];
   
   CRS1 *Matrix = malloc(sizeof(CRS1));
   
   Matrix->max_val = elem_num;
   Matrix->max_row = dim+1;
   
   Matrix->row_dim = dim;
   Matrix->col_dim = dim;
   Matrix->Val = GET_ARRAY_DOUBLE1(elem_num);
   Matrix->Col = GET_ARRAY_INT1(elem_num);
   Matrix->Row = GET_ARRAY_LINT1(dim+1);
   
   COPY_CRS1(M, Matrix, 1);
   
   return Matrix;
   
}
