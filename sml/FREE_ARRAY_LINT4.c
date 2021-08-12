#include <stdlib.h>

void FREE_ARRAY_LINT4(long ****Matrix, long row, long col, long col2) {
   
   long i,j,k;
   
   for (i = 0; i < row; i++) {
      for (j = 0; j < col; j++) {
         for (k = 0; k < col2; k++) {
            free(Matrix[i][j][k]);
         }
      }
   }
   
   for (i = 0; i < row; i++) {
      for (j = 0; j < col; j++) {
         free(Matrix[i][j]);
      }
   }
   
   for (i = 0; i < row; i++) {
      free(Matrix[i]);
   }
   
   free(Matrix);
}
