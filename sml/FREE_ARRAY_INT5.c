#include <stdlib.h>

void FREE_ARRAY_INT5(int *****Matrix, long row, long col, long col_2, long col_3) {
   
   long i,j,k,l;
   
   for (i = 0; i < row; i++) {
      for (j = 0; j < col; j++) {
         for (k = 0; k < col_2; k++) {
            for (l = 0; l < col_3; l++) {
               free(Matrix[i][j][k][l]);
            }
         }
      }
   }
   
   for (i = 0; i < row; i++) {
      for (j = 0; j < col; j++) {
         for (k = 0; k < col_2; k++) {
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
