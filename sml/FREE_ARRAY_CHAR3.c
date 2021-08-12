#include <stdlib.h>

void FREE_ARRAY_CHAR3(char ***Matrix, long row, long col) {
   
   long i,j;
   
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
