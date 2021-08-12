#include <stdlib.h>
#include "SML.h"

void FREE_CRS4(CRS1 ****Matrix, int row1, int row2, int row3) {
   
   int i,j,k;
   
   for (i = 0; i < row1; i++) {
      for (j = 0; j < row2; j++) {
         for( k = 0; k < row3; k++) {
            FREE_CRS1(Matrix[i][j][k]);
         }
      }
   }
   
   for (i = 0; i < row1; i++) {
      for (j = 0; j < row2; j++) {
         free(Matrix[i][j]);
      }
   }
   
   for (i = 0; i < row1; i++) {
      free(Matrix[i]);
   }
   
   free(Matrix);
   
}
