#include <stdlib.h>
#include "SML.h"

void FREE_CRS3(CRS1 ***Matrix, int row1, int row2) {
   
   int i,j;
   
   for (i = 0; i < row1; i++) {
      for (j = 0; j < row2; j++) {
         FREE_CRS1(Matrix[i][j]);
      }
   }
   
   for (i = 0; i < row1; i++) {
      free(Matrix[i]);
   }
   
   free(Matrix);
   
}
