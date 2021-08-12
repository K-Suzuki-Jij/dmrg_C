#include <stdlib.h>
#include "SML.h"

void FREE_CRS2(CRS1 **Matrix, int row) {
   
   long i;
   
   for (i = 0; i < row; i++) {
      FREE_CRS1(Matrix[i]);
   }
   
   free(Matrix);
   
}
