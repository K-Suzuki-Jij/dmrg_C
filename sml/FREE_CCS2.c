#include <stdlib.h>
#include "SML.h"

void FREE_CCS2(CCS1 **Matrix, int row) {
   
   long i;
   
   for (i = 0; i < row; i++) {
      FREE_CCS1(Matrix[i]);
   }
   
   free(Matrix);
   
}
