#include <stdlib.h>
#include "SML.h"

void FREE_CRS1(CRS1 *Matrix) {
   
   free(Matrix->Val);
   free(Matrix->Row);
   free(Matrix->Col);
   free(Matrix);
   
}

