#include <stdlib.h>
#include "SML.h"

void FREE_LCRS1(LCRS1 *Matrix) {
   
   free(Matrix->Val);
   free(Matrix->Row);
   free(Matrix->Col);
   free(Matrix);
   
}

