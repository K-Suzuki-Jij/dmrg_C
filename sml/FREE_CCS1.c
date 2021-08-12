#include <stdlib.h>
#include "SML.h"

void FREE_CCS1(CCS1 *Matrix) {
   
   free(Matrix->Val);
   free(Matrix->Row);
   free(Matrix->Col);
   free(Matrix);
   
}
