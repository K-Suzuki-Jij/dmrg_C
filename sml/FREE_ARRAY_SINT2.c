#include <stdlib.h>

void FREE_ARRAY_SINT2(short int **Matrix, long row) {
   
   long i;
   
   for (i = 0; i < row; i++) {
      free(Matrix[i]);
   }
   
   free(Matrix);
   
}
