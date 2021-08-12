#include <stdlib.h>

void FREE_ARRAY_INT2(int **Matrix, long row) {
   
   long i;
   
   for (i = 0; i < row; i++) {
      free(Matrix[i]);
   }
   free(Matrix);
}
