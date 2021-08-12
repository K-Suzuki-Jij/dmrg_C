#include <stdio.h>
#include <stdlib.h>
#include "SML.h"

void PRINT_DOUBLE1(double *Vec, long dim, char Name[]) {
   
   long i;
   
   for (i = 0; i < dim; i++) {
      printf("%s[%ld]=%.15lf\n", Name, i, Vec[i] );
   }
   
}
