#include <stdio.h>
#include "SML.h"

void PRINT_CCS1(CCS1 *Matrix, char Name[]) {
   
   long i,j;
   
   for (i = 0; i < Matrix->col_dim; i++) {
      for (j = Matrix->Col[i]; j < Matrix->Col[i+1]; j++) {
         printf("%s[%d][%ld]=%.15lf\n", Name, Matrix->Row[j], i, Matrix->Val[j]);
      }
   }
}
