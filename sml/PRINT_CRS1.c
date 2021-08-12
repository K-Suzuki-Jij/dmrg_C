#include <stdio.h>
#include "SML.h"

void PRINT_CRS1(CRS1 *Matrix, char Name[]) {
   
   long i,j;
   
   for (i = 0; i < Matrix->row_dim; i++) {
      for (j = Matrix->Row[i]; j < Matrix->Row[i+1]; j++) {
         printf("%s[%-3ld][%-3d]=%+.15lf\n", Name, i, Matrix->Col[j], Matrix->Val[j]);
      }
   }
}
