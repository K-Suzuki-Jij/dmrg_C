#include <stdio.h>
#include <stdlib.h>
#include "SML.h"
#include <omp.h>

void MATRIX_VECTOR_PRODUCT(CRS1 *M, double *V1, double *Out, int p_threads) {
   
   
   long i,j;
   double temp;
   
#pragma omp parallel for private (j,temp) num_threads (p_threads)
   for (i = 0; i < M->row_dim; i++) {
      temp = 0;
      for (j = M->Row[i]; j < M->Row[i+1]; j++) {
         temp = temp + M->Val[j]*V1[M->Col[j]];
      }
      Out[i] = temp;
   }
   
}
