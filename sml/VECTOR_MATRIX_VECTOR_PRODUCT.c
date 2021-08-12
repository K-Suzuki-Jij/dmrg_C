#include <stdio.h>
#include <stdlib.h>
#include "SML.h"
#include <omp.h>

double VECTOR_MATRIX_VECTOR_PRODUCT(double *V1, CRS1 *M, double *V2, int p_threads) {
   
   long i,j;
   double temp = 0;
   double v1;
   
#pragma omp parallel for private (v1,j) reduction (+:temp) num_threads (p_threads)
   for (i = 0; i < M->row_dim; i++) {
      v1 = V1[i];
      for (j = M->Row[i]; j < M->Row[i+1]; j++) {
         temp = temp + v1*M->Val[j]*V2[M->Col[j]];
      }
   }
   
   return temp;
   
}
