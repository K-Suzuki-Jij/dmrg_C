#include "SML.h"

void MATRIX_CONSTAN_MULTIPLICATION_CRS1(CRS1 *M, double c, int p_threads) {
   
   long i,j;
   
#pragma omp parallel for private (j) num_threads (p_threads)
   for (i = 0; i < M->row_dim; i++) {
      for (j = M->Row[i]; j < M->Row[i+1]; j++) {
         M->Val[j] = c*M->Val[j];
      }
   }
   
}

















