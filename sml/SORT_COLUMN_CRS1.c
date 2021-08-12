#include <stdio.h>
#include <omp.h>
#include "SML.h"

void SORT_COLUMN_CRS1(CRS1 *M, int p_threads) {
   
   long i;
   int dim = M->row_dim;

#pragma omp parallel for num_threads (p_threads)
   for (i = 0; i < dim; i++) {
      if (M->Row[i+1] - M->Row[i] >= 2) {
         QUICK_SORT_INT1_DOUBLE1(M->Col, M->Val, M->Row[i], M->Row[i+1]);
      }
   }
   
}
