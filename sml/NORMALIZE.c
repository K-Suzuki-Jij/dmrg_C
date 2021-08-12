#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void NORMALIZE(double *V, long dim, int p_threads) {
   
   long i;
   double norm = 0;
   
#pragma omp parallel for reduction(+:norm) num_threads (p_threads)
   for (i = 0; i < dim; i++) {
      norm = norm + V[i]*V[i];
   }
   
   norm = 1.0/sqrt(norm);
   
#pragma omp parallel for num_threads (p_threads)
   for (i = 0; i < dim; i++) {
      V[i] = norm*V[i];
   }
   
}
