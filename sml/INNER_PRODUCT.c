#include <omp.h>

double INNER_PRODUCT(double *V1, double *V2, long dim, int p_threads) {
   
   long i;
   double out = 0;
   
#pragma omp parallel for reduction (+:out) num_threads (p_threads)
   for (i = 0; i < dim; i++) {
      out = out + V1[i]*V2[i];
   }
   
   return out;
   
}
