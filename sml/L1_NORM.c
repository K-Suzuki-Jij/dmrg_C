#include <omp.h>
#include <math.h>

double L1_NORM(double *V, long dim, int p_threads) {
  
  long i;
  double norm = 0;
  
  #pragma omp parallel for reduction(+:norm) num_threads (p_threads)
  for (i = 0; i < dim; i++) {
    norm = norm + fabs(V[i]);
  }
  
  return norm;
  
}
