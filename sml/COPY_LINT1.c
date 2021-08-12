#include <omp.h>

void COPY_LINT1(long *V, long *Copyed, long dim, int p_threads) {
   
   long i;
   
#pragma omp parallel for num_threads (p_threads)
   for (i = 0; i < dim; i++) {
      Copyed[i] = V[i];
   }
   
}
