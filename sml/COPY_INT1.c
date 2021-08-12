#include <omp.h>

void COPY_INT1(int *V, int *Copyed, long dim, int p_threads) {
   
   long i;
   
#pragma omp parallel for num_threads (p_threads)
   for (i = 0; i < dim; i++) {
      Copyed[i] = V[i];
   }
   
}
