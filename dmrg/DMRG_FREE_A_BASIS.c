#include <stdlib.h>
#include "dmrg.h"

void DMRG_FREE_A_BASIS(DMRG_A_BASIS **A_Basis, int p_threads) {
   
   int thread_num;
   
   for (thread_num = 0; thread_num < p_threads; thread_num++) {
      FREE_ARRAY_SINT1(A_Basis[thread_num]->LL_LLLRRRRL);
      FREE_ARRAY_SINT1(A_Basis[thread_num]->LR_LLLRRRRL);
      FREE_ARRAY_SINT1(A_Basis[thread_num]->RL_LLLRRRRL);
      FREE_ARRAY_SINT1(A_Basis[thread_num]->RR_LLLRRRRL);
      FREE_ARRAY_DOUBLE1(A_Basis[thread_num]->Val_LLLRRRRL);
      FREE_ARRAY_INT1(A_Basis[thread_num]->Inv_LLLRRRRL);
      free(A_Basis[thread_num]);
   }
   
   free(A_Basis);
   
}
